#!/usr/bin/env python3
"""
Continuation tracker: follow each inscribed-square branch as roughness h varies.

At each h step, every live branch's (t1,t2) is re-polished (warm start) on the
new curve. A branch dies when the polished residual no longer certifies a
square (after perturbation retries). Deaths in the down-pass = annihilations
as roughness grows; deaths in the up-pass = branches that only exist at high
roughness (births, seen in reverse).

Usage:
  python3 -u continuation.py seedpass SEED H_START H_END N_STEPS state.json  # init + steps
  python3 -u continuation.py steps state.json N                              # continue N steps
  python3 -u continuation.py plot state_down.json,state_up.json out.png
"""

import sys
import json
import os
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inscribed_squares import weierstrass_curve, find_squares, square_corners

SEED_DEFAULT = 3


def polish(curve, x0):
    def obj(x):
        c = square_corners(curve, x[0], x[1])
        d = curve.dist_exact(c[[1, 3]])
        return float(d @ d)
    res = minimize(obj, x0=np.asarray(x0), method="Nelder-Mead",
                   options=dict(xatol=1e-12, fatol=1e-20, maxiter=300))
    cs = square_corners(curve, res.x[0], res.x[1])
    side = float(np.linalg.norm(cs[1] - cs[0]))
    err = float(curve.dist_exact(cs[[1, 3]]).max())
    return res.x, side, err


def certify(curve, side, err):
    return err <= min(1e-7 * curve.diam, 0.01 * side)


def step_state(state):
    """Advance the state by one h step, re-polishing every live branch."""
    hs = state["h_schedule"]
    idx = state["h_index"] + 1
    if idx >= len(hs):
        return False
    h = hs[idx]
    curve = weierstrass_curve(h=h, seed=state["seed"])
    rng = np.random.default_rng(1234 + idx)
    n_dead = 0
    for br in state["branches"]:
        if not br["alive"]:
            continue
        x0 = br["hist"][-1][1:3]
        x, side, err = polish(curve, x0)
        if not certify(curve, side, err):
            for _ in range(3):                        # perturbation retries
                xr, sr, er = polish(curve, np.asarray(x0) + rng.normal(0, 2e-3, 2))
                if certify(curve, sr, er):
                    x, side, err = xr, sr, er
                    break
            else:
                br["alive"] = False
                br["died_at"] = h
                n_dead += 1
                continue
        br["hist"].append([h, float(x[0]), float(x[1]),
                           side / curve.diam, err])
    state["h_index"] = idx
    alive = sum(b["alive"] for b in state["branches"])
    print(f"  h={h:.3f}: {alive} alive, {n_dead} died")
    return True


def seedpass(seed, h_start, h_end, n_steps, out_json):
    curve = weierstrass_curve(h=h_start, seed=seed)
    squares = find_squares(curve, n_grid=520, side_floor_rel=0.006, verbose=True)
    hs = np.linspace(h_start, h_end, n_steps + 1).tolist()
    state = dict(seed=seed, h_schedule=hs, h_index=0,
                 branches=[dict(id=k, alive=True, died_at=None,
                                hist=[[h_start, s["x"][0], s["x"][1],
                                       s["side"] / curve.diam, s["err"]]])
                           for k, s in enumerate(squares)])
    print(f"seeded {len(squares)} branches at h={h_start}")
    json.dump(state, open(out_json, "w"))


def steps(state_json, n):
    state = json.load(open(state_json))
    for _ in range(n):
        if not step_state(state):
            print("schedule complete")
            break
        json.dump(state, open(state_json, "w"))   # checkpoint every step
    done = state["h_index"] >= len(state["h_schedule"]) - 1
    print(f"progress: {state['h_index']}/{len(state['h_schedule'])-1} steps"
          + (" (DONE)" if done else ""))


def plot(state_jsons, out_png):
    fig, ax = plt.subplots(figsize=(11, 6))
    styles = ["-", "--"]
    labels = ["down-pass (0.70 -> 0.40)", "up-pass (0.40 -> 0.70)"]
    for si, sj in enumerate(state_jsons):
        state = json.load(open(sj))
        cmap = plt.get_cmap("turbo")
        nb = len(state["branches"])
        for br in state["branches"]:
            hist = np.array(br["hist"])
            col = cmap(0.05 + 0.9 * br["id"] / max(1, nb - 1))
            ax.semilogy(hist[:, 0], hist[:, 3], styles[si], color=col,
                        lw=1.1, alpha=0.75)
            if br["died_at"] is not None:
                ax.semilogy(hist[-1, 0], hist[-1, 3], "x", color=col, ms=8,
                            mew=2)
        n_died = sum(b["died_at"] is not None for b in state["branches"])
        print(f"{sj}: {len(state['branches'])} branches, {n_died} died")
    ax.invert_xaxis()
    ax.set_xlabel("roughness h  (smooth -> rough)")
    ax.set_ylabel("square side / curve diameter (log)")
    ax.set_title("Square branches under continuation in h "
                 "(x = branch death; solid = down-pass, dashed = up-pass)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "seedpass":
        seedpass(int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]),
                 int(sys.argv[5]), sys.argv[6])
    elif mode == "steps":
        steps(sys.argv[2], int(sys.argv[3]))
    elif mode == "plot":
        plot(sys.argv[2].split(","), sys.argv[3])
