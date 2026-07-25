#!/usr/bin/env python3
"""
Seed-robustness experiment for the degeneration onset.

Question: the smallest-inscribed-square collapse appeared near h = 0.5 for
seed 7. Structural, or a seed artifact? (h = 1/2 is Brownian regularity, so
a pinned onset would be interesting.)

Metric per (seed, h): min inscribed-square side / curve diameter.
"Small-square population present" := min_side/diam < 0.1
(the two populations are well separated: big ~0.38-0.46, small <0.08).

Usage:
  python3 -u seed_sweep.py run SEED "0.65,0.55,0.5" seeds.json   # appends
  python3 -u seed_sweep.py summarize seeds.json
"""

import sys
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inscribed_squares import weierstrass_curve, find_squares

SMALL_THRESH = 0.1


def run(seed, hs, out_json):
    rows = []
    if os.path.exists(out_json):
        rows = json.load(open(out_json))
    for h in hs:
        curve = weierstrass_curve(h=h, seed=seed)
        squares = find_squares(curve, n_grid=520, side_floor_rel=0.008, verbose=True)
        sides = np.array([s["side"] for s in squares])
        rows.append(dict(seed=seed, h=h, count=len(squares),
                         min_side=float(sides.min() / curve.diam) if len(sides) else None,
                         med_side=float(np.median(sides) / curve.diam) if len(sides) else None))
        json.dump(rows, open(out_json, "w"))   # checkpoint after every run
    print(f"checkpointed {out_json} ({len(rows)} rows)")


def summarize(json_files, outdir="."):
    rows = []
    for jf in json_files:
        rows.extend(json.load(open(jf)))
    seeds = sorted({r["seed"] for r in rows})
    hs = sorted({r["h"] for r in rows}, reverse=True)
    get = {(r["seed"], r["h"]): r for r in rows}

    print(f"{'seed':>4} | " + " | ".join(f"h={h:<4}" for h in hs) + " | onset h")
    onsets = {}
    M = np.full((len(seeds), len(hs)), np.nan)
    for i, s in enumerate(seeds):
        cells = []
        onset = None
        for j, h in enumerate(hs):          # hs descending: first hit = onset
            r = get.get((s, h))
            if r is None or r["min_side"] is None:
                cells.append("  -   ")
                continue
            M[i, j] = r["min_side"]
            small = r["min_side"] < SMALL_THRESH
            cells.append(f"{r['min_side']:.3f}{'*' if small else ' '}")
            if small and onset is None:
                onset = h
        onsets[s] = onset
        print(f"{s:>4} | " + " | ".join(cells) + f" | {onset}")
    print(f"(* = small-square population present, min_side/diam < {SMALL_THRESH})")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    im = ax1.imshow(np.log10(M), aspect="auto", cmap="viridis_r")
    ax1.set_xticks(range(len(hs)), [str(h) for h in hs])
    ax1.set_yticks(range(len(seeds)), [str(s) for s in seeds])
    ax1.set_xlabel("roughness h (smooth -> rough)")
    ax1.set_ylabel("seed")
    ax1.set_title("log10(min square side / diameter)")
    fig.colorbar(im, ax=ax1)
    for i in range(len(seeds)):
        for j in range(len(hs)):
            if not np.isnan(M[i, j]) and M[i, j] < SMALL_THRESH:
                ax1.plot(j, i, "r*", ms=10)

    for s in seeds:
        hh = [h for h in hs if (s, h) in get and get[(s, h)]["min_side"] is not None]
        vv = [get[(s, h)]["min_side"] for h in hh]
        ax2.semilogy(hh, vv, "o-", label=f"seed {s}", alpha=0.8)
    ax2.axhline(SMALL_THRESH, color="gray", ls=":", lw=1)
    ax2.invert_xaxis()
    ax2.set_xlabel("roughness h")
    ax2.set_ylabel("min square side / diameter")
    ax2.set_title("Degeneration onset per seed")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    path = f"{outdir}/seed_sweep.png"
    fig.savefig(path, dpi=160)
    print(f"wrote {path}")
    print("onsets:", onsets)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "run":
        run(int(sys.argv[2]), [float(x) for x in sys.argv[3].split(",")], sys.argv[4])
    elif mode == "summarize":
        summarize(sys.argv[2].split(","), sys.argv[3] if len(sys.argv) > 3 else ".")
