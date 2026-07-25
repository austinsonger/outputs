#!/usr/bin/env python3
"""
Arc-inflation ladder for Tao's Conjecture 6.

Idea: random search cannot reach axiom-(iii)-valid configurations at k >= 7,
but valid configurations can be GROWN. Inserting an adjacent nested pair
(u, u+delta) into a valid sequence adds a small arc to the meander (axiom ii
preserved for suitable positions/values), and if u sits close to an existing
value the new triple-sum signs copy the neighbor's, usually preserving (iii).
All candidates are verified against the exact axioms; survivors' sign-cells
are settled exactly by LP. Ladder: k=5 -> 7 -> 9 -> ...

Usage:
  python3 -u inflate.py ladder IN1.json,IN2.json OUT.json TIME_BUDGET_S
      IN*.json: files with {"ys": [[...],[...],[...]]} valid configs
      OUT.json: {"cells": [...records...], "best": record}
  python3 -u inflate.py walk IN.json OUT.json TIME_BUDGET_S   # neighbor-cell walk
"""

import sys
import json
import time
import numpy as np

from conj6_search import (check_i, check_ii, check_iii, F_obj, cell_lp_max,
                          iii_violations)


def valid(ys):
    return (all(check_ii(y) for y in ys) and check_i(*ys) and check_iii(*ys))


def pattern_of(y):
    return tuple(int(r) for r in (np.argsort(np.argsort(y)) + 1))


def cell_key(ys):
    S = np.sign(ys[0][:, None, None] + ys[1][None, :, None]
                + ys[2][None, None, :]).astype(int)
    return (tuple(pattern_of(y) for y in ys),
            tuple(S.ravel().tolist()))


def candidate_insertions(ys, rng):
    """Yield inflated configs: insert adjacent pair into the SHORTEST sequence
    (balanced growth: (5,5,5) -> (7,5,5) -> (7,7,5) -> (7,7,7) -> ...)."""
    ks = [len(y) for y in ys]
    for i in [int(np.argmin(ks))]:
        y = ys[i]
        k = len(y)
        for p in range(k + 1):                      # insertion gap after pos p
            bases = []
            if p > 0:
                bases.append(y[p - 1])              # near left neighbor
            if p < k:
                bases.append(y[p])                  # near right neighbor
            bases.append(y.min() - 0.05)            # new bottom arc
            for w in bases:
                for d1, d2 in ((1e-3, 2e-3), (2e-3, 1e-3)):
                    u, v = w + d1, w + d2
                    y_new = np.concatenate([y[:p], [u, v], y[p:]])
                    trial = [ys[0].copy(), ys[1].copy(), ys[2].copy()]
                    trial[i] = y_new
                    if valid(trial):
                        yield trial


def settle(ys, cells, best):
    key = cell_key(ys)
    if key in cells:
        return best
    pats = [pattern_of(y) for y in ys]
    Fmax, x = cell_lp_max(ys, pats)
    if Fmax is None:
        return best
    ks = [len(y) for y in ys]
    xs = [v.tolist() for v in np.split(x, np.cumsum(ks)[:-1])]
    rec = dict(F=float(Fmax), ys=xs, k=ks)
    cells[key] = rec
    if Fmax > best["F"]:
        best.update(rec)
        if Fmax > 0:
            print("  *** F > 0: COUNTEREXAMPLE CANDIDATE ***")
    return best


def neighbor_walk(ys, cells, best, rng, steps=400, sigma=0.02):
    """Random walk within the valid region; settle every new cell touched."""
    cur = [y.copy() for y in ys]
    for _ in range(steps):
        i = rng.integers(3)
        j = rng.integers(len(cur[i]))
        trial = [y.copy() for y in cur]
        trial[i][j] += rng.normal(0, sigma)
        if valid(trial):
            cur = trial
            best = settle(cur, cells, best)
    return best


def ladder(in_files, out_json, budget_s):
    rng = np.random.default_rng(0)
    t_end = time.time() + budget_s
    seeds = []
    for f in in_files:
        d = json.load(open(f))
        ys = [np.array(v, dtype=float) for v in d["ys"]]
        if valid(ys):
            seeds.append(ys)
        else:
            print(f"  note: config in {f} not valid as stored (LP boundary "
                  f"rounding), nudging values apart")
            # nudge: spread coincident values slightly, retry
            for i in range(3):
                y = ys[i]
                for a in range(len(y)):
                    y[a] += 1e-5 * a
            if valid(ys):
                seeds.append(ys)
    print(f"{len(seeds)} valid seed configs loaded")

    cells, best = {}, dict(F=-np.inf)
    frontier = seeds
    level = 0
    while frontier and time.time() < t_end:
        k_now = [len(y) for y in frontier[0]]
        print(f"level {level}: {len(frontier)} configs at k={k_now}")
        next_frontier = []
        for ys in frontier:
            if time.time() > t_end:
                break
            best = settle(ys, cells, best)
            n_inf = 0
            for infl in candidate_insertions(ys, rng):
                best = settle(infl, cells, best)
                next_frontier.append(infl)
                n_inf += 1
                if n_inf >= 12 or time.time() > t_end:   # cap fan-out
                    break
        # dedupe frontier by cell, keep at most 10 configs for next level
        seen, uniq = set(), []
        for ys in next_frontier:
            key = cell_key(ys)
            if key not in seen:
                seen.add(key)
                uniq.append(ys)
        frontier = uniq[:10]
        level += 1

    by_k = {}
    for key, rec in cells.items():
        kk = tuple(rec["k"])
        by_k.setdefault(kk, []).append(rec["F"])
    for kk in sorted(by_k):
        Fs = by_k[kk]
        print(f"k={kk}: {len(Fs)} cells settled, sup F = {max(Fs):.6f}")
    json.dump(dict(best=best,
                   summary={str(kk): dict(cells=len(v), supF=max(v))
                            for kk, v in by_k.items()},
                   cells=list(cells.values())),
              open(out_json, "w"))
    print(f"wrote {out_json}; global best F = {best['F']:.6f}")


def walk_mode(in_file, out_json, budget_s):
    rng = np.random.default_rng(1)
    d = json.load(open(in_file))
    ys = [np.array(v, dtype=float) for v in d["best"]["ys"]]
    if not valid(ys):
        for i in range(3):
            for a in range(len(ys[i])):
                ys[i][a] += 1e-5 * a
    assert valid(ys), "walk seed invalid"
    cells, best = {}, dict(F=-np.inf)
    t_end = time.time() + budget_s
    while time.time() < t_end:
        best = neighbor_walk(ys, cells, best, rng, steps=150,
                             sigma=10 ** rng.uniform(-2.3, -0.7))
    print(f"walk: {len(cells)} cells settled, best F = {best['F']:.6f}")
    json.dump(dict(best=best, n_cells=len(cells)), open(out_json, "w"))


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "ladder":
        ladder(sys.argv[2].split(","), sys.argv[3], float(sys.argv[4]))
    elif mode == "walk":
        walk_mode(sys.argv[2], sys.argv[3], float(sys.argv[4]))
