#!/usr/bin/env python3
"""
Exhaustive verification of Tao's Conjecture 6 base case k=(5,5,5) via MILP.

For each triple of axiom-(ii) order patterns (8^3 = 512; 120 up to permuting
the three sequences, which is a symmetry of the axioms and of F):

  variables:  y (15 continuous in [-1,1]),  s_pqr (125 binaries = sign bits
              of the triple sums y1p + y2q + y3r)
  constraints:
    - orderings: consecutive rank gaps >= EPS       (axiom ii via pattern)
    - big-M linking: s=1 <=> sum >= EPS, s=0 <=> sum <= -EPS   (axiom i)
    - axiom (iii): for every same-parity window (p,q,r), the alternating sum
      of signs over the 8 corners equals 0 (linear equality in s; corners
      involving sentinel indices 0 or 6 contribute fixed sign -1)
  objective:  maximize F = sum_i (y_i1 - y_i2 + y_i3 - y_i4 + y_i5)

If the optimum is < 0 for every pattern triple, no counterexample to
Conjecture 6 exists at k=(5,5,5) with relative margins >= EPS.
(Full certificate would repeat the argument in rational arithmetic; cells are
open cones, so any counterexample admits a margin-normalized representative.)

Usage:
  python3 -u exhaustive_k5.py run START END OUT.json   # triples index range
  python3 -u exhaustive_k5.py combine OUT1.json,OUT2.json,...
"""

import sys
import json
import time
import itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

from conj6_search import enumerate_patterns

EPS = 1e-5
BIG_M = 4.0
K = 5


def triple_list():
    pats = enumerate_patterns(K)
    # symmetry reduction: sequences are interchangeable -> multisets of patterns
    seen, triples = set(), []
    for t in itertools.product(range(len(pats)), repeat=3):
        key = tuple(sorted(t))
        if key not in seen:
            seen.add(key)
            triples.append(t)
    return pats, triples


def solve_triple(pats, tri, time_limit=20.0):
    """Max F over configurations with the given pattern triple. Returns dict."""
    P = [pats[i] for i in tri]
    ny, ns = 3 * K, K ** 3
    n = ny + ns

    def yvar(i, p):          # i in 0..2, p in 0..K-1  (value y_{i,p+1})
        return i * K + p

    def svar(p, q, r):       # 0-based finite indices
        return ny + (p * K + q) * K + r

    cons = []

    # ordering constraints per pattern: consecutive ranks gap >= EPS
    for i, pat in enumerate(P):
        order = np.argsort(np.array(pat))
        for a, b in zip(order[:-1], order[1:]):    # y[a] < y[b]
            row = np.zeros(n)
            row[yvar(i, a)] = -1.0
            row[yvar(i, b)] = 1.0
            cons.append(LinearConstraint(row, EPS, np.inf))

    # big-M sign linking
    for p in range(K):
        for q in range(K):
            for r in range(K):
                rows = np.zeros(n)
                rows[yvar(0, p)] = rows[yvar(1, q)] = rows[yvar(2, r)] = 1.0
                rows[svar(p, q, r)] = -BIG_M
                # sum - M s <= -EPS  (s=0 -> sum <= -EPS)
                cons.append(LinearConstraint(rows, -np.inf, -EPS + 0.0))
                rows2 = np.zeros(n)
                rows2[yvar(0, p)] = rows2[yvar(1, q)] = rows2[yvar(2, r)] = 1.0
                rows2[svar(p, q, r)] = -BIG_M
                # sum >= EPS - M(1-s)  ->  sum - M s >= EPS - M
                cons.append(LinearConstraint(rows2, EPS - BIG_M, np.inf))

    # axiom (iii): same-parity windows, alternating corner-sign sum = 0
    for par in (0, 1):
        Ps = range(par, K + 1, 2)
        for p in Ps:
            for q in Ps:
                for r in Ps:
                    row = np.zeros(n)
                    const = 0.0
                    for a in (p, p + 1):
                        for b in (q, q + 1):
                            for c in (r, r + 1):
                                sign = (-1.0) ** (a + b + c)
                                if (1 <= a <= K and 1 <= b <= K
                                        and 1 <= c <= K):
                                    # sgn = 2 s - 1
                                    row[svar(a - 1, b - 1, c - 1)] += 2 * sign
                                    const += -sign
                                else:
                                    const += -sign        # sentinel: sgn = -1
                    cons.append(LinearConstraint(row, -const, -const))

    # objective: maximize F -> minimize -F
    c = np.zeros(n)
    for i in range(3):
        for p in range(K):
            c[yvar(i, p)] = -((-1.0) ** p)

    integrality = np.zeros(n)
    integrality[ny:] = 1
    bounds = Bounds(np.concatenate([-np.ones(ny), np.zeros(ns)]),
                    np.ones(n))
    res = milp(c, constraints=cons, integrality=integrality, bounds=bounds,
               options=dict(time_limit=time_limit))
    if res.status == 2:            # infeasible: no valid config for this triple
        return dict(tri=list(tri), status="infeasible", F=None)
    if not res.success:
        return dict(tri=list(tri), status=f"status{res.status}", F=None)
    return dict(tri=list(tri), status="ok", F=float(-res.fun),
                y=res.x[:ny].tolist())


def run(start, end, out_json):
    pats, triples = triple_list()
    print(f"{len(pats)} patterns, {len(triples)} pattern triples up to symmetry")
    results = []
    t0 = time.time()
    for idx in range(start, min(end, len(triples))):
        r = solve_triple(pats, triples[idx])
        results.append(r)
        f = "inf" if r["F"] is None else f"{r['F']:.6f}"
        print(f"  [{idx}] tri={r['tri']} {r['status']} F={f} "
              f"({time.time()-t0:.0f}s)")
        if r["F"] is not None and r["F"] >= 0:
            print("  *** F >= 0: COUNTEREXAMPLE CANDIDATE, verify exactly ***")
    json.dump(results, open(out_json, "w"))
    print(f"wrote {out_json}")


def combine(files):
    results = []
    for f in files:
        results.extend(json.load(open(f)))
    ok = [r for r in results if r["status"] == "ok"]
    bad = [r for r in results if r["status"].startswith("status")]
    infeas = [r for r in results if r["status"] == "infeasible"]
    Fs = [r["F"] for r in ok]
    print(f"{len(results)} triples: {len(ok)} solved, {len(infeas)} infeasible,"
          f" {len(bad)} solver issues")
    if Fs:
        print(f"max F over all solved triples = {max(Fs):.6f}")
        print(f"min F = {min(Fs):.6f}")
        print("VERDICT:", "NO counterexample at k=(5,5,5) (margins >= EPS)"
              if max(Fs) < 0 else "COUNTEREXAMPLE CANDIDATE FOUND")
    if bad:
        print("retry needed for:", [r["tri"] for r in bad])


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "run":
        run(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    elif mode == "combine":
        combine(sys.argv[2].split(","))
