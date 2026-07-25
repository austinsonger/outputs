#!/usr/bin/env python3
"""
Solve the MILP for specific unbalanced pattern triples that had no
aligned-edge-only ladder cell, then test whether their GLOBAL optimum has one.
"""

import sys
import json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

from conj6_search import enumerate_patterns
from cert_test_higher import has_aligned_certificate

EPS = 1e-5
BIG_M = 6.0


def solve_triple(P, time_limit=30.0):
    """Maximize F over configs with pattern triple P; return (F, ys) or (None,None)."""
    ks = [len(p) for p in P]
    ny = sum(ks)
    ns = ks[0] * ks[1] * ks[2]
    n = ny + ns
    off = [0, ks[0], ks[0] + ks[1]]

    def yvar(i, p):
        return off[i] + p

    def svar(p, q, r):
        return ny + (p * ks[1] + q) * ks[2] + r

    cons = []
    for i, pat in enumerate(P):
        order = np.argsort(np.array(pat))
        for a, b in zip(order[:-1], order[1:]):
            row = np.zeros(n)
            row[yvar(i, a)] = -1.0
            row[yvar(i, b)] = 1.0
            cons.append(LinearConstraint(row, EPS, np.inf))
    for p in range(ks[0]):
        for q in range(ks[1]):
            for r in range(ks[2]):
                row = np.zeros(n)
                row[yvar(0, p)] = row[yvar(1, q)] = row[yvar(2, r)] = 1.0
                row[svar(p, q, r)] = -BIG_M
                cons.append(LinearConstraint(row, -np.inf, -EPS))
                cons.append(LinearConstraint(row, EPS - BIG_M, np.inf))
    for par in (0, 1):
        for p in range(par, ks[0] + 1, 2):
            for q in range(par, ks[1] + 1, 2):
                for r in range(par, ks[2] + 1, 2):
                    row = np.zeros(n)
                    const = 0.0
                    for a in (p, p + 1):
                        for b in (q, q + 1):
                            for c in (r, r + 1):
                                sign = (-1.0) ** (a + b + c)
                                if (1 <= a <= ks[0] and 1 <= b <= ks[1]
                                        and 1 <= c <= ks[2]):
                                    row[svar(a - 1, b - 1, c - 1)] += 2 * sign
                                    const += -sign
                                else:
                                    const += -sign
                    cons.append(LinearConstraint(row, -const, -const))
    c = np.zeros(n)
    for i in range(3):
        for p in range(ks[i]):
            c[yvar(i, p)] = -((-1.0) ** p)
    integrality = np.zeros(n)
    integrality[ny:] = 1
    bounds = Bounds(np.concatenate([-np.ones(ny), np.zeros(ns)]), np.ones(n))
    res = milp(c, constraints=cons, integrality=integrality, bounds=bounds,
               options=dict(time_limit=time_limit))
    if not res.success:
        return None, None
    y = res.x[:ny]
    ys = [y[off[i]:off[i] + ks[i]] for i in range(3)]
    return float(-res.fun), ys


def main(no_ok_file):
    results = json.load(open(no_ok_file))
    print(f"testing {len(results)} triples")
    all_ok = True
    for rec in results:
        ks = rec["ks"]
        pats = [tuple(p) for p in rec["pats"]]
        F, ys = solve_triple(pats)
        if F is None:
            print(f"  k={ks} {pats}: MILP failed")
            all_ok = False
            continue
        ok, rows, target = has_aligned_certificate(ys, time_limit=10.0)
        print(f"  k={ks} {pats}: F={F:.6e}, aligned-cert ok={ok}, rows={rows}/{target}")
        if not ok:
            all_ok = False
    print("all global optima have aligned-edge cert:", all_ok)


if __name__ == "__main__":
    # build list of no-ok triples from aligned_higher_results.json
    res = json.load(open("aligned_higher_results.json"))
    from collections import defaultdict
    groups = defaultdict(list)
    for r in res:
        key = (tuple(r['ks']), tuple(tuple(p) for p in sorted(r['pats'])))
        groups[key].append(r)
    no_ok = []
    for key, vals in groups.items():
        if not any(r['ok'] for r in vals):
            ks, pats = key
            no_ok.append(dict(ks=list(ks), pats=[list(p) for p in pats]))
    json.dump(no_ok, open("no_ok_triples.json", "w"))
    print(f"found {len(no_ok)} triples with no ok ladder cell")
    main("no_ok_triples.json")
