#!/usr/bin/env python3
"""
Margin-law follow-ups.

Part 1 (mixed sizes): generalized MILP for arbitrary (k1,k2,k3) to test the
conjectured law  sup F = -((k1+k2+k3-1)/2) * EPS.
  Predictions: (5,5,7) -> -8e-5,  (5,7,9) -> -10e-5.

Part 2 (certificates): for an optimal k=5 cell, fix its sign-cell and solve
the within-cell LP, then extract the DUAL multipliers. The dual is a
nonnegative combination of the cell's margin constraints that proves
F <= -c*EPS for that cell: a machine-written certificate. Compare certificate
structure across pattern triples to test for uniformity.

Usage:
  python3 -u dual_cert.py mixed 5 5 7
  python3 -u dual_cert.py dual IDX          # cell from ex_a.json record IDX
"""

import sys
import json
import numpy as np
from scipy.optimize import milp, linprog, LinearConstraint, Bounds

from conj6_search import enumerate_patterns

EPS = 1e-5
BIG_M = 6.0


def solve_mixed(P, time_limit=30.0):
    """Max F over configs with pattern triple P = (pat1, pat2, pat3),
    sequence lengths implied by pattern lengths."""
    ks = [len(p) for p in P]
    ny = sum(ks)
    off = [0, ks[0], ks[0] + ks[1]]
    ns = ks[0] * ks[1] * ks[2]
    n = ny + ns

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
    return (float(-res.fun) if res.success else None), res.status


def mixed(k1, k2, k3):
    pats = {k: enumerate_patterns(k) for k in {k1, k2, k3}}
    pred = -((k1 + k2 + k3 - 1) / 2) * EPS
    print(f"(k1,k2,k3)=({k1},{k2},{k3}), prediction sup F = {pred:.6e}")
    rng = np.random.default_rng(0)
    for trial in range(3):
        P = [pats[k][rng.integers(len(pats[k]))] for k in (k1, k2, k3)]
        F, status = solve_mixed(P)
        print(f"  trial {trial}: F = {F}   (status {status})")


# ---------------------------------------------------------------------------
# Part 2: dual certificates
# ---------------------------------------------------------------------------

def cell_lp_with_duals(ys, pats, eps=EPS, box=1.0):
    """Within-cell LP (signs fixed by ys) + dual multipliers.
    Constraint list (A_ub x <= b_ub):
      type 'ord':  ordering gap  (y_lo - y_hi <= -eps)
      type 'sum':  triple-sum sign  (-s*(sum) <= -eps)
    Returns dict with F_max, x, duals (list of (type, meta, lambda))."""
    ks = [len(y) for y in ys]
    n = sum(ks)
    off = [0, ks[0], ks[0] + ks[1]]
    A, b, meta = [], [], []
    for i, (y, pat) in enumerate(zip(ys, pats)):
        order = np.argsort(np.array(pat))
        for a_, b_ in zip(order[:-1], order[1:]):
            row = np.zeros(n)
            row[off[i] + a_] = 1.0
            row[off[i] + b_] = -1.0
            A.append(row)
            b.append(-eps)
            meta.append(("ord", (i, int(a_), int(b_))))
    for p in range(ks[0]):
        for q in range(ks[1]):
            for r in range(ks[2]):
                s = 1.0 if (ys[0][p] + ys[1][q] + ys[2][r]) > 0 else -1.0
                row = np.zeros(n)
                row[off[0] + p] = -s
                row[off[1] + q] = -s
                row[off[2] + r] = -s
                A.append(row)
                b.append(-eps)
                meta.append(("sum", (p, q, r, int(s))))
    c = np.concatenate([-((-1.0) ** np.arange(k)) for k in ks])
    res = linprog(c, A_ub=np.array(A), b_ub=np.array(b),
                  bounds=[(-box, box)] * n, method="highs")
    if not res.success:
        return None
    lam = -res.ineqlin.marginals          # nonnegative duals
    duals = [(meta[j][0], meta[j][1], float(lam[j]))
             for j in range(len(meta)) if abs(lam[j]) > 1e-9]
    # bound duals (nonzero means the box is active: certificate artifact)
    bnd = [(int(j), float(res.lower.marginals[j]), float(res.upper.marginals[j]))
           for j in range(n)
           if abs(res.lower.marginals[j]) > 1e-9 or abs(res.upper.marginals[j]) > 1e-9]
    return dict(F=float(-res.fun), x=res.x, duals=duals, bound_duals=bnd)


def dual_report(idx):
    import exhaustive_k5 as ex
    pats_all, triples = ex.triple_list()
    recs = json.load(open("ex_a.json")) + json.load(open("ex_b.json")) \
        + json.load(open("ex_c.json"))
    rec = recs[idx]
    tri = rec["tri"]
    P = [pats_all[i] for i in tri]
    y = np.array(rec["y"])
    ys = [y[0:5], y[5:10], y[10:15]]
    out = cell_lp_with_duals(ys, P)
    print(f"record {idx}, tri={tri}, cell LP F_max = {out['F']:.6e} "
          f"(expected -7e-5)")
    n_ord = sum(1 for d in out["duals"] if d[0] == "ord")
    n_sum = sum(1 for d in out["duals"] if d[0] == "sum")
    tot = sum(d[2] for d in out["duals"])
    print(f"  active dual constraints: {n_ord} ordering + {n_sum} triple-sum, "
          f"sum of multipliers = {tot:.4f}")
    for d in sorted(out["duals"], key=lambda d: -d[2]):
        print(f"    {d[0]} {d[1]}  lambda={d[2]:.4f}")
    if out["bound_duals"]:
        print(f"  WARNING box bounds active: {out['bound_duals']}")
    else:
        print("  box bounds inactive (certificate is purely cone-internal)")
    return out


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "mixed":
        mixed(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    elif mode == "dual":
        dual_report(int(sys.argv[2]))
