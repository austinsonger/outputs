#!/usr/bin/env python3
"""
Adversarial search for counterexamples to Tao's Conjecture 6 (blog, Nov 2016;
step 1 of his roadmap to the square peg problem).

Statement. Odd k1,k2,k3; for each i, distinct reals y_{i,1..k_i}; convention
y_{i,0} = y_{i,k_i+1} = -infinity. Axioms:
 (i)   all triple sums y_{1,p}+y_{2,q}+y_{3,r} != 0  (1<=p,q,r<=k_i)
 (ii)  non-crossing: for 0 <= p < q <= k_i, same parity:
         sum_{a in {p,p+1}} sum_{b in {q,q+1}} (-1)^{a+b} sgn(y_{i,a}-y_{i,b}) = 0
 (iii) non-crossing sums: for 0 <= p <= k1, 0 <= q <= k2, 0 <= r <= k3, same parity:
         sum over the 8 corners with alternating sign of sgn(y_{1,a}+y_{2,b}+y_{3,c}) = 0
Conclusion:  F := sum_i sum_p (-1)^(p-1) y_{i,p} < 0.

A counterexample = configuration satisfying (i),(ii),(iii) with F >= 0.

Search design:
 - Axiom (ii) constrains only the ORDER pattern of each sequence (a meander
   permutation). Enumerate all valid order patterns per k by brute force once.
 - Annealing: state = (pattern_1, pattern_2, pattern_3, values sorted into
   patterns). Moves perturb values or swap patterns. Axiom (iii) enforced as
   penalty -> hard filter for records. Objective: maximize F.

Usage:
  python3 -u conj6_search.py patterns 7          # count meander patterns for k
  python3 -u conj6_search.py search 5 SECONDS OUT.json
  python3 -u conj6_search.py search 7 SECONDS OUT.json
"""

import sys
import json
import time
import itertools
import numpy as np

NEG_INF = -1e15


def seq_with_sentinels(y):
    return np.concatenate([[NEG_INF], y, [NEG_INF]])


def check_ii(y):
    """Axiom (ii) for one sequence y (finite values, length k)."""
    z = seq_with_sentinels(y)
    k = len(y)
    for p in range(0, k):
        for q in range(p + 2, k + 1, 2):
            s = 0
            for a in (p, p + 1):
                for b in (q, q + 1):
                    d = z[a] - z[b]
                    s += (-1) ** (a + b) * (0 if d == 0 else (1 if d > 0 else -1))
            if s != 0:
                return False
    return True


def enumerate_patterns(k):
    """All permutation patterns of 1..k passing axiom (ii)."""
    pats = []
    for perm in itertools.permutations(range(1, k + 1)):
        if check_ii(np.array(perm, dtype=float)):
            pats.append(perm)
    return pats


def check_iii(y1, y2, y3):
    z1, z2, z3 = map(seq_with_sentinels, (y1, y2, y3))
    k1, k2, k3 = len(y1), len(y2), len(y3)
    for p in range(0, k1 + 1):
        for q in range(p % 2, k2 + 1, 2):
            for r in range(p % 2, k3 + 1, 2):
                s = 0
                for a in (p, p + 1):
                    for b in (q, q + 1):
                        for c in (r, r + 1):
                            v = z1[a] + z2[b] + z3[c]
                            s += (-1) ** (a + b + c) * (1 if v > 0 else -1)
                if s != 0:
                    return False
    return True


def check_i(y1, y2, y3, eps=1e-9):
    S = y1[:, None, None] + y2[None, :, None] + y3[None, None, :]
    return bool(np.min(np.abs(S)) > eps)


def F_obj(y1, y2, y3):
    return sum(((-1) ** np.arange(len(y))) @ y for y in (y1, y2, y3))


def iii_violations(y1, y2, y3):
    z1, z2, z3 = map(seq_with_sentinels, (y1, y2, y3))
    k1, k2, k3 = len(y1), len(y2), len(y3)
    n = 0
    for p in range(0, k1 + 1):
        for q in range(p % 2, k2 + 1, 2):
            for r in range(p % 2, k3 + 1, 2):
                s = 0
                for a in (p, p + 1):
                    for b in (q, q + 1):
                        for c in (r, r + 1):
                            v = z1[a] + z2[b] + z3[c]
                            s += (-1) ** (a + b + c) * (1 if v > 0 else -1)
                n += (s != 0)
    return n


def values_for_pattern(pattern, rng, scale=1.0):
    """Random distinct values arranged so the sequence has the given order pattern."""
    vals = np.sort(rng.uniform(-scale, scale, len(pattern)))
    y = np.empty(len(pattern))
    for pos, rank in enumerate(pattern):
        y[pos] = vals[rank - 1]
    return y


def cell_lp_max(ys, pats, eps=1e-4, box=1.0):
    """
    Exact supremum of F within the sign-cell of configuration ys.
    The cell = {orderings fixed per pattern} x {sign of every triple sum fixed}.
    Keeping every individual sgn unchanged preserves axioms (i)-(iii), so any
    point in the cell is a valid configuration. Cells are open cones: if the
    LP finds F > 0 anywhere in the cell, that is a genuine counterexample.
    Returns (F_max, x_opt) with margins eps and box bound [-box, box].
    """
    from scipy.optimize import linprog
    k1, k2, k3 = (len(y) for y in ys)
    n = k1 + k2 + k3
    off = [0, k1, k1 + k2]

    A_ub, b_ub = [], []

    def add(le_row, rhs):        # le_row . x <= rhs
        A_ub.append(le_row)
        b_ub.append(rhs)

    # ordering constraints: consecutive ranks per pattern, gap >= eps
    for i, (y, pat) in enumerate(zip(ys, pats)):
        order = np.argsort(np.array(pat))          # positions sorted by rank
        for a, b in zip(order[:-1], order[1:]):    # y[a] < y[b]
            row = np.zeros(n)
            row[off[i] + a] = 1.0
            row[off[i] + b] = -1.0
            add(row, -eps)
    # triple-sum sign constraints: s * (y1p + y2q + y3r) >= eps
    y1, y2, y3 = ys
    for p in range(k1):
        for q in range(k2):
            for r in range(k3):
                s = 1.0 if (y1[p] + y2[q] + y3[r]) > 0 else -1.0
                row = np.zeros(n)
                row[off[0] + p] = -s
                row[off[1] + q] = -s
                row[off[2] + r] = -s
                add(row, -eps)
    # objective: maximize F = sum_i sum_p (-1)^p_index y  -> minimize -F
    c = np.concatenate([-((-1.0) ** np.arange(kk)) for kk in (k1, k2, k3)])
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  bounds=[(-box, box)] * n, method="highs")
    if not res.success:
        return None, None
    return -res.fun, res.x


def search(k, seconds, out_json, seed=0):
    rng = np.random.default_rng(seed)
    print(f"enumerating axiom-(ii) patterns for k={k} ...")
    t0 = time.time()
    pats = enumerate_patterns(k)
    print(f"  {len(pats)} patterns ({time.time()-t0:.1f}s)")

    best = dict(F=-np.inf, ys=None, pats=None)
    t_end = time.time() + seconds
    n_restart = n_valid = 0
    seen_cells = set()
    while time.time() < t_end:
        n_restart += 1
        pat = [pats[rng.integers(len(pats))] for _ in range(3)]
        ys = [values_for_pattern(p, rng) for p in pat]
        # local annealing on values, pattern fixed
        T = 0.3
        cur_pen = iii_violations(*ys)
        cur_F = F_obj(*ys)
        for it in range(3500):
            T *= 0.9988
            i = rng.integers(3)
            j = rng.integers(k)
            trial = [y.copy() for y in ys]
            trial[i][j] += rng.normal(0, 0.15 * T + 1e-3)
            if not check_ii(trial[i]):
                continue
            pen = iii_violations(*trial)
            Ft = F_obj(*trial)
            # lexicographic-ish annealing: first drive violations down, then F up
            dE = (pen - cur_pen) * 10.0 - (Ft - cur_F)
            if dE < 0 or rng.random() < np.exp(-dE / max(T, 1e-6)):
                ys, cur_pen, cur_F = trial, pen, Ft
            if cur_pen == 0 and check_i(*ys):
                key = (tuple(map(tuple, pat)),
                       tuple(np.sign(ys[0][:, None, None] + ys[1][None, :, None]
                                     + ys[2][None, None, :]).astype(int).ravel()))
                if key in seen_cells:
                    break                      # cell already settled by LP
                seen_cells.add(key)
                n_valid += 1
                Fmax, x = cell_lp_max(ys, pat)
                if Fmax is not None and Fmax > best["F"]:
                    k123 = [len(y) for y in ys]
                    xs = np.split(x, np.cumsum(k123)[:-1])
                    best = dict(F=float(Fmax), ys=[v.tolist() for v in xs],
                                pats=[list(p) for p in pat])
                    if Fmax > 0:
                        print("  *** F > 0 cell found: COUNTEREXAMPLE CANDIDATE ***")
                break                          # restart for cell diversity
    print(f"k={k}: {n_restart} restarts, distinct valid cells settled by LP: "
          f"{len(seen_cells)}")
    if best["ys"] is not None:
        print(f"  best cell-supremum F = {best['F']:.6f} "
              f"(counterexample requires F >= 0)")
        json.dump(best, open(out_json, "w"))
        print(f"  wrote {out_json}")
    else:
        print("  no fully valid configuration found in budget")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "patterns":
        k = int(sys.argv[2])
        pats = enumerate_patterns(k)
        print(f"k={k}: {len(pats)} axiom-(ii) patterns")
    elif mode == "search":
        k = int(sys.argv[2])
        secs = float(sys.argv[3])
        out = sys.argv[4]
        seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
        search(k, secs, out, seed)
