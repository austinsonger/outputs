#!/usr/bin/env python3
"""
Test the aligned-edge-only certificate conjecture on cells from inflate.py.

Usage:
  python3 -u cert_test_higher.py test CELLS.json [MAX_K] [TIME_LIMIT_S]

For each cell with max sequence length <= MAX_K, solve an ILP that chooses
only aligned ordering edges plus triple sums. Report success rate and
whether the optimum hits the predicted (k1+k2+k3-1)/2 row count.
"""

import sys
import json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds


def objective_vec(k):
    return np.array([(-1) ** j for j in range(k)])


def aligned_edges_for(pat):
    """Return list of aligned ordering edges (a,b) for pattern pat."""
    k = len(pat)
    obj = objective_vec(k)
    ro = list(np.argsort(np.array(pat)))
    edges = []
    for j in range(k - 1):
        a, b = ro[j], ro[j + 1]
        if obj[a] == 1 and obj[b] == -1:
            edges.append((a, b))
    return edges


def has_aligned_certificate(ys, time_limit=5.0):
    """ILP over aligned ordering edges + triple sums. Returns (ok, n_rows, target)."""
    ks = [len(y) for y in ys]
    n = sum(ks)
    obj_seqs = [objective_vec(k) for k in ks]
    offs = [0, ks[0], ks[0] + ks[1]]

    # variables: one per aligned edge, then one per triple sum
    edge_meta = []      # (seq_index, a, b)
    edge_cols = []      # n-dim coefficient column
    for i, y in enumerate(ys):
        pat = tuple(np.argsort(np.argsort(y)) + 1)
        for a, b in aligned_edges_for(pat):
            col = np.zeros(n)
            col[offs[i] + a] = 1
            col[offs[i] + b] = -1
            edge_cols.append(col)
            edge_meta.append((i, a, b))

    S = np.sign(ys[0][:, None, None] + ys[1][None, :, None] + ys[2][None, None, :])
    triple_meta = []
    triple_cols = []
    for p in range(ks[0]):
        for q in range(ks[1]):
            for r in range(ks[2]):
                sgn = int(S[p, q, r])
                col = np.zeros(n)
                col[p] = -sgn
                col[offs[1] + q] = -sgn
                col[offs[2] + r] = -sgn
                triple_cols.append(col)
                triple_meta.append((p, q, r, sgn))

    A_eq = np.array(edge_cols + triple_cols).T  # n x m
    target = np.concatenate(obj_seqs).astype(float)
    m = A_eq.shape[1]
    con = LinearConstraint(A_eq, target, target)
    res = milp(np.ones(m), constraints=con, integrality=np.ones(m),
               bounds=Bounds(np.zeros(m), np.ones(m)),
               options=dict(time_limit=time_limit))
    if not res.success:
        return False, None, (sum(ks) - 1) // 2
    rows = int(round(res.fun))
    return rows == (sum(ks) - 1) // 2, rows, (sum(ks) - 1) // 2


def pattern_triple(ys):
    return tuple(tuple(int(r) for r in (np.argsort(np.argsort(y)) + 1)) for y in ys)


def test_file(cells_json, max_k=None, time_limit=5.0, out_json=None):
    data = json.load(open(cells_json))
    cells = data.get("cells", [])
    print(f"loaded {len(cells)} cells")
    results = []
    for rec in cells:
        ys = [np.array(y, dtype=float) for y in rec["ys"]]
        ks = [len(y) for y in ys]
        if max_k is not None and max(ks) > max_k:
            continue
        ok, rows, target = has_aligned_certificate(ys, time_limit)
        pats = pattern_triple(ys)
        results.append(dict(ks=ks, ok=bool(ok), rows=rows, target=target,
                            F=rec.get("F"), pats=[list(p) for p in pats]))
    print(f"tested {len(results)} cells with max k <= {max_k}")
    ok_count = sum(1 for r in results if r["ok"])
    print(f"aligned-edge-only certificates: {ok_count}/{len(results)}")
    by_k = {}
    for r in results:
        by_k.setdefault(tuple(r["ks"]), []).append(r["ok"])
    for ks in sorted(by_k):
        vals = by_k[ks]
        print(f"  k={ks}: {sum(vals)}/{len(vals)} ok")
    if out_json:
        json.dump(results, open(out_json, "w"))
        print(f"wrote per-cell results to {out_json}")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "test":
        cells_json = sys.argv[2]
        max_k = int(sys.argv[3]) if len(sys.argv) > 3 else None
        tlim = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
        test_file(cells_json, max_k, tlim)
    elif mode == "export":
        cells_json = sys.argv[2]
        out_json = sys.argv[3]
        max_k = int(sys.argv[4]) if len(sys.argv) > 4 else None
        tlim = float(sys.argv[5]) if len(sys.argv) > 5 else 5.0
        test_file(cells_json, max_k, tlim, out_json)
    else:
        print("usage: cert_test_higher.py {test|export} CELLS.json [OUT.json] [MAX_K] [TIME_LIMIT_S]")
