#!/usr/bin/env python3
"""
Deterministic aligned-edge certificate via two-stage MILP.

Usage:
  python3 -u cert_deterministic.py build SAMPLE.json OUT.json [ORDER]

ORDER is 'asc' (default) or 'desc'.  Aligned edges are assigned distinct
small tie-breaking weights in that order; the two-stage MILP first minimizes
total rows, then picks the unique optimal certificate consistent with the
chosen order.

This is faster than the explicit even-drop search in cert_greedy.py while
still being fully deterministic.
"""

import sys
import json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from cert_test_higher import aligned_edges_for


def objective_vec(k):
    return np.array([(-1) ** j for j in range(k)])


def rank_positions(pat):
    pos = [0] * len(pat)
    for i, r in enumerate(pat):
        pos[r - 1] = i
    return pos


def edge_rank_pair(pat, a, b):
    pos = rank_positions(pat)
    return tuple(sorted((pos[a] + 1, pos[b] + 1)))


def build_certificate(rec, order="asc", time_limit=10.0):
    ys = [np.array(y, dtype=float) for y in rec["ys"]]
    pats = [tuple(p) for p in rec["patterns"]]
    ks = [len(y) for y in ys]
    n = sum(ks)
    offs = [0, ks[0], ks[0] + ks[1]]
    target = np.concatenate([objective_vec(k) for k in ks])

    # Build aligned-edge columns.
    edge_meta = []   # (seq, a, b)
    edge_cols = []
    for seq, pat in enumerate(pats):
        for a, b in aligned_edges_for(pat):
            a = int(a); b = int(b)
            col = np.zeros(n)
            col[offs[seq] + a] = 1
            col[offs[seq] + b] = -1
            edge_cols.append(col)
            edge_meta.append((seq, a, b, edge_rank_pair(pat, a, b)))
    n_edges = len(edge_cols)

    # Build triple-sum columns.
    S = np.sign(ys[0][:, None, None] + ys[1][None, :, None] + ys[2][None, None, :])
    triple_meta = []
    triple_cols = []
    for p in range(ks[0]):
        for q in range(ks[1]):
            for r in range(ks[2]):
                s = int(S[p, q, r])
                col = np.zeros(n)
                col[p] = -s
                col[offs[1] + q] = -s
                col[offs[2] + r] = -s
                triple_cols.append(col)
                triple_meta.append((p, q, r, s))
    n_triples = len(triple_cols)

    # Full constraint matrix: columns = [edge_vars, triple_vars].
    A = np.array(edge_cols + triple_cols).T  # n x (n_edges + n_triples)
    m = A.shape[1]
    con = LinearConstraint(A, target, target)
    integrality = np.ones(m)
    lb = np.zeros(m)
    ub = np.ones(m)

    # Deterministic tie-breaking weights for edge drops.
    idx_order = list(range(n_edges))
    idx_order.sort(key=lambda i: (edge_meta[i][0], edge_meta[i][3]),
                   reverse=(order == "desc"))
    # Use powers of two so every weighted sum is unique.
    edge_weights = np.zeros(n_edges)
    for rank, i in enumerate(idx_order):
        edge_weights[i] = 2.0 ** rank

    # Stage 1: minimize total rows (kept edges count as used rows too).
    c1 = np.concatenate([np.ones(n_edges), np.ones(n_triples)])
    res1 = milp(c1, constraints=con, integrality=integrality,
                bounds=Bounds(lb, ub), options=dict(time_limit=time_limit))
    if not res1.success:
        return None
    n_rows_opt = int(round(res1.fun))

    # Stage 2: among minimum-row solutions, pick the one with the chosen order.
    row_sum = np.ones(m)
    con2 = [con, LinearConstraint(row_sum, n_rows_opt, n_rows_opt)]
    c2 = np.concatenate([edge_weights, np.zeros(n_triples)])
    res2 = milp(c2, constraints=con2, integrality=integrality,
                bounds=Bounds(lb, ub), options=dict(time_limit=time_limit))
    if not res2.success:
        return None
    x = res2.x
    keep_edges = [(seq, int(a), int(b)) for i, (seq, a, b, _) in enumerate(edge_meta)
                  if x[i] > 0.5]
    triples = [triple_meta[i] for i in range(n_triples) if x[n_edges + i] > 0.5]
    return {
        "tri": rec["tri"],
        "patterns": rec["patterns"],
        "F": rec["F"],
        "edges": keep_edges,
        "sums": [{"p": int(p), "q": int(q), "r": int(r), "actual_sign": int(s)}
                 for (p, q, r, s) in triples],
        "n_ord": len(keep_edges),
        "n_sum": len(triples),
    }


def main(sample_json, out_json, order="asc"):
    data = json.load(open(sample_json))
    if isinstance(data, dict):
        data = data.get("records", [data])
    results = []
    for rec in data:
        if rec.get("ys") is None:
            continue
        cert = build_certificate(rec, order=order, time_limit=10.0)
        if cert is None:
            print(f"FAIL tri={rec['tri']}")
        else:
            results.append(cert)
    json.dump(results, open(out_json, "w"))
    print(f"wrote {len(results)} certificates to {out_json}")


if __name__ == "__main__":
    order = sys.argv[4] if len(sys.argv) > 4 else "asc"
    main(sys.argv[2], sys.argv[3], order)
