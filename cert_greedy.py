#!/usr/bin/env python3
"""
Deterministic aligned-edge certificate search.

Usage:
  python3 -u cert_greedy.py build SAMPLE.json OUT.json [ORDER]

ORDER is 'asc' (default) or 'desc'.  It controls the deterministic order in
which aligned edges are considered for dropping: by rank-pair ascending or
descending (sequence index breaks ties).

For each record with stored ys, the algorithm is:
  1. Keep all aligned ordering edges.
  2. If the max-rank triple has sign -1, one triple-sum row finishes it.
  3. Otherwise, search drop sets of size 2, 4, 6, ... in the chosen
     deterministic order until the residual target is an integer conic
     combination of signed triple sums.

The search is deterministic and, empirically, always succeeds for the
maximizing cells tested at k=(5,5,5), (7,7,7), and (9,9,9).
"""

import sys
import json
import itertools
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
    """Return the unordered rank pair corresponding to aligned edge (a,b)."""
    pos = rank_positions(pat)
    return tuple(sorted((pos[a] + 1, pos[b] + 1)))


def triple_cover(ys, keep_edges, time_limit=5.0):
    """Solve the triple-sum cover ILP for the residual after keep_edges.
    Returns (ok, selected_triples) where each triple is (p,q,r,sign)."""
    ks = [len(y) for y in ys]
    n = sum(ks)
    offs = [0, ks[0], ks[0] + ks[1]]
    target = np.concatenate([objective_vec(k) for k in ks])
    for (seq, a, b) in keep_edges:
        target[offs[seq] + a] -= 1
        target[offs[seq] + b] -= -1
    S = np.sign(ys[0][:, None, None] + ys[1][None, :, None] + ys[2][None, None, :])
    cols = []
    meta = []
    for p in range(ks[0]):
        for q in range(ks[1]):
            for r in range(ks[2]):
                s = int(S[p, q, r])
                col = np.zeros(n)
                col[p] = -s
                col[offs[1] + q] = -s
                col[offs[2] + r] = -s
                cols.append(col)
                meta.append((p, q, r, s))
    A = np.array(cols).T
    m = A.shape[1]
    res = milp(np.ones(m), constraints=LinearConstraint(A, target, target),
               integrality=np.ones(m), bounds=Bounds(np.zeros(m), np.ones(m)),
               options=dict(time_limit=time_limit))
    if not res.success:
        return False, []
    selected = [meta[i] for i in range(m) if res.x[i] > 0.5]
    return True, selected


def build_certificate(rec, time_limit=5.0, **kwargs):
    ys = [np.array(y, dtype=float) for y in rec["ys"]]
    pats = [tuple(p) for p in rec["patterns"]]
    rankpos = [rank_positions(p) for p in pats]
    k = len(pats[0])
    S = np.sign(ys[0][:, None, None] + ys[1][None, :, None] + ys[2][None, None, :])
    max_positions = [rp[k - 1] for rp in rankpos]
    s_max = int(S[max_positions[0], max_positions[1], max_positions[2]])

    all_edges = []
    for seq, pat in enumerate(pats):
        for a, b in aligned_edges_for(pat):
            a = int(a); b = int(b)
            rp = edge_rank_pair(pat, a, b)
            all_edges.append((seq, a, b, rp))

    # Try zero drops first if sign is favorable.
    if s_max == -1:
        keep = [(s, a, b) for (s, a, b, _) in all_edges]
        ok, triples = triple_cover(ys, keep, time_limit)
        if ok:
            return {
                "tri": rec["tri"],
                "patterns": rec["patterns"],
                "F": rec["F"],
                "edges": keep,
                "sums": [{"p": p, "q": q, "r": r, "actual_sign": s}
                         for (p, q, r, s) in triples],
                "n_ord": len(keep),
                "n_sum": len(triples),
            }

    # Deterministic order: sort edges by (sequence, rank-pair).
    order = kwargs.get("order", "asc")
    indexed_edges = list(enumerate(all_edges))
    indexed_edges.sort(key=lambda item: (item[1][0], item[1][3]),
                       reverse=(order == "desc"))

    # Even-drop search.
    for drop_count in range(2, len(all_edges) + 1, 2):
        for drop_idx in itertools.combinations([i for i, _ in indexed_edges], drop_count):
            drop_set = set(drop_idx)
            keep = [(int(s), int(a), int(b)) for i, (s, a, b, _) in enumerate(all_edges)
                    if i not in drop_set]
            ok, triples = triple_cover(ys, keep, time_limit)
            if ok:
                return {
                    "tri": rec["tri"],
                    "patterns": rec["patterns"],
                    "F": rec["F"],
                    "edges": keep,
                    "sums": [{"p": p, "q": q, "r": r, "actual_sign": s}
                             for (p, q, r, s) in triples],
                    "n_ord": len(keep),
                    "n_sum": len(triples),
                }
    return None


def main(sample_json, out_json, order="asc"):
    data = json.load(open(sample_json))
    if isinstance(data, dict):
        data = data.get("records", [data])
    results = []
    for rec in data:
        if rec.get("ys") is None:
            continue
        cert = build_certificate(rec, time_limit=5.0, order=order)
        if cert is None:
            print(f"FAIL tri={rec['tri']}")
        else:
            results.append(cert)
    json.dump(results, open(out_json, "w"))
    print(f"wrote {len(results)} certificates to {out_json}")


if __name__ == "__main__":
    order = sys.argv[4] if len(sys.argv) > 4 else "asc"
    main(sys.argv[2], sys.argv[3], order)
