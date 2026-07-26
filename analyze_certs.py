#!/usr/bin/env python3
"""
Analyze aligned-edge certificates.

Usage:
  python3 -u analyze_certs.py min_drops SAMPLE.json

For each record (must contain ys), compute the minimum number of aligned
ordering edges that must be dropped so that the residual target is an
integer conic combination of signed triple-sum rows.  Reports the
distribution and the correlation with the sign of the max-rank triple.
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
    """Return the position of each rank 1..k in pattern pat."""
    pos = [0] * len(pat)
    for i, r in enumerate(pat):
        pos[r - 1] = i
    return pos


def residual_coverable(ys, keep_edges, time_limit=5.0):
    ks = [len(y) for y in ys]
    n = sum(ks)
    offs = [0, ks[0], ks[0] + ks[1]]
    target = np.concatenate([objective_vec(k) for k in ks])
    for (seq, a, b) in keep_edges:
        target[offs[seq] + a] -= 1
        target[offs[seq] + b] -= -1
    S = np.sign(ys[0][:, None, None] + ys[1][None, :, None] + ys[2][None, None, :])
    cols = []
    for p in range(ks[0]):
        for q in range(ks[1]):
            for r in range(ks[2]):
                s = int(S[p, q, r])
                col = np.zeros(n)
                col[p] = -s
                col[offs[1] + q] = -s
                col[offs[2] + r] = -s
                cols.append(col)
    A = np.array(cols).T
    m = A.shape[1]
    res = milp(np.ones(m), constraints=LinearConstraint(A, target, target),
               integrality=np.ones(m), bounds=Bounds(np.zeros(m), np.ones(m)),
               options=dict(time_limit=time_limit))
    return res.success


def min_drops_and_max_sign(rec, time_limit=5.0):
    ys = [np.array(y, dtype=float) for y in rec["ys"]]
    pats = [tuple(p) for p in rec["patterns"]]
    rankpos = [rank_positions(p) for p in pats]
    k = len(pats[0])
    S = np.sign(ys[0][:, None, None] + ys[1][None, :, None] + ys[2][None, None, :])
    s_max = int(S[rankpos[0][k - 1], rankpos[1][k - 1], rankpos[2][k - 1]])

    all_edges = []
    for seq, pat in enumerate(pats):
        for a, b in aligned_edges_for(pat):
            all_edges.append((seq, a, b))

    for drop in range(len(all_edges) + 1):
        for drop_set in itertools.combinations(all_edges, drop):
            keep = [e for e in all_edges if e not in drop_set]
            if residual_coverable(ys, keep, time_limit):
                return drop, s_max, keep
    return None, s_max, None


def main(path):
    data = json.load(open(path))
    if isinstance(data, dict):
        data = data.get("records", [data])
    print(f"loaded {len(data)} records")
    dist = {}
    sgn_dist = {}
    for rec in data:
        drop, s_max, keep = min_drops_and_max_sign(rec, time_limit=5.0)
        dist[(drop, s_max)] = dist.get((drop, s_max), 0) + 1
        sgn_dist[s_max] = sgn_dist.get(s_max, 0) + 1
    print("sign of max-rank triple:", sgn_dist)
    print("(min_drops, sgn_max) distribution:")
    for k in sorted(dist):
        print(f"  {k}: {dist[k]}")


if __name__ == "__main__":
    main(sys.argv[2])
