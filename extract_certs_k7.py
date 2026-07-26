#!/usr/bin/env python3
"""
Extract the actual aligned-edge certificates for the 200 sampled k=(7,7,7)
pattern triples, then look for structural patterns.

Usage:
  python3 -u extract_certs_k7.py extract IN.json OUT.json
  python3 -u extract_certs_k7.py analyze OUT.json
"""

import sys
import json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

from cert_test_higher import aligned_edges_for
from solve_unbalanced import solve_triple


def objective_vec(k):
    return np.array([(-1) ** j for j in range(k)])


def find_certificate(ys, prefer_edges=True, time_limit=10.0):
    """Find an aligned-edge-only certificate for a cell; return edges, sums, n_rows."""
    ks = [len(y) for y in ys]
    n = sum(ks)
    obj_seqs = [objective_vec(k) for k in ks]
    offs = [0, ks[0], ks[0] + ks[1]]

    edge_meta = []
    edge_cols = []
    for i, y in enumerate(ys):
        pat = tuple(int(r) for r in (np.argsort(np.argsort(y)) + 1))
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

    A_eq = np.array(edge_cols + triple_cols).T
    target = np.concatenate(obj_seqs).astype(float)
    m = A_eq.shape[1]
    c = np.ones(m)
    if prefer_edges:
        c[:len(edge_meta)] = 0.999  # prefer aligned edges among min-row solutions
    con = LinearConstraint(A_eq, target, target)
    res = milp(c, constraints=con, integrality=np.ones(m),
               bounds=Bounds(np.zeros(m), np.ones(m)),
               options=dict(time_limit=time_limit))
    if not res.success:
        return None
    x = np.round(res.x).astype(int)
    edges = [edge_meta[j] for j in range(len(edge_meta)) if x[j] == 1]
    sums = [triple_meta[j] for j in range(len(triple_meta)) if x[len(edge_meta) + j] == 1]
    return edges, sums


def extract(in_json, out_json):
    recs = json.load(open(in_json))
    out = []
    for idx, rec in enumerate(recs):
        if "ys" in rec:
            ys = [np.array(y, dtype=float) for y in rec["ys"]]
        else:
            P = [tuple(p) for p in rec["patterns"]]
            F, ys = solve_triple(P, time_limit=30.0)
            if F is None:
                print(f"[{idx}] tri={rec['tri']}: MILP fail")
                continue
        cert = find_certificate(ys)
        if cert is None:
            print(f"[{idx}] tri={rec['tri']}: no aligned-edge cert found")
            continue
        edges, sums = cert
        out.append(dict(
            tri=rec["tri"],
            patterns=rec["patterns"],
            F=rec["F"],
            edges=[[int(i) for i in e] for e in edges],
            sums=[dict(p=int(t[0]), q=int(t[1]), r=int(t[2]), actual_sign=int(t[3])) for t in sums],
            n_ord=len(edges),
            n_sum=len(sums),
        ))
        if (idx + 1) % 25 == 0:
            print(f"processed {idx+1}/{len(recs)}")
    json.dump(out, open(out_json, "w"))
    print(f"wrote {out_json} ({len(out)} records)")


def analyze(in_json):
    recs = json.load(open(in_json))
    from collections import Counter
    comp = Counter((r["n_ord"], r["n_sum"]) for r in recs)
    print(f"records: {len(recs)}")
    print("compositions (n_ord, n_sum):")
    for c, cnt in sorted(comp.items()):
        print(f"  {c}: {cnt}")
    # distribution of number of aligned edges per sequence
    edge_counts = [0, 0, 0]
    for r in recs:
        by_seq = [0, 0, 0]
        for e in r["edges"]:
            by_seq[e[0]] += 1
        for i in range(3):
            edge_counts[i] += by_seq[i]
    print("average aligned edges per sequence:", [ec / len(recs) for ec in edge_counts])
    # total edges used distribution
    total_edges = Counter(r["n_ord"] for r in recs)
    print("total aligned edges distribution:", dict(sorted(total_edges.items())))


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "extract":
        extract(sys.argv[2], sys.argv[3])
    elif mode == "analyze":
        analyze(sys.argv[2])
    else:
        print("usage: extract_certs_k7.py {extract|analyze} IN.json [OUT.json]")
