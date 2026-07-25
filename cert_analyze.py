#!/usr/bin/env python3
"""
Human-readable certificate anatomy for the k=(5,5,5) exhaustive run.

Usage:
  python3 -u cert_analyze.py all [OUT.txt]     # all 120 records
  python3 -u cert_analyze.py detail IDX        # one record
  python3 -u cert_analyze.py summary           # compact table
"""

import sys
import json
import numpy as np
from collections import Counter

import exhaustive_k5 as ex
from dual_cert import cell_lp_with_duals

K = 5
OBJ = np.array([1, -1, 1, -1, 1])


def load_records():
    return json.load(open("ex_a.json")) + json.load(open("ex_b.json")) \
         + json.load(open("ex_c.json"))


def rank_order(pat):
    return list(np.argsort(np.array(pat)))


def is_aligned_edge(a, b):
    """Edge (a,b) contributes +e_a - e_b; aligned if this matches the
    objective signs at a and b (i.e. a even, b odd)."""
    return OBJ[a] == 1 and OBJ[b] == -1


def residual_for_seq(i, pat, active_ords):
    r = OBJ.copy()
    for kind, meta, lam in active_ords:
        if kind != "ord":
            continue
        si, a, b = meta
        if si != i:
            continue
        r[a] -= 1
        r[b] += 1
    return r


def analyze_record(idx):
    pats_all, triples = ex.triple_list()
    recs = load_records()
    rec = recs[idx]
    tri = rec["tri"]
    P = [pats_all[i] for i in tri]
    y = np.array(rec["y"])
    ys = [y[0:5], y[5:10], y[10:15]]
    out = cell_lp_with_duals(ys, P)
    duals = out["duals"]
    active_ords = [d for d in duals if d[0] == "ord"]
    active_sums = [d for d in duals if d[0] == "sum"]
    comp = (len(active_ords), len(active_sums))

    per_seq = []
    for i, pat in enumerate(P):
        ro = rank_order(pat)
        edges = []
        for d in active_ords:
            si, a, b = d[1]
            if si == i:
                edges.append((a, b, is_aligned_edge(a, b)))
        r = residual_for_seq(i, pat, active_ords)
        residual_pos = [j for j in range(K) if r[j] != 0]
        per_seq.append(dict(
            pattern=pat,
            rank_order=ro,
            edges=edges,
            residual=r,
            residual_pos=residual_pos,
        ))

    sum_rows = []
    for d in active_sums:
        _, (p, q, r, sgn), lam = d
        sum_rows.append(dict(
            p=p, q=q, r=r,
            actual_sign=sgn,
            eff_sign=-sgn,
        ))

    return dict(
        idx=idx,
        tri=tri,
        patterns=P,
        comp=comp,
        per_seq=per_seq,
        sums=sum_rows,
        F=out["F"],
    )


def format_record(a):
    lines = []
    lines.append(f"--- record {a['idx']:3d}  tri={a['tri']}  comp={a['comp']}  F={a['F']:.6e} ---")
    for i, s in enumerate(a["per_seq"]):
        lines.append(f"  seq{i}: pattern={s['pattern']}  rank_order={s['rank_order']}")
        edge_strs = []
        for a_e, b_e, aligned in s["edges"]:
            edge_strs.append(f"({a_e},{b_e}){'*' if aligned else ''}")
        lines.append(f"         active ord edges: " + ", ".join(edge_strs))
        r = s["residual"]
        rpos = [(j, int(r[j])) for j in s["residual_pos"]]
        lines.append(f"         residual vector={r.tolist()}  positions={rpos}")
    lines.append("  active triple sums (eff. sign):")
    for t in a["sums"]:
        lines.append(f"      ({t['p']},{t['q']},{t['r']})  actual={t['actual_sign']:+d}  eff={t['eff_sign']:+d}")
    lines.append("")
    return "\n".join(lines)


def compact_summary():
    recs = load_records()
    rows = []
    comps = Counter()
    for idx in range(len(recs)):
        a = analyze_record(idx)
        comps[a["comp"]] += 1
        rows.append(a)
    print(f"Total records: {len(recs)}")
    print(f"Compositions: {dict(comps)}")
    print("\nIdx | tri                | comp | residual runs (seq0/seq1/seq2)")
    print("-" * 80)
    for a in rows:
        runs = "/".join(str(s["residual_pos"]) for s in a["per_seq"])
        print(f"{a['idx']:3d} | {str(a['tri']):18s} | {a['comp']}  | {runs}")


def main():
    mode = sys.argv[1]
    if mode == "detail":
        idx = int(sys.argv[2])
        a = analyze_record(idx)
        print(format_record(a))
    elif mode == "all":
        out = sys.stdout if len(sys.argv) < 3 else open(sys.argv[2], "w")
        for idx in range(len(load_records())):
            out.write(format_record(analyze_record(idx)))
        if out is not sys.stdout:
            out.close()
            print(f"wrote {sys.argv[2]}")
    elif mode == "summary":
        compact_summary()
    else:
        print("unknown mode")


if __name__ == "__main__":
    main()
