#!/usr/bin/env python3
"""
Synthetic certificate tests for k=(5,5,5).

Modes:
  min IDX       - solve a MILP choosing any unit-coefficient subset of
                  ordering + triple-sum rows that reconstructs F.
  min_all       - run `min` for all 120 records; report distribution.
  aligned IDX   - fix ordering rows to ALL aligned adjacent pairs, then
                  solve a MILP for the remaining triple-sum rows.
  aligned_all   - run `aligned` for all records.
"""

import sys
import json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

import exhaustive_k5 as ex

K = 5
OBJ = np.array([1, -1, 1, -1, 1])


def load_records():
    return json.load(open("ex_a.json")) + json.load(open("ex_b.json")) \
         + json.load(open("ex_c.json"))


def rank_order(pat):
    return list(np.argsort(np.array(pat)))


def allowed_edges(pat):
    """List of ordering edges (a,b) where a,b are consecutive in rank order
    and rank(a) < rank(b)."""
    ro = rank_order(pat)
    return [(ro[j], ro[j + 1]) for j in range(len(ro) - 1)]


def is_aligned(a, b):
    return OBJ[a] == 1 and OBJ[b] == -1


def prepare_record(idx):
    pats_all, triples = ex.triple_list()
    recs = load_records()
    rec = recs[idx]
    tri = rec["tri"]
    P = [pats_all[i] for i in tri]
    y = np.array(rec["y"])
    ys = [y[0:5], y[5:10], y[10:15]]
    # actual sign tensor
    S = np.sign(ys[0][:, None, None] + ys[1][None, :, None] + ys[2][None, None, :])
    return P, ys, S


def build_min_milp(P, S):
    """Variables: one per allowed ordering edge + one per triple sum.
    Choose a subset whose signed sum equals the objective vector.
    Minimize number of selected rows."""
    rows = []          # each row: (type, meta, 15-dim coefficient vector)
    # ordering rows
    edge_index = []
    for i, pat in enumerate(P):
        for a, b in allowed_edges(pat):
            vec = np.zeros(15)
            off = 5 * i
            vec[off + a] = 1
            vec[off + b] = -1
            rows.append(("ord", (i, a, b), vec))
            edge_index.append((i, a, b))
    n_ord = len(rows)
    # triple-sum rows
    triple_index = []
    for p in range(K):
        for q in range(K):
            for r in range(K):
                sgn = int(S[p, q, r])
                # actual sign sgn; effective contribution to objective is -sgn
                vec = np.zeros(15)
                vec[p] = -sgn
                vec[5 + q] = -sgn
                vec[10 + r] = -sgn
                rows.append(("sum", (p, q, r, sgn), vec))
                triple_index.append((p, q, r))
    n = len(rows)
    A_eq = np.array([r[2] for r in rows]).T   # 15 x n
    c = np.ones(n)
    integrality = np.ones(n)
    bounds = Bounds(np.zeros(n), np.ones(n))
    # equality A_eq x = OBJ_flat
    target = np.tile(OBJ, 3)
    con = LinearConstraint(A_eq, target, target)
    res = milp(c, constraints=con, integrality=integrality, bounds=bounds,
               options=dict(time_limit=10.0))
    if not res.success:
        return None, None, res.status
    x = np.round(res.x).astype(int)
    selected = [rows[j] for j in range(n) if x[j] == 1]
    return selected, int(res.fun), res.status


def build_aligned_milp(P, S):
    """Fix all aligned ordering edges; choose triple sums to cover residual."""
    target = np.tile(OBJ, 3).astype(float)
    fixed_edges = []
    for i, pat in enumerate(P):
        for a, b in allowed_edges(pat):
            if is_aligned(a, b):
                off = 5 * i
                target[off + a] -= 1
                target[off + b] += 1
                fixed_edges.append((i, a, b))
    # variables: one per triple sum
    cols = []
    meta = []
    for p in range(K):
        for q in range(K):
            for r in range(K):
                sgn = int(S[p, q, r])
                vec = np.zeros(15)
                vec[p] = -sgn
                vec[5 + q] = -sgn
                vec[10 + r] = -sgn
                cols.append(vec)
                meta.append((p, q, r, sgn))
    n = len(cols)
    A_eq = np.array(cols).T
    c = np.ones(n)
    con = LinearConstraint(A_eq, target, target)
    res = milp(c, constraints=con, integrality=np.ones(n),
               bounds=Bounds(np.zeros(n), np.ones(n)),
               options=dict(time_limit=10.0))
    if not res.success:
        return None, fixed_edges, res.status
    x = np.round(res.x).astype(int)
    selected = [meta[j] for j in range(n) if x[j] == 1]
    return selected, fixed_edges, int(res.fun)


def check_certificate(P, S, selected_edges, selected_sums):
    """Verify vector identity exactly for a synthetic certificate."""
    target = np.tile(OBJ, 3).astype(int)
    v = np.zeros(15, dtype=int)
    for i, a, b in selected_edges:
        off = 5 * i
        v[off + a] += 1
        v[off + b] -= 1
    for p, q, r, sgn in selected_sums:
        t = -sgn
        v[p] += t
        v[5 + q] += t
        v[10 + r] += t
    return np.array_equal(v, target)


def report_min(idx):
    P, ys, S = prepare_record(idx)
    sel, val, status = build_min_milp(P, S)
    print(f"record {idx}: status={status}")
    if sel is None:
        print("  no certificate found")
        return
    n_ord = sum(1 for r in sel if r[0] == "ord")
    n_sum = len(sel) - n_ord
    print(f"  selected {n_ord} ord + {n_sum} sum = {val} rows")
    edges = []
    for r in sel:
        if r[0] == "ord":
            i, a, b = r[1]
            edges.append((i, a, b))
        else:
            p, q, rr, sgn = r[1]
            print(f"    sum ({p},{q},{rr}) actual={sgn:+d}")
    for i in range(3):
        e = [(a, b) for (ii, a, b) in edges if ii == i]
        print(f"    seq{i} ord edges: {e}")
    edges = [(i, a, b) for kind, (i, a, b), _ in sel if kind == "ord"]
    sums = [(p, q, r, s) for kind, (p, q, r, s), _ in sel if kind == "sum"]
    ok = check_certificate(P, S, edges, sums)
    print(f"  vector identity ok: {ok}")


def report_aligned(idx):
    P, ys, S = prepare_record(idx)
    sel, fixed, val = build_aligned_milp(P, S)
    print(f"record {idx}: aligned edges fixed = {len(fixed)}")
    if sel is None:
        print("  residual not coverable")
        return
    print(f"  triple sums selected = {val} (total rows = {len(fixed)+val})")
    for i in range(3):
        e = [(a, b) for (ii, a, b) in fixed if ii == i]
        print(f"    seq{i} ord edges: {e}")
    for p, q, r, sgn in sel:
        print(f"    sum ({p},{q},{r}) actual={sgn:+d}")
    ok = check_certificate(P, S, fixed, sel)
    print(f"  vector identity ok: {ok}")


def solve_with_edge_subset(P, S, edges_subset):
    """Given a chosen set of aligned ordering edges, solve for triple sums."""
    target = np.tile(OBJ, 3).astype(float)
    for i, a, b in edges_subset:
        off = 5 * i
        target[off + a] -= 1
        target[off + b] += 1
    cols, meta = [], []
    for p in range(K):
        for q in range(K):
            for r in range(K):
                sgn = int(S[p, q, r])
                vec = np.zeros(15)
                vec[p] = -sgn
                vec[5 + q] = -sgn
                vec[10 + r] = -sgn
                cols.append(vec)
                meta.append((p, q, r, sgn))
    n = len(cols)
    A_eq = np.array(cols).T
    c = np.ones(n)
    con = LinearConstraint(A_eq, target, target)
    res = milp(c, constraints=con, integrality=np.ones(n),
               bounds=Bounds(np.zeros(n), np.ones(n)),
               options=dict(time_limit=5.0))
    if not res.success:
        return None
    x = np.round(res.x).astype(int)
    return [meta[j] for j in range(n) if x[j] == 1]


def aligned_subset_search(P, S):
    """Search over subsets of aligned edges; return first feasible certificate."""
    aligned = []
    for i, pat in enumerate(P):
        edges_i = [(i, a, b) for a, b in allowed_edges(pat) if is_aligned(a, b)]
        aligned.append(edges_i)
    # iterate subsets per sequence (small: at most 2^2=4 each)
    import itertools
    for e0 in itertools.chain.from_iterable(itertools.combinations(aligned[0], r) for r in range(len(aligned[0]) + 1)):
        for e1 in itertools.chain.from_iterable(itertools.combinations(aligned[1], r) for r in range(len(aligned[1]) + 1)):
            for e2 in itertools.chain.from_iterable(itertools.combinations(aligned[2], r) for r in range(len(aligned[2]) + 1)):
                subset = list(e0) + list(e1) + list(e2)
                sums = solve_with_edge_subset(P, S, subset)
                if sums is not None and len(subset) + len(sums) == 7:
                    return subset, sums
    return None, None


def run_all(mode):
    recs = load_records()
    stats = []
    for idx in range(len(recs)):
        P, ys, S = prepare_record(idx)
        if mode == "min":
            sel, val, status = build_min_milp(P, S)
            if sel is None:
                stats.append((idx, None, status))
                continue
            n_ord = sum(1 for r in sel if r[0] == "ord")
            n_sum = len(sel) - n_ord
            edges = []
            sums = []
            for kind, meta, _ in sel:
                if kind == "ord":
                    edges.append(meta)
                else:
                    sums.append(meta)
            ok = check_certificate(P, S, edges, sums)
            stats.append((idx, (n_ord, n_sum), ok))
        elif mode == "aligned":
            sel, fixed, val = build_aligned_milp(P, S)
            ok = check_certificate(P, S, fixed, sel) if sel is not None else False
            stats.append((idx, (len(fixed), val if sel is not None else None), ok))
        elif mode == "aligned_subsets":
            edges, sums = aligned_subset_search(P, S)
            ok = check_certificate(P, S, edges, sums) if edges is not None else False
            n_edges = len(edges) if edges is not None else None
            n_sums = len(sums) if sums is not None else None
            stats.append((idx, (n_edges, n_sums), ok))
    # report
    ok_count = sum(1 for s in stats if s[2] is True)
    fail = [s for s in stats if not s[2]]
    print(f"mode={mode}: ok={ok_count}/{len(recs)}")
    if fail:
        print("failures:", fail[:20])
    else:
        comp_counter = {}
        for s in stats:
            comp = s[1]
            comp_counter[comp] = comp_counter.get(comp, 0) + 1
        print("compositions:", comp_counter)


def export(out_json):
    recs = load_records()
    pats_all, triples = ex.triple_list()
    data = []
    for idx in range(len(recs)):
        P, ys, S = prepare_record(idx)
        edges, sums = aligned_subset_search(P, S)
        data.append(dict(
            idx=idx,
            tri=triples[idx],
            patterns=[list(p) for p in P],
            edges=[[int(i) for i in e] for e in edges],
            sums=[dict(p=t[0], q=t[1], r=t[2], actual_sign=int(t[3])) for t in sums],
            n_ord=len(edges),
            n_sum=len(sums),
        ))
    json.dump(data, open(out_json, "w"))
    print(f"wrote {out_json} ({len(data)} records)")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "export":
        export(sys.argv[2])
    elif mode == "min":
        report_min(int(sys.argv[2]))
    elif mode == "min_all":
        run_all("min")
    elif mode == "aligned":
        report_aligned(int(sys.argv[2]))
    elif mode == "aligned_all":
        run_all("aligned")
    elif mode == "aligned_subsets_all":
        run_all("aligned_subsets")
    else:
        print("usage: cert_synth.py {min|min_all|aligned|aligned_all|aligned_subsets_all|export} [IDX|OUT]")
