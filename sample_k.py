#!/usr/bin/env python3
"""
Random spot-check of the aligned-edge conjecture at arbitrary odd k.

Usage:
  python3 -u sample_k.py run K N OUT.json [SEED] [TIME_LIMIT_S]

Solves the global MILP for N random symmetry-reduced pattern triples at size k
and tests whether each optimum has an aligned-edge-only unit certificate.
Stores the optimal ys so certificates can be re-extracted later.
"""

import sys
import json
import itertools
import numpy as np

from conj6_search import enumerate_patterns
from solve_unbalanced import solve_triple
from cert_test_higher import has_aligned_certificate


def main(k, N, out_json, seed=0, time_limit=30.0):
    if k % 2 == 0:
        raise ValueError("k must be odd")
    pats = enumerate_patterns(k)
    print(f"k={k}: {len(pats)} patterns", flush=True)
    total = len(list(itertools.combinations_with_replacement(range(len(pats)), 3)))
    print(f"symmetry-reduced triples: {total}", flush=True)
    N = min(N, total)
    rng = np.random.default_rng(seed)
    all_indices = list(itertools.combinations_with_replacement(range(len(pats)), 3))
    chosen = rng.choice(len(all_indices), size=N, replace=False)
    results = []
    for idx in chosen:
        tri = all_indices[idx]
        P = [pats[i] for i in tri]
        F, ys = solve_triple(P, time_limit=time_limit)
        if F is None:
            results.append(dict(k=k, tri=list(tri), patterns=[list(p) for p in P],
                                F=None, ys=None, ok=False, reason="milp_fail"))
            print(f"  tri={tri}: MILP fail", flush=True)
            continue
        ok, rows, target = has_aligned_certificate(ys, time_limit=10.0)
        results.append(dict(k=k, tri=list(tri), patterns=[list(p) for p in P],
                            F=F, ys=[y.tolist() for y in ys],
                            ok=bool(ok), rows=rows, target=target))
        print(f"  tri={tri}: F={F:.6e} ok={ok}", flush=True)
    ok_count = sum(1 for r in results if r["ok"])
    print(f"\naligned-edge certs at global optimum: {ok_count}/{len(results)}", flush=True)
    json.dump(results, open(out_json, "w"))
    print(f"wrote {out_json}", flush=True)


if __name__ == "__main__":
    k = int(sys.argv[2])
    N = int(sys.argv[3])
    out = sys.argv[4]
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    time_limit = float(sys.argv[6]) if len(sys.argv) > 6 else 30.0
    main(k, N, out, seed, time_limit)
