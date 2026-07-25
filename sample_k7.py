#!/usr/bin/env python3
"""
Systematic spot-check of the aligned-edge conjecture at k=(7,7,7).

Usage:
  python3 -u sample_k7.py run N OUT.json [SEED]

Solves the global MILP for N random symmetry-reduced pattern triples and
tests whether the optimum has an aligned-edge-only unit certificate.
"""

import sys
import json
import itertools
import numpy as np

from conj6_search import enumerate_patterns
from solve_unbalanced import solve_triple
from cert_test_higher import has_aligned_certificate


def main(N, out_json, seed=0):
    pats = enumerate_patterns(7)
    print(f"k=7: {len(pats)} patterns")
    # symmetry-reduced triples as sorted indices
    total = len(list(itertools.combinations_with_replacement(range(len(pats)), 3)))
    print(f"symmetry-reduced triples: {total}")
    N = min(N, total)
    rng = np.random.default_rng(seed)
    all_indices = list(itertools.combinations_with_replacement(range(len(pats)), 3))
    chosen = rng.choice(len(all_indices), size=N, replace=False)
    results = []
    for idx in chosen:
        tri = all_indices[idx]
        P = [pats[i] for i in tri]
        F, ys = solve_triple(P, time_limit=10.0)
        if F is None:
            results.append(dict(tri=list(tri), patterns=[list(p) for p in P],
                                F=None, ok=False, reason="milp_fail"))
            print(f"  tri={tri}: MILP fail")
            continue
        ok, rows, target = has_aligned_certificate(ys, time_limit=5.0)
        results.append(dict(tri=list(tri), patterns=[list(p) for p in P],
                            F=F, ok=bool(ok), rows=rows, target=target))
        print(f"  tri={tri}: F={F:.6e} ok={ok}")
    ok_count = sum(1 for r in results if r["ok"])
    print(f"\naligned-edge certs at global optimum: {ok_count}/{len(results)}")
    json.dump(results, open(out_json, "w"))
    print(f"wrote {out_json}")


if __name__ == "__main__":
    N = int(sys.argv[2])
    out = sys.argv[3]
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    main(N, out, seed)
