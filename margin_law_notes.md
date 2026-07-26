# The margin law and its certificates: what is now PROVEN vs observed

Status of the sharp form of Tao's Conjecture 6:

    F  <=  -((k1+k2+k3-1)/2) * margin        (margin = min strict quantity)

## Proven (hand proofs, general odd k1,k2,k3, any pattern)

### Lemma 1 (boundary lemma, bottom). Axiom (iii) forces
    sgn(y_{1,1} + y_{2,1} + y_{3,1}) = -1.
Proof. Apply the (iii) window (p,q,r) = (0,0,0) (even parity, always allowed).
Corners (a,b,c) in {0,1}^3; index 0 is the sentinel -infty. Contributions:
all-sentinel corner: factor +1, sgn = -1, gives -1; one finite index (two
sentinels): factor -1, sgn -1, three corners, gives +3; two finite (one
sentinel): factor +1, sgn -1, three corners, gives -3; all finite (1,1,1):
factor -1, sgn = tau, gives -tau. Total: -1 + 3 - 3 - tau = -1 - tau.
Axiom (iii) demands 0, so tau = -1. QED.

### Lemma 2 (boundary lemma, top). Axiom (iii) forces
    sgn(y_{1,k1} + y_{2,k2} + y_{3,k3}) = -1.
Proof. Window (p,q,r) = (k1,k2,k3): all k_i odd, so same parity; indices
k_i + 1 are sentinels. Corner (k1,k2,k3): factor (-1)^(odd sum) = -1,
sgn = sigma. One sentinel: factor +1, sgn -1, x3 = -3. Two sentinels:
factor -1, sgn -1, x3 = +3. Three sentinels: factor +1, sgn -1 = -1.
Total: -sigma - 3 + 3 - 1 = -sigma - 1 = 0, so sigma = -1. QED.

### Theorem (identity patterns). If all three sequences are increasing,
Conjecture 6 holds for all odd k1,k2,k3, with the sharp constant:
    F <= -((k1+k2+k3-1)/2) * margin.
Proof. F = sum_i [ (y_{i,1}-y_{i,2}) + (y_{i,3}-y_{i,4}) + ...
        + (y_{i,k_i-2}-y_{i,k_i-1}) ] + (y_{1,k1} + y_{2,k2} + y_{3,k3}).
Each of the sum_i (k_i - 1)/2 bracketed descents is <= -margin because the
sequence increases; the final triple sum is <= -margin by Lemma 2 plus
axiom (i)'s strictness. Term count: (k1-1)/2 + (k2-1)/2 + (k3-1)/2 + 1
= (k1+k2+k3-1)/2. QED.

This matches the LP-observed optimum -7 EPS at k=(5,5,5) exactly and proves
that constant is SHARP for identity patterns (take equality-approaching
configurations).

### Proposition (parity of the certificate). In any unit-coefficient
certificate, the number of triple-sum terms is odd, and (number with sign -1)
= (number with sign +1) + 1.
Proof. Sum the entries of the vector identity "objective = sum of active
rows". The objective's entries sum to 3 (each sequence contributes
+1-1+...+1 = 1). Ordering rows sum to 0. A triple-sum row with cell sign s
sums to -3s. Hence 3 = -3 * sum(s_j), so sum(s_j) = -1. QED.
(Observed compositions at k=5: n_sum in {1,3,5} only. Consistent.)

## Verified exactly by machine (no floating point in the final check)

For ALL 120 symmetry-reduced pattern triples at k=(5,5,5): the LP dual of the
optimal cell has all multipliers EXACTLY 1, exactly 7 active constraints, box
bounds inactive, and the INTEGER vector identity
    (objective vector) = (sum of the 7 active constraint rows)
holds exactly (entries in {-1,0,1}; checked in integer arithmetic;
`cert_class.json`). Each active row is <= -margin inside the cell, so
F <= -7 margin follows exactly for every point of each of those cells.

Scope caveat: this settles the OPTIMAL cell of each pattern triple exactly.
The full base-case verdict over all cells rests on the float MILP; the
certificates above remove float doubt only where they exist.

### Aligned-edge restriction (newly verified)

All 120 optimal-cell certificates can be chosen to use ONLY **aligned**
ordering edges: consecutive-rank pairs (a,b) with rank(a) < rank(b),
OBJ[a] = +1, OBJ[b] = -1. No LP dual uses a misaligned ordering edge.
The identity-pattern theorem is the special case that uses all aligned
edges plus one boundary triple sum. For non-identity patterns, a subset
of aligned edges is kept and the residual positions are covered by a
telescoping chain of triple-sum terms.

A brute-force subset search over aligned edges finds a valid certificate
for every k=(5,5,5) record, reproducing the (6,1)/(4,3)/(2,5)
composition counts. Higher-k tests refine this:

- Among 1,103 inflation cells with max sequence length ≤ 15, 98.5% have
  aligned-edge-only certificates. The few failures are suboptimal cells
  inside their pattern triple.
- For every pattern triple where no ladder cell had an aligned-edge cert,
  solving the global MILP for that triple produces an optimum that DOES
  have an aligned-edge-only certificate.
- 200 random k=(7,7,7) pattern triples solved to global optima all have
  aligned-edge-only certificates.

This motivates the refined conjecture:

**Refined certificate-existence conjecture.** For every pattern triple,
the cell that maximizes F admits a unit-coefficient certificate using
only aligned ordering edges plus triple sums.

This is the form needed to clear a whole pattern triple, because only
the maximizing cell controls the sign of sup F.

### Parity theorem for minimal aligned-edge drops

Let k be odd and consider a balanced cell (k,k,k). Write the objective
sign vector as `OBJ[j] = (-1)^j`. An aligned ordering edge (a,b) contributes
`+1` at the even position a and `-1` at the odd position b. If `e_i`
aligned edges are kept in sequence i, the residual signed difference

    D_i = (sum over even j of residual[j]) - (sum over odd j of residual[j])

is `1 - 2·e_i`, because each kept edge contributes equally and oppositely
to the even and odd sums.

Every triple-sum row contributes `±1` at exactly three positions, one in each
sequence. Modulo 2, its contribution to `D_i` is therefore `1`, regardless of
sign or of which positions are hit. Hence, if `N` triple-sum rows are used,
`D_i ≡ N (mod 2)` for each i. Since `D_i = 1 - 2·e_i` is odd, `N` must be odd.

At size (k,k,k) the margin law gives target row count `(3k-1)/2`. The total
number of aligned edges available is `3·(k-1)/2`, so

    N = (3k-1)/2 - (3·(k-1)/2 - drops) = 1 + drops.

Thus `N` is odd **iff the number of dropped aligned edges is even**.
Empirically this is sharp: minimal-drop counts are always 0, 2, 4, … .

**Sign condition for zero drops.** With all aligned edges kept, the only
residual position in each sequence is the position of the largest rank k.
The residual is coverable by a single triple-sum row exactly when the sign
of `y^1_{p_1}+y^2_{p_2}+y^3_{p_3}` (with `p_i` the position of rank k in
sequence i) is negative. This happens in about 1/3 of k=5 and 1/6 of k=7
random optima; otherwise an even number of edges must be dropped.

**Constructive selection rule.** The parity theorem gives a deterministic
search: keep all aligned edges, use the max-rank triple if its sign is
negative, otherwise try dropping 2, 4, 6, ... aligned edges in a fixed order
until the residual becomes an integer conic combination of signed triple
sums. This rule succeeds on every tested maximizing cell (k=5, 7, and 9).

A faster implementation uses a two-stage MILP: stage 1 minimizes the total
number of active rows, and stage 2 breaks ties with distinct powers-of-two
weights on the aligned-edge drop variables. The resulting drop set is unique
and deterministic, and the builder is 10–50x faster than the explicit
even-drop search while producing the same minimal drop counts.

**Unimodal sign-tensor observation.** For every tested maximizing cell, the
sign tensor `S(p,q,r) = sgn(y^1_p + y^2_q + y^3_r)` is unimodal along each
coordinate when positions are ordered by *rank* (not by index): for fixed
values of the other two coordinates, `S` takes one sign and then, at most
once, switches to the other sign as the rank increases. This is verified for
all k=5, k=7, and k=9 records. It is the 3D analogue of the Monge property
and is the likely source of the residual tileability: a signed 3D matching in
a unimodal tensor should admit a greedy decomposition.

## Observed, not yet proven

- The margin law itself across all patterns: verified exhaustively at
  (5,5,5), spot-checked (7,7,7), (9,9,9), (5,5,7), (5,7,9). No misses.
- Certificate existence with unit coefficients for EVERY valid cell (not
  just optimal ones), i.e. the combinatorial decomposition conjecture:
  every valid cell admits a partition of F into (k1+k2+k3-1)/2 terms, each
  an ordering descent or a (iii)-forced signed triple sum, each negative.
- Insertion adds exactly one term to the certificate (the induction step;
  consistent with the ladder's -1 EPS drift per inserted pair).

## Attack plan for the general case

Non-identity patterns replace some descents with telescoping chains of
triple sums (observed: sums sharing indices, alternating signs). The
natural object: a path/matching in a graph whose vertices are the sequence
positions and whose edges are descents and (iii)-forced sums. Integrality
of all 120 duals suggests the constraint matrix restricted to valid cells
is totally unimodular or the certificate is a min-max dual of a
combinatorial packing. Next concrete steps:
1. For the (4,3) and (2,5) composition classes, write out the telescopes and
   find the general rule mapping pattern -> certificate.
2. Prove the interior sum-sign lemmas generalizing Lemmas 1-2: which
   (iii) windows force which interior triple-sum signs, given the pattern.
3. Induction: inserted arc contributes its own descent term (or extends a
   telescope by one link).
