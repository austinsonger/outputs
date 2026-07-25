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
composition counts. The same restriction also holds for a tested
k=(7,7,7) cell from the inflation ladder. This motivates the conjecture:

certificate existence = for every valid pattern triple, some subset of
aligned edges leaves a residual that can be tiled by triple sums
consistent with axiom (iii).

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
