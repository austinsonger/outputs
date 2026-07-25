# Inscribed Square Problem: Research Log

Purpose: permanent record of everything attempted on the Toeplitz / Square Peg Problem so no approach, parameter choice, or dead end gets repeated. Append per session, never rewrite history.

## Problem status (context, do not re-derive)

- Conjecture (Toeplitz 1911): every Jordan curve inscribes a square. Open for general continuous curves.
- Proven: convex (Toeplitz), analytic (Jerrard 1961), locally monotone (Stromquist 1989), smooth via every-aspect-ratio rectangles (Greene-Lobb 2020, symplectic). Rectangle existence holds for ALL Jordan curves (Vaughan Mobius-band argument).
- The entire difficulty: on rough curves, inscribed squares in approximating curves can degenerate (side -> 0) in the limit. Any proof must keep one square uniformly large.
- We are NOT attempting a proof of the open case. Experimental/numerical track.

---

## Session 2026-07-25: numerical square-finder + degeneration probe

### Direction chosen
Build a numerical inscribed-square finder and probe degeneration experimentally. (Options not taken this session: proof walkthrough, frontier/lit study, interactive demo.)

### Method (final, works, keep)
Diagonal-pair formulation: square inscribed in curve G <=> parameters (t1,t2) with p1=G(t1), p2=G(t2), c=(p1+p2)/2, w=rot90((p2-p1)/2), such that q1=c+w and q2=c-w both lie on G. Corners (p1,q1,p2,q2) form an exact square by construction; only corners-on-curve needs verifying.

Pipeline in `inscribed_squares.py`:
1. Dense polyline sample (24k pts smooth, 60k rough) + cKDTree.
2. Vectorized grid scan of (t1,t2) torus (460-520 per axis), objective d(q1,G)^2 + d(q2,G)^2, mask |t1-t2| < 0.004 (fake zeros at p1=p2) and tiny diagonals.
3. 8-neighbor torus local minima below (0.03 diam)^2, cap 400, Nelder-Mead refine on exact parametrization.
4. Dedupe by corner sets (0.02 diam tolerance), THEN polish survivors against segment-exact distance.
5. Accept: side >= floor (0.01 diam default) AND residual <= min(1e-7 diam, 0.01 side).

### Failed attempts (do not repeat)
| # | Attempt | Failure | Lesson |
|---|---------|---------|--------|
| 1 | Acceptance tol 1.5e-3·diam, KD distance only | Accepted near-misses: tiny "squares" (side 0.04) with corner error 8% of side | Acceptance must be scale-aware: cap at fraction of SIDE, not just diameter |
| 2 | NM restart + scale-aware filter, still KD point-cloud distance | Residuals stuck ~2e-3 on big squares | KD nearest-sample distance has a resolution floor ~ sample spacing; cannot certify below it |
| 3 | Segment-exact distance polish inside the main NM loop, per-point Python loops | Timeout (>42s for 400 candidates); Python-loop geometry inside an optimizer is 100x too slow | Vectorize distance (batch einsum segment projection); polish AFTER dedupe (40 squares, not 400 candidates) |
| 4 | `nohup ... &` background run in sandbox | Process killed at call boundary, log never written | Sandbox bash calls are independent; background jobs die. Chunk long runs + JSON checkpoints, `python -u`, <=42s per call |

### Verified results (trust these, don't rerun)
- Ellipse a=1.6, b=1.0: exactly 1 square, corners (±0.847998, ±0.847998) = analytic s=ab/sqrt(a²+b²). Residual 3e-12.
- Smooth random blob (seed 7, 6 modes): 1 square. Odd count as parity theorem predicts.
- 5-point star polygon: 5 congruent squares (side 1.122), verified to map onto each other under 72° rotation. Odd.
- Residual distribution is cleanly bimodal: true squares at ~1e-12, near-misses at >=1.7e-5, gap of 7 orders. Cutoff 1e-7·diam sits safely inside the gap.
- 147/147 reported squares across all runs pass exact squareness checks (equal sides, equal diagonals, right angle).
- Grid independence: n_grid 333 vs 460 give identical counts on smooth curves.

### Roughness sweep (Weierstrass-type radial curve, base 3, K=7, seed 7)
r(theta) = 1 + 0.5·norm·sum lam^k cos(3^k theta + phi_k), lam = 3^(-h).

| h | squares found | min side/diam | median side/diam |
|-----|----|--------|--------|
| 1.0 | 4 | 0.4504 | 0.4620 |
| 0.8 | 4 | 0.4386 | 0.4536 |
| 0.65 | 6 | 0.4131 | 0.4324 |
| 0.5 | 21 | 0.0125 | 0.4215 |
| 0.4 | 42 | 0.0080 | 0.3874 |
| 0.3 | 70 | 0.0086 | 0.3820 |

Interpretation: below h≈0.5 a second population of squares appears at the scale of the fine wiggles; the smallest found square sits AT the detection floor (side_floor_rel=0.008), meaning the cascade continues below our resolution. Count explosion 4 -> 70. This is the degeneration mechanism of the open problem made visible: median square stays large, but nothing forbids the relevant continuum square from living arbitrarily deep in the cascade.

### Known caveats / honesty notes
- Counts on rough curves are lower bounds: candidate cap (400) and grid resolution miss basins; parity claims only reliable where grid fully resolves geometry (ellipse, blob, star).
- The "curve" is numerically its 60k-point polyline; squares are exact for the polyline (1e-12), approximate for the underlying truncated analytic curve.
- Truncated Weierstrass (K=7) is still C-infinity; genuine nowhere-differentiability needs K -> infinity, unreachable numerically. We only see the beginning of the cascade.

### Environment notes (sandbox)
- scipy not preinstalled: `pip install scipy --break-system-packages`.
- bash timeout hard cap 45s per call; chunk work, checkpoint to JSON, `python3 -u` for progress visibility.
- numpy 2.2.6, scipy 1.15.3, matplotlib 3.10.9.

### Files
- `inscribed_squares.py`: finder + gallery + chunked sweep (modes: gallery | sweep_part "h1,h2" out.json | sweep_plot "a.json,b.json,...").
- `squares_gallery.png`: 4-curve gallery with found squares.
- `roughness_sweep.png`: sweep panels + degeneration plot.
- `p10.json ... p03.json`: raw sweep data (corners, sides, residuals per h). Regenerable in ~1.5 min.

### Open next steps (queued, not started)
1. Track ONE square continuously as h varies (homotopy continuation on (t1,t2)): does the big-square branch persist, or annihilate with a partner (parity event)?
2. Push side floor down (adaptive/multiscale grid near small-diagonal region) to measure cascade depth vs h; estimate scaling law min_side(h, K).
3. Vary seed: is the h≈0.5 onset seed-dependent or structural? (h=1/2 is Brownian regularity, suspicious coincidence worth checking.)
4. Read Tao 2017 (integration approach, works for curves = union of two Lipschitz graphs with constant < 1) and map his obstruction onto our cascade picture.
5. Interactive draw-a-curve demo (not started).

---

## Session 2026-07-25 (later): seed sweep on degeneration onset

### Question
Is the small-square onset at h≈0.5 (seed 7) structural (h=1/2 = Brownian regularity?) or a seed artifact?

### Setup
Seeds 1-7, h in {0.65, 0.55, 0.5, 0.45, 0.4}, same pipeline (n_grid=520, side_floor_rel=0.008). Onset := first h (descending) with min_side/diam < 0.1. Two populations remain well separated in every run (big ~0.40-0.44, small <0.016). Script: `seed_sweep.py`, data: `seeds.json`, figure: `seed_sweep.png`.

### Results
| seed | onset h |
|------|---------|
| 1, 4, 5 | 0.50 |
| 2, 6, 7 | 0.55 |
| 3 | 0.65 |

### Findings
1. ANSWERED: onset is NOT pinned at h=1/2. It scatters over 0.5-0.65 by seed. The Brownian-regularity coincidence is dead; do not chase it again.
2. Once the small population appears, min_side sits AT the detection floor (0.008-0.016) for every seed. The cascade itself looks universal; only its entry point depends on the specific wiggle phases. By h=0.5 all 7 seeds have it.
3. Anomaly: seed 3 is non-monotone (small squares at h=0.65, none found at 0.55, back at 0.5). Either genuine (tiny squares need specific local wiggle configurations, which shift with h) or a grid miss (counts are lower bounds; absence is not evidence of absence). Follow-up queued: rerun seed 3, h=0.55 at finer grid before citing this cell.

---

## Session 2026-07-25 (later still): seed 3 recheck + branch continuation

### Seed 3 anomaly: RESOLVED, it was a detection miss
Rerun at n_grid=800, side_floor_rel=0.006: seed 3 at h=0.55 has 13 squares including two tiny ones (side/diam 0.0061, 0.0063), BELOW the standard 0.008 floor. So seed 3 is monotone after all. Methodological consequence, now standing policy: absence cells near the detection floor are unreliable; reported per-seed onsets are UPPER bounds (finer detection can only move onsets earlier). The "onset not pinned at h=1/2" conclusion survives (seed 3 onset is 0.65 or earlier regardless).

### Branch continuation (`continuation.py`, seed 3, h stepped 0.70 <-> 0.40 in 0.01 steps)
Warm-start re-polish of every square's (t1,t2) at each h step; death = residual fails certification after 3 perturbed retries. Data: `state_down.json`, `state_up.json`; figure: `branch_continuation.png`.

Results:
- Down-pass (5 branches seeded at h=0.70): 4 die, 1 reaches h=0.40.
- Up-pass (47 branches seeded at h=0.40): 44 die, 3 reach h=0.70.
- Cross-validation: continuation attrition (47 -> 3 upward) matches the static count collapse (47 vs 5) between h=0.4 and 0.7. Two independent methods agree.

Findings:
1. Individual squares are mortal; the population is not. Branch deaths (annihilations) happen throughout the big-square band, but at every h in [0.40, 0.70] at least one big branch is alive. This is the shape of the open problem in miniature: any proof must attach an invariant to the population (parity/degree), not to any individual square, and must forbid the LAST large branch dying in the rough limit.
2. Small-square branches are ephemeral churn, not a smooth cascade. They wander an order of magnitude in side (0.001-0.02 of diameter) over small h windows, die young, and are created/destroyed by local wiggle rearrangements. The degenerate limit is not "one square shrinking to zero" but continuous turnover at ever-smaller scales.
3. Big-band branches keep near-constant side (~0.40-0.47 of diameter) until they abruptly die. Death is not preceded by gradual shrinking. (Caveat: deaths are detection events; 0.01 h steps + 3 retries could lose a branch a more careful homotopy would keep. Partner-matching of annihilating pairs not implemented.)

---

## Session 2026-07-25 (cont.): Tao 2017 read-and-map

Full notes in `tao_2017_notes.md`. Executive version:
- Tao's conserved integral (alternating sum of int y dx over the four vertex trajectories = exact differential of (a^2-b^2)/2) is a FAMILY-level invariant. It independently confirms our continuation finding: nothing protects individual squares; proofs must ride population invariants.
- His Lipschitz < 1 condition is a uniqueness/single-valuedness condition. Outside it, the moving square becomes multivalued and the fourth-vertex trace self-crosses, which is precisely the churn we measured. Our count explosion (4 -> 70) is the numerical face of his multivaluedness.
- The small-square blow-up limit is his periodic cylinder problem, where the homological count is EVEN. Our ephemeral small-square band is the shadow of that regime: parity-neutral churn, existence never certified. Explains why tiny squares always sit at our detection floor.
- Tao's roadmap: prove combinatorial Conjecture 6 (sign patterns of y1+y2+y3 over finite sets) -> area formulation (Conj 4) -> Toeplitz.
- Dead ends per Tao, added to never-repeat: perturbative counterexamples (killed to 6th order by the integral), pure parity in the periodic setting.

New queued experiments from the mapping (testable with existing tooling):
A. Partner-match branch deaths (expect fold/pair annihilation, corners converging).
B. Max-slope-vs-h per seed; count-explosion onset should track slope crossing ~1. Would explain seed-dependent onset.
C. Adversarial search for Conjecture 6 counterexamples at k=9-13 (constraint solver / annealing). Wagner only did ~500 random (7,7,7) instances in 2016.

---

## Session 2026-07-25 (cont.): experiments B (slope threshold) and C (Conjecture 6 search)

### B: slope-threshold hypothesis FALSIFIED
Metric: S = max |r'(theta)|/r(theta) per (seed, h). Figure `slope_onset.png`.
- S at measured onset scatters 10.7-18.5 across seeds (mean 15.5, std 2.8): no common threshold.
- At fixed h=0.55, S does not separate presence/absence of small squares (seed 3: S=16.1 without; seed 2: S=12.2 with).
- All curves in the studied range have S >> 1, so Tao's slope-1 boundary cannot be the onset criterion in this family.
Conclusion: global max steepness is the wrong statistic; small-square formation depends on specific local wiggle configurations. If revisited, use LOCAL statistics (e.g., co-occurrence of steep segments at ~90 deg relative orientation at matched scales), not global maxima. Do not re-test global slope statistics.

### C: adversarial Conjecture 6 search (`conj6_search.py`)
Validation: axiom-(ii) order patterns per sequence = 2, 8, 42 for k=3,5,7 = the open meandric numbers. Confirms (ii) encodes meander structure; reusable as a correctness check.

Method: annealing over (meander pattern triple, values) with axiom (iii) as penalty; every valid configuration's sign-cell is then settled EXACTLY by LP (within a cell all sgn quantities are constant, so the cell is an open cone and the within-cell supremum of F is a linear program; F > 0 in any cell = genuine counterexample).

Results:
- k=(5,5,5): multiple distinct cells settled across seeds; best cell supremum F = -7e-4 with margin eps=1e-4, i.e. sup F -> 0^- as margins shrink. NO counterexample; the inequality F < 0 is empirically TIGHT (boundary of validity approaches F = 0 through degenerating configurations). This sharpens Wagner's 2016 random sampling: per-cell suprema are exact, not sampled.
- k=(7,7,7): random-start annealing reached ZERO (iii)-valid configurations in ~170 restarts even with 3500-step schedules. Axiom (iii) is combinatorially rigid at k=7; local repair from random meander triples fails.

Never repeat: random-start annealing for k >= 7.
Queued instead (next session): structural generators for valid configs, two candidates:
  a. Arc inflation: insert a nested adjacent pair (u, u+delta) into a valid k config (meander insertion preserves (ii); delta small keeps all sgn sums, so (iii) survives) -> ladder k=5 -> 7 -> 9.
  b. Build configs from Tao's curve dictionary directly (three winding-1 curve systems -> y sequences), guaranteeing validity by construction.

### Session file inventory (added today)
`slope_onset.png`, `conj6_search.py`, `conj6_k5_lp.json` (best settled cell), `state_down.json`, `state_up.json`, `branch_continuation.png`, `tao_2017_notes.md`, `seed_sweep.py`, `seeds.json`, `seed_sweep.png`.

---

## Session 2026-07-25 (cont.): arc-inflation ladder (`inflate.py`)

### Method
Grow valid Conjecture 6 configurations instead of searching for them: insert an adjacent
nested pair (w+1e-3, w+2e-3) next to an existing value w. Meander insertion preserves
axiom (ii); copied signs usually preserve (iii); every candidate is verified exactly and
every new sign-cell settled by LP. This bypasses the k>=7 wall that killed random search.

### Results
- Inflation preserves validity extremely robustly. First (unbalanced) run inflated one
  sequence to k=107, settling ~600 cells along the way, all valid.
- Balanced ladder (insert into shortest sequence): reached k=(23,23,21) in 30s,
  ~90 cells settled at EVERY balanced level: (7,7,7), (9,9,9), ..., (21,21,21).
- Neighbor walk from a (7,7,7) config (crossing cell walls, i.e. sampling cells NOT
  on the inflation tree): 6 further cells, best F = -0.001.
- ZERO counterexamples anywhere. Sup F per cell sits just below 0 (margin-limited) at
  every k. Files: `ladder1.json`, `ladder2.json`, `k777_walk.json`, `k777_seed.json`.

### Where this stands historically
Prior art (per Tao's comment thread, 2016): Wagner checked ~500 random instances at
(7,7,7); Tao suggested (5,5,5) "fairly straightforward" numerically. We now have exact
per-cell LP suprema for hundreds of cells at sizes up to k=23 per sequence. Conjecture 6
looks true and TIGHT (sup F -> 0^-) across the entire explored region.

### Structural insight worth pursuing (possibly the real prize)
Arc insertion changes F only marginally, and inflated cells inherit F < 0 from parents.
This suggests an induction strategy for PROVING Conjecture 6:
  (a) settle the finite base case k=(5,5,5) exhaustively (8^3 pattern triples; MILP over
      sign tensors per triple, symmetry-reduced; plausibly a weekend of compute),
  (b) prove F < 0 is preserved under arc insertion (looks like a finite local computation
      on how the LP optimum moves under nesting),
  (c) prove every valid configuration reduces to a small base by arc DELETIONS
      (meander/temperley-lieb structure theory; this is the hard part).
If (c) fails, the non-reducible "prime" configurations are exactly the interesting
objects, and the search should target them.

### Queued next
- Exhaustive k=(5,5,5): enumerate all 512 pattern triples, MILP or LP-per-sign-cell with
  symmetry reduction. Would make the base case a THEOREM (modulo floating point; use
  rational arithmetic for the final certificate).
- Characterize which insertions/deletions preserve validity (toward step (c)).
- Prime-configuration search at k=7: cells whose arc-deletions all break validity.

---

## Session 2026-07-25 (cont.): EXHAUSTIVE k=(5,5,5) base case + margin scaling law

### Base case settled (`exhaustive_k5.py`, data ex_a/b/c.json)
MILP per pattern triple: 15 continuous y vars, 125 sign binaries, big-M linking
(M=4, EPS=1e-5), axiom (iii) as linear equalities on signs, maximize F.
All 120 pattern triples up to sequence-permutation symmetry (= all 512):
solved, none infeasible, max F = -7e-5 < 0 in EVERY triple.

VERDICT: no counterexample to Conjecture 6 exists at k=(5,5,5) with relative
margins >= 1e-5. First exhaustive treatment of the base case we are aware of
(Tao 2016 suggested it was doable; Wagner only sampled randomly).
Caveats for theorem-grade: HiGHS floating point, big-M formulation; a rational
arithmetic certificate pass is queued. Cells are open cones, so margin
normalization loses no generality up to the EPS threshold.

### DISCOVERY: uniform sharp optimum + margin scaling law
1. Every one of the 120 triples has the SAME optimum: F_max = -7 EPS exactly.
   Not just F < 0: a single universal constant across all meander structures.
2. Spot checks at larger k (3 triples each): k=7 gives F_max = -10 EPS,
   k=9 gives F_max = -13 EPS. Perfect fit to:

       sup F = -((3k-1)/2) * EPS,   i.e.   F <= -((k1+k2+k3-1)/2) * margin

   (equal-k form verified at k=5,7,9; mixed-size form is the natural guess,
   UNTESTED - predicts -8 EPS at (5,5,7).)
This is a sharp quantitative strengthening of Conjecture 6: F is not merely
negative, it is bounded by a linear function of the configuration's minimal
margin. If provable, it implies Conjecture 6 outright.

### Path to proof (next session, in order)
1. Extract LP dual multipliers at the optimum for one cell (scipy linprog dual
   or HiGHS duals): the dual is a nonnegative combination of margin constraints
   certifying F <= -7 EPS. Inspect its structure.
2. Check dual uniformity across triples. A pattern-independent dual = a single
   human-readable inequality chain proving the law for all cells at k=5.
3. Induction over k via arc insertion (we already know insertion preserves
   validity): does the dual certificate extend by +3/2 per inserted pair?
   The +3 EPS per k-step strongly suggests each new arc contributes exactly
   3 margin constraints to the chain.
4. Test mixed sizes (5,5,7) to pin the general constant.
5. Rational arithmetic re-verification of the k=5 base (exact certificate).

### Files added
`exhaustive_k5.py`, `ex_a.json`, `ex_b.json`, `ex_c.json`, `ex_probe.json`,
`ladder1.json`, `ladder2.json`, `k777_seed.json`, `k777_walk.json`.

---

## Session 2026-07-25 (cont.): margin law pinned + dual certificates extracted (`dual_cert.py`)

### Mixed-size law: CONFIRMED, 5-for-5
Generalized MILP for unequal sequence lengths. Blind predictions hit exactly:
- (5,5,7): predicted -8e-5, got -8e-5 (3 pattern triples)
- (5,7,9): predicted -10e-5, got -10e-5 (3 pattern triples)
The law stands as:  sup F = -((k1+k2+k3-1)/2) * EPS
Verified: (5,5,5) exhaustive; (7,7,7), (9,9,9), (5,5,7), (5,7,9) spot-checked. No misses.

### Dual certificates: UNIT MULTIPLIERS, INVARIANT COUNT
For each sampled optimal cell (records 0, 7, 44, 87, 119 = pattern triples
[0,0,0], [0,0,7], [1,2,3], [3,3,5], [7,7,7]):
- ALL dual multipliers are EXACTLY 1 (integral certificate).
- Active constraint count is EXACTLY 7 = (k1+k2+k3-1)/2 in every cell.
- Composition varies: 6 ordering + 1 triple-sum, or 4 ordering + 3 triple-sum.
- Box bounds never active: certificates are purely cone-internal.

Identity-pattern certificate, written out (record 0): F literally decomposes as
  F = sum_i [(y_i1 - y_i2) + (y_i3 - y_i4)]  +  (y_15 + y_25 + y_35)
6 ordering descents + 1 triple sum that axiom (iii) forces negative. Human-readable.

### Emerging theorem shape
For every valid cell, F equals a sum of exactly (k1+k2+k3-1)/2 quantities, each
individually forced negative (ordering descents and (iii)-forced signed triple sums),
with unit coefficients. Integrality of ALL observed duals suggests total unimodularity
or a min-max/matching theorem underneath. Arc insertion adds +1 to the count
(consistent with the ladder's observed -1 EPS drift per inserted pair and the +3 per
equal-k step). Induction target: every insertion adds exactly one negative term to
the decomposition.

### Next (in order)
1. Classify all 120 k=5 certificates by composition (n_ord, n_sum) and by which
   triple sums appear; look for the telescoping chain structure (record 7's sums
   telescope through shared indices).
2. Prove the identity-pattern case by hand: show (iii) forces the last-elements
   triple sum negative for increasing sequences. Small, self-contained lemma.
3. Formalize the insertion step: one inserted arc = one new descent term.
4. Rational re-verification pass (certificates are integral, so exact rational
   checking is trivial once extracted: just verify 7 inequalities sum to F).

### Ruled out / never repeat
- Proving the general conjecture in-session. Not a candidate activity.
- Naive KD-only acceptance without scale-aware residual bounds (attempt 1-2).
- Python-loop exact distance inside optimizers (attempt 3).
- Background processes in the sandbox (attempt 4).
- accept_tol looser than ~1e-6·diam: near-miss population starts at 1.7e-5, would contaminate counts.
