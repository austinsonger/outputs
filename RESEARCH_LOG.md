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

### Ruled out / never repeat
- Proving the general conjecture in-session. Not a candidate activity.
- Naive KD-only acceptance without scale-aware residual bounds (attempt 1-2).
- Python-loop exact distance inside optimizers (attempt 3).
- Background processes in the sandbox (attempt 4).
- accept_tol looser than ~1e-6·diam: near-miss population starts at 1.7e-5, would contaminate counts.
