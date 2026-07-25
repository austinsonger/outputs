# Inscribed Square Problem: Numerical Toolkit

Experimental code for the Toeplitz / Square Peg Problem: does every Jordan curve inscribe a square? This repo does not attempt a proof. It builds numerical instruments to make the open case's failure mode visible, and it attacks the finite combinatorial statement (Tao's Conjecture 6) that sits at the front of Tao's 2017 roadmap.

Two independent tracks:

1. **Geometry track.** A certified inscribed-square finder for arbitrary parametrized Jordan curves, plus roughness sweeps and homotopy continuation that track individual squares as a curve is made rougher.
2. **Combinatorics track.** An adversarial search and a structural generator ("arc inflation") for counterexamples to Tao's Conjecture 6, with every sign-cell's supremum settled exactly by linear programming.

Full experimental history, including every dead end, lives in [RESEARCH_LOG.md](RESEARCH_LOG.md). Read it before rerunning anything.

## Table of Contents

- [Results at a glance](#results-at-a-glance)
- [Prerequisites](#prerequisites)
- [Getting started](#getting-started)
- [Repository layout](#repository-layout)
- [Command reference](#command-reference)
- [How the square finder works](#how-the-square-finder-works)
- [How the Conjecture 6 search works](#how-the-conjecture-6-search-works)
- [Data formats](#data-formats)
- [Python API](#python-api)
- [Reproducing every figure](#reproducing-every-figure)
- [Caveats and honesty notes](#caveats-and-honesty-notes)
- [Known dead ends](#known-dead-ends)
- [Open next steps](#open-next-steps)

## Results at a glance

| Finding | Where |
|---|---|
| Finder validated: ellipse gives exactly 1 square matching the analytic value `s = ab/sqrt(a²+b²)` to 3e-12; star polygon gives 5 congruent squares related by 72° rotation | `squares_gallery.png` |
| Square count explodes 4 → 70 as roughness increases (Weierstrass exponent `h` from 1.0 to 0.3); a second population appears at the scale of the fine wiggles | `roughness_sweep.png`, `p10.json`…`p03.json` |
| The degeneration onset is NOT pinned at `h = 1/2`. It scatters over 0.50 to 0.65 across seeds, killing the Brownian-regularity coincidence | `seed_sweep.png`, `seeds.json` |
| Individual squares are mortal, the population is not. Big branches keep near-constant side then die abruptly; at every `h` in [0.40, 0.70] at least one big branch is alive | `branch_continuation.png`, `state_up.json`, `state_down.json` |
| Global max-slope threshold hypothesis falsified: steepness at onset scatters 10.7 to 18.5 with no common threshold | `slope_onset.png` |
| Axiom (ii) order-pattern counts come out as 2, 8, 42 for k = 3, 5, 7, exactly the open meandric numbers (self-validation of the encoding) | `conj6_search.py patterns k` |
| Zero counterexamples to Conjecture 6 anywhere. Exact per-cell LP suprema for hundreds of cells up to k = 23 per sequence, all with `sup F → 0⁻`. The inequality is empirically true and tight | `conj6_k5_lp.json`, `ladder1.json`, `ladder2.json`, `k777_walk.json` |
| **Base case settled exhaustively.** All 512 pattern triples at `k=(5,5,5)` (120 up to symmetry) solved by MILP: none infeasible, every optimum negative. No counterexample with relative margins ≥ 1e-5 | `exhaustive_k5.py`, `ex_a.json`, `ex_b.json`, `ex_c.json` |
| **Margin scaling law (new, conjectural).** Every pattern triple has the SAME optimum, `F_max = -7·EPS` at k = 5, and spot checks give `-10·EPS` at k = 7 and `-13·EPS` at k = 9: a perfect three-point fit to `sup F = -((3k-1)/2)·EPS`, i.e. `F ≤ -((k1+k2+k3-1)/2)·margin`. A sharp quantitative strengthening of Conjecture 6; provable versions would imply it outright | `ex_*.json`, RESEARCH_LOG session entries |

Prior art for comparison: Wagner checked roughly 500 random `(7,7,7)` instances in 2016, and Tao called `(5,5,5)` "fairly straightforward" numerically. This repo settles `(5,5,5)` exhaustively and reaches k = 23 with exact cell suprema via arc inflation.

## Prerequisites

- Python 3.10+
- numpy, scipy, matplotlib

Verified versions: numpy 2.2.6, scipy 1.15.3, matplotlib 3.10.9.

## Getting started

```bash
git clone <repo-url>
cd outputs
```

Install dependencies:

```bash
pip install numpy scipy matplotlib
```

If you are on a Debian-style system with an externally managed Python (the sandbox this was developed in), scipy is not preinstalled and pip refuses to touch the system environment. Either use a venv:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install numpy scipy matplotlib
```

or override:

```bash
pip install scipy --break-system-packages
```

Smoke test, roughly 30 seconds, writes `squares_gallery.png`:

```bash
python3 -u inscribed_squares.py gallery .
```

Expected console output: 4 curves, with the ellipse reporting exactly 1 square and the 5-point star reporting 5.

Second smoke test, instant, validates the Conjecture 6 encoding:

```bash
python3 -u conj6_search.py patterns 7
# k=7: 42 axiom-(ii) patterns
```

42 is the third open meandric number. If you get anything else, the axiom (ii) implementation is broken.

## Repository layout

```
.
├── RESEARCH_LOG.md          # Append-only experiment history. Source of truth.
├── tao_2017_notes.md        # Notes on Tao 2017 + the Conjecture 6 roadmap
│
├── inscribed_squares.py     # Core: curves, Curve class, find_squares, gallery, roughness sweep
├── seed_sweep.py            # Seed-robustness of the degeneration onset
├── continuation.py          # Homotopy continuation of individual square branches over h
├── conj6_search.py          # Conjecture 6: axioms, meander patterns, annealing, per-cell LP
├── inflate.py               # Conjecture 6: arc-inflation ladder + neighbor-cell walk
├── exhaustive_k5.py         # Conjecture 6: exhaustive k=(5,5,5) base case via MILP
│
├── squares_gallery.png      # 4-curve gallery with found squares
├── roughness_sweep.png      # Sweep panels + degeneration plot
├── seed_sweep.png           # Onset vs seed
├── slope_onset.png          # Max-slope-vs-onset (falsified hypothesis)
├── branch_continuation.png  # Branch lifetimes over h
│
├── p10.json p08.json p065.json p05.json p04.json p03.json   # Roughness sweep raw data
├── seeds.json               # Seed sweep raw data (35 rows)
├── state_down.json          # Continuation state, h 0.70 -> 0.40
├── state_up.json            # Continuation state, h 0.40 -> 0.70
├── conj6_k5*.json           # Annealing / LP results at k=(5,5,5)
├── ladder1.json ladder2.json  # Arc-inflation ladders (unbalanced to k=107, balanced to k=23)
├── k777_seed.json k777_walk.json  # (7,7,7) seed config + neighbor-cell walk
└── ex_a.json ex_b.json ex_c.json ex_probe.json  # Exhaustive k=(5,5,5) MILP results (all 120 triples)
```

Every script is a plain CLI with `mode` as `argv[1]`. There is no package, no build step, and no config file. Run scripts from the repo root, since `seed_sweep.py`, `continuation.py`, and `inflate.py` import from their siblings.

## Command reference

Use `python3 -u` for all long runs so progress is visible when output is piped or captured.

### inscribed_squares.py

| Command | Description |
|---|---|
| `python3 -u inscribed_squares.py gallery [OUTDIR]` | Find squares on ellipse, smooth blob, star polygon, and a Weierstrass curve. Writes `squares_gallery.png` to `OUTDIR` (default `.`). |
| `python3 -u inscribed_squares.py sweep_part "H1,H2,..." OUT.json` | Run the roughness sweep for the listed `h` values, write raw square data to `OUT.json`. Chunk this: one or two `h` values per call. |
| `python3 -u inscribed_squares.py sweep_plot "a.json,b.json,..." [OUTDIR]` | Combine sweep part files into `roughness_sweep.png`. |

Full sweep reproduction, roughly 1.5 minutes total:

```bash
python3 -u inscribed_squares.py sweep_part "1.0"  p10.json
python3 -u inscribed_squares.py sweep_part "0.8"  p08.json
python3 -u inscribed_squares.py sweep_part "0.65" p065.json
python3 -u inscribed_squares.py sweep_part "0.5"  p05.json
python3 -u inscribed_squares.py sweep_part "0.4"  p04.json
python3 -u inscribed_squares.py sweep_part "0.3"  p03.json
python3 -u inscribed_squares.py sweep_plot "p10.json,p08.json,p065.json,p05.json,p04.json,p03.json" .
```

### seed_sweep.py

| Command | Description |
|---|---|
| `python3 -u seed_sweep.py run SEED "H1,H2,..." OUT.json` | Run the finder for one seed across several `h`. Appends to `OUT.json` if it exists, so it is safe to call repeatedly. |
| `python3 -u seed_sweep.py summarize "a.json,..." [OUTDIR]` | Build `seed_sweep.png` and print the onset table. |

Reproduce `seeds.json` (35 rows, 7 seeds x 5 h values). Run one seed per call to stay inside short timeouts:

```bash
for s in 1 2 3 4 5 6 7; do
  python3 -u seed_sweep.py run $s "0.65,0.55,0.5,0.45,0.4" seeds.json
done
python3 -u seed_sweep.py summarize seeds.json .
```

Settings used for the published numbers: `n_grid=520`, `side_floor_rel=0.008`. "Small-square population present" is defined as `min_side/diam < 0.1`; the two populations are cleanly separated (big ~0.40 to 0.46, small < 0.016), so the threshold choice is not load-bearing.

### continuation.py

| Command | Description |
|---|---|
| `python3 -u continuation.py seedpass SEED H_START H_END N_STEPS state.json` | Seed branches at `H_START` with a full static find, then walk `N_STEPS` toward `H_END`, checkpointing to `state.json`. |
| `python3 -u continuation.py steps state.json N` | Continue an existing state by `N` more steps. Resumable; this is how long runs are chunked. |
| `python3 -u continuation.py plot "state_down.json,state_up.json" OUT.png` | Draw branch lifetimes vs `h`. |

Reproduce the published passes (seed 3, `h` stepped in 0.01 increments):

```bash
python3 -u continuation.py seedpass 3 0.70 0.40 5 state_down.json
python3 -u continuation.py steps state_down.json 25   # repeat until h_index reaches the end
python3 -u continuation.py seedpass 3 0.40 0.70 5 state_up.json
python3 -u continuation.py steps state_up.json 25
python3 -u continuation.py plot "state_down.json,state_up.json" branch_continuation.png
```

A branch is declared dead when its warm-started re-polish fails certification after 3 perturbed retries. Death is a detection event, not a proof of annihilation. See caveats.

### conj6_search.py

| Command | Description |
|---|---|
| `python3 -u conj6_search.py patterns K` | Enumerate and count axiom-(ii)-valid order patterns for sequence length `K`. Should print the open meandric numbers: 2, 8, 42 for k = 3, 5, 7. |
| `python3 -u conj6_search.py search K SECONDS OUT.json [SEED]` | Anneal over (pattern triple, values) maximizing `F` under axioms (i)-(iii), then settle every valid configuration's sign-cell exactly by LP. Writes best record to `OUT.json`. |

```bash
python3 -u conj6_search.py patterns 5           # k=5: 8 axiom-(ii) patterns
python3 -u conj6_search.py search 5 40 conj6_k5.json 1
```

Do not run `search 7` from random starts. 170 restarts with 3500-step schedules produced zero axiom-(iii)-valid configurations. Axiom (iii) is combinatorially rigid at k = 7. Use `inflate.py` instead.

### inflate.py

| Command | Description |
|---|---|
| `python3 -u inflate.py ladder "IN1.json,IN2.json" OUT.json BUDGET_S` | Grow valid configurations by inserting nested adjacent pairs, climbing k = 5 → 7 → 9 → …. Each input file supplies `{"ys": [[...],[...],[...]]}`. Writes `{"cells": [...], "best": record}`. |
| `python3 -u inflate.py walk IN.json OUT.json BUDGET_S` | Perturb a valid configuration across sign-cell walls to sample cells that are not on the inflation tree. |

```bash
python3 -u inflate.py ladder "conj6_k5_lp.json,conj6_k5_s1.json" ladder2.json 30
python3 -u inflate.py walk k777_seed.json k777_walk.json 20
```

The balanced ladder (always insert into the currently shortest sequence) reached `k=(23,23,21)` in 30 seconds, settling roughly 90 cells at every balanced level from `(7,7,7)` to `(21,21,21)`.

### exhaustive_k5.py

| Command | Description |
|---|---|
| `python3 -u exhaustive_k5.py run START END OUT.json` | Solve pattern triples `START..END-1` (of 120 up to symmetry) as MILPs: 15 continuous `y` vars, 125 sign binaries with big-M linking (M=4, EPS=1e-5), axiom (iii) as linear equalities on the sign bits, maximize `F`. Prints per-triple optima; screams if any `F ≥ 0`. |
| `python3 -u exhaustive_k5.py combine "a.json,b.json,..."` | Merge chunk results and print the verdict. |

Published run (three chunks, ~15 s each):

```bash
python3 -u exhaustive_k5.py run 0 45 ex_a.json
python3 -u exhaustive_k5.py run 45 90 ex_b.json
python3 -u exhaustive_k5.py run 90 120 ex_c.json
python3 -u exhaustive_k5.py combine ex_a.json,ex_b.json,ex_c.json
# 120 triples: 120 solved, 0 infeasible, 0 solver issues
# max F = min F = -0.000070  ->  NO counterexample at k=(5,5,5)
```

Larger `k` spot checks: set `exhaustive_k5.K = 7` (or 9) after import and pass patterns from `enumerate_patterns(K)` to `solve_triple`. Observed optima: `-10·EPS` at k = 7, `-13·EPS` at k = 9.

## How the square finder works

### Formulation

A square inscribed in curve `G` corresponds exactly to a parameter pair `(t1, t2)`:

```
p1 = G(t1),  p2 = G(t2)
c  = (p1 + p2) / 2
w  = rot90((p2 - p1) / 2)
q1 = c + w,  q2 = c - w
```

`(p1, q1, p2, q2)` has equal perpendicular diagonals sharing a midpoint, so it is an exact square by construction. The only thing left to verify is that `q1` and `q2` also lie on `G`. That reduces the search to a 2D problem on the `(t1, t2)` torus with objective `d(q1, G)² + d(q2, G)²`.

### Pipeline

1. **Dense sample.** Polyline sample of the curve (24k points smooth, 60k rough) plus a `cKDTree`.
2. **Grid scan.** Fully vectorized scan of the `(t1, t2)` torus at 460 to 520 points per axis. Two bands are masked to infinity: `|t1 - t2| < dt_min` (where `p1 ≈ p2` produces a fake zero) and geometrically tiny diagonals.
3. **Local minima.** 8-neighbor minima with torus wraparound, below `(0.03 · diam)²`, capped at 400 candidates.
4. **Refine.** Nelder-Mead on the exact parametrization, not on the sampled polyline.
5. **Dedupe, then polish.** Dedupe by corner sets at `0.02 · diam` tolerance FIRST, then polish survivors against segment-exact point-to-polyline distance. Order matters: polishing before dedupe means polishing 400 candidates instead of 40 and blows the time budget.
6. **Accept.** `side >= side_floor_rel · diam` AND `residual <= min(accept_tol_rel · diam, 0.01 · side)`.

### Why the acceptance rule looks like that

KD nearest-sample distance has a resolution floor at roughly the sample spacing. It cannot certify anything below that, which is why the final polish uses `dist_exact` (batch einsum segment projection over the 6 nearest samples and their adjacent segments).

The residual distribution is cleanly bimodal: true squares land at ~1e-12, near-misses start at ~1.7e-5. A cutoff of `1e-7 · diam` sits inside a 7-order-of-magnitude gap. Do not loosen `accept_tol_rel` past ~1e-6, or the near-miss population contaminates counts.

The scale-aware half of the rule (`0.01 · side`) is what stops the finder reporting a "square" of side 0.04 with 8% corner error. Diameter-relative tolerance alone is not enough.

### Curves provided

| Constructor | Description |
|---|---|
| `ellipse(a=1.6, b=1.0)` | Convex control case with a known analytic answer. |
| `fourier_blob(seed, n_modes=6, amp=0.16)` | Smooth generic Jordan curve, radial Fourier series with fast decay. |
| `star_polygon(points=5, r_out=1.0, r_in=0.45)` | Non-convex, piecewise linear, corners where they matter. |
| `weierstrass_curve(h, K=7, base=3, amp=0.5, seed, n_dense=60000)` | Roughness dial. `r(θ) = 1 + amp·norm·Σ λ^k cos(base^k θ + φ_k)`, `λ = base^(-h)`. Radial form guarantees the curve is simple. Lower `h` means rougher. |

## How the Conjecture 6 search works

### The statement

Odd `k1, k2, k3`. For each `i`, distinct reals `y_{i,1..k_i}`, with convention `y_{i,0} = y_{i,k_i+1} = -∞`. Axioms:

- **(i)** All triple sums `y_{1,p} + y_{2,q} + y_{3,r} ≠ 0`.
- **(ii)** Non-crossing, for `0 <= p < q <= k_i` of the same parity: `Σ_{a∈{p,p+1}} Σ_{b∈{q,q+1}} (-1)^{a+b} sgn(y_{i,a} - y_{i,b}) = 0`.
- **(iii)** Non-crossing sums, for `0 <= p <= k1, 0 <= q <= k2, 0 <= r <= k3` of the same parity: the alternating sum over the 8 corners of `sgn(y_{1,a} + y_{2,b} + y_{3,c})` vanishes.

Conclusion: `F := Σ_i Σ_p (-1)^{p-1} y_{i,p} < 0`.

A counterexample is any configuration satisfying (i), (ii), (iii) with `F >= 0`.

### Search design

Axiom (ii) constrains only the ORDER pattern of each sequence, which is a meander permutation. Those are enumerated by brute force once per `k`. The count reproducing the open meandric numbers (2, 8, 42) is the correctness check for the whole encoding.

Given a pattern triple and concrete values, all `sgn` quantities are constant inside a sign-cell, so the cell is an open cone and the within-cell supremum of the linear functional `F` is a **linear program**. This is the key move: each cell is settled exactly rather than sampled. `cell_lp_max` runs that LP with margin `eps` (default 1e-4) and box bound (default 1.0).

Two ways to reach valid configurations:

- **Annealing** (`conj6_search.py search`): state is `(pattern_1, pattern_2, pattern_3, values)`, moves perturb values or swap patterns, axiom (iii) is a penalty and a hard filter for records. Works at k = 5. Fails completely at k = 7.
- **Arc inflation** (`inflate.py ladder`): take a valid configuration and insert an adjacent nested pair `(u + 1e-3, u + 2e-3)` next to an existing value `u`. Meander insertion preserves (ii), and because the new values sit next to an existing one, the new triple-sum signs copy the neighbor's, so (iii) usually survives. Every candidate is verified against the exact axioms anyway, and every survivor's cell is settled by LP. This bypasses the k >= 7 wall entirely.

### The structural opening

Arc insertion changes `F` only marginally, and inflated cells inherit `F < 0` from their parents. That suggests an induction proof strategy:

1. ~~Settle the base case `k=(5,5,5)` exhaustively.~~ **DONE** (`exhaustive_k5.py`): all 120 symmetry-reduced pattern triples solved, every optimum negative, none infeasible.
2. Prove `F < 0` is preserved under arc insertion (looks like a finite local computation on how the LP optimum moves under nesting).
3. Prove every valid configuration reduces to a small base by arc DELETIONS (meander / Temperley-Lieb structure theory; this is the hard part).

If step 3 fails, the non-reducible "prime" configurations are exactly the interesting objects, and search should target them.

### The margin scaling law (conjectural, strongest finding in the repo)

The exhaustive run did not just verify the base case; it exposed structure. Every one of the 120 pattern triples has the SAME MILP optimum, `F_max = -7·EPS` exactly. Spot checks at three pattern triples each: k = 7 gives `-10·EPS`, k = 9 gives `-13·EPS` (the k = 9 value was predicted before it was computed). Three exact hits on:

```
sup F = -((3k-1)/2) · EPS        (equal sequence lengths k)
F ≤ -((k1+k2+k3-1)/2) · margin   (natural mixed-size form, UNTESTED)
```

where `margin` is the minimal strict quantity (ordering gaps and `|triple sums|`). This is sharper than Conjecture 6's `F < 0`, and the `+3` per `k`-step suggests each inserted arc contributes exactly three constraints to the underlying inequality chain, dovetailing with the arc-insertion induction. The LP dual multipliers at any optimum literally write down a certificate for that cell; a pattern-independent dual would be a human-readable proof at k = 5. See RESEARCH_LOG "path to proof" for the ordered attack.

## Data formats

### Roughness sweep parts (`p10.json` … `p03.json`)

List of records, one per `h`:

```json
[{"h": 0.5, "diam": 2.83, "squares": [{"corners": [[x,y],...], "side": 1.19, "err": 3.1e-13}, ...]}]
```

### Seed sweep (`seeds.json`)

Flat list, 35 rows:

```json
[{"seed": 3, "h": 0.55, "count": 13, "min_side": 0.0061, "med_side": 0.4318}]
```

`min_side` and `med_side` are relative to curve diameter. `min_side` is `null` when no squares passed the floor.

### Continuation state (`state_down.json`, `state_up.json`)

```json
{"seed": 3, "h_schedule": [0.70, 0.69, ...], "h_index": 30, "branches": [...]}
```

Resumable by design: `continuation.py steps` reads `h_index`, advances, and rewrites the file.

### Conjecture 6 records

`conj6_k5_lp.json` is a single settled cell:

```json
{"F": -0.0007, "ys": [[...],[...],[...]], "pats": [[...],[...],[...]]}
```

`ladder*.json` and `k777_walk.json` wrap the same record shape:

```json
// ladder: summary is keyed by (k1,k2,k3), one entry per rung
{"best": {...record...},
 "summary": {"(7, 7, 7)":   {"cells": 70, "supF": -0.0010},
             "(21, 21, 21)": {"cells": 90, "supF": -0.0031},
             "(23, 23, 21)": {"cells": 86, "supF": -0.0033}}}

// walk
{"best": {...record...}, "n_cells": 6}
```

Note the pattern in `supF`: it drifts down by exactly one `eps` per inserted pair. That is the margin term, not a real trend. Every rung is equally tight against 0.

`F` is the LP supremum over the cell, not the value at the sampled point. Negative `F` means that entire cell is counterexample-free.

### Exhaustive base case (`ex_a.json`, `ex_b.json`, `ex_c.json`)

List of per-triple MILP results:

```json
[{"tri": [0, 0, 1], "status": "ok", "F": -7e-05, "y": [ ...15 values... ]}]
```

`tri` indexes into the symmetry-reduced pattern-triple list (`triple_list()` in `exhaustive_k5.py`); `status` is `ok`, `infeasible`, or `statusN` (solver issue, retry). `F` is the maximum of the objective over ALL sign-cells with that pattern triple, so a negative value clears the whole triple at once.

## Python API

Import directly from `inscribed_squares`:

```python
from inscribed_squares import ellipse, weierstrass_curve, find_squares, square_corners

curve = ellipse(a=1.6, b=1.0)
squares = find_squares(curve, n_grid=460, side_floor_rel=0.01, verbose=True)
for s in squares:
    print(s["side"] / curve.diam, s["err"], s["corners"])
```

### `find_squares(curve, n_grid=460, dt_min=0.004, side_floor_rel=0.01, accept_tol_rel=1e-7, max_refine=400, verbose=False)`

| Parameter | Meaning | Notes |
|---|---|---|
| `n_grid` | Grid points per torus axis | 460 default, 520 for seed sweeps, 800 for near-floor rechecks. Counts are lower bounds at any resolution. |
| `dt_min` | Mask width around `t1 ≈ t2` | Suppresses the degenerate fake zero at `p1 = p2`. |
| `side_floor_rel` | Reject squares below this fraction of diameter | This is the numerical degeneracy floor. It is exactly the leak that the open problem's limiting argument suffers from. |
| `accept_tol_rel` | Residual cap, relative to diameter | Do not exceed ~1e-6. Effective tolerance is `min(accept_tol_rel·diam, 0.01·side)`. |
| `max_refine` | Candidate cap before refinement | 400. Raising it costs time roughly linearly. |

Returns a list of dicts with `corners` (4x2 array in diagonal order `p1, q1, p2, q2`), `side`, and `err` (max corner distance to the curve).

### `Curve`

| Member | Purpose |
|---|---|
| `Curve(fn, name, n_dense=24000)` | Wraps a vectorized `t ∈ [0,1) → (n,2)` map. |
| `.diam` | Bounding-box diagonal, the length scale everything is relative to. |
| `.spacing` | Mean segment length, the KD distance resolution floor. |
| `.point(t)` | Exact evaluation, wraps `t` mod 1. |
| `.dist_fast(pts)` | KD nearest-sample distance. Fast, slight overestimate, floored at `.spacing`. |
| `.dist_exact(pts)` | Segment-projected point-to-polyline distance, vectorized over a batch. Use this for certification. |

To use your own curve, pass any vectorized parametrization:

```python
import numpy as np
from inscribed_squares import Curve, find_squares

def fn(t):
    th = 2*np.pi*np.asarray(t)
    r = 1 + 0.3*np.cos(3*th)
    return np.column_stack([r*np.cos(th), r*np.sin(th)])

squares = find_squares(Curve(fn, "trefoil-ish blob", n_dense=40000))
```

The parametrization must be vectorized over an array of `t` and must trace a simple closed curve. Radial forms `r(θ) > 0` are simple automatically.

## Reproducing every figure

| Figure | Command |
|---|---|
| `squares_gallery.png` | `python3 -u inscribed_squares.py gallery .` |
| `roughness_sweep.png` | The six `sweep_part` calls, then `sweep_plot` (see above). ~1.5 min. |
| `seed_sweep.png` | The seven `seed_sweep.py run` calls, then `summarize`. |
| `branch_continuation.png` | `continuation.py seedpass` down and up, then `plot`. |
| `slope_onset.png` | Generated during the experiment B analysis of `seeds.json`. Hypothesis was falsified; see RESEARCH_LOG. |

All plotting uses the `Agg` backend, so everything works headless.

## Caveats and honesty notes

Read these before citing any number from this repo.

- **Counts on rough curves are lower bounds.** The candidate cap (400) and grid resolution miss basins. Parity claims are only reliable where the grid fully resolves the geometry: ellipse, blob, star.
- **Absence near the detection floor is not evidence of absence.** This is standing policy after the seed 3 anomaly, where a "missing" cell at `h = 0.55` turned out to hold two squares at side/diam 0.0061 and 0.0063, just below the 0.008 floor. Reported per-seed onsets are therefore UPPER bounds; finer detection can only move onsets earlier.
- **The "curve" is numerically its polyline.** Squares are exact for the 60k-point polyline (1e-12) and approximate for the underlying analytic curve.
- **Truncated Weierstrass (K = 7) is still C-infinity.** Genuine nowhere-differentiability needs K → ∞, which is unreachable numerically. What is visible is the beginning of the cascade, not the limit.
- **Branch deaths are detection events.** 0.01 `h` steps plus 3 perturbed retries could lose a branch that a more careful homotopy would keep. Partner-matching of annihilating pairs is not implemented.
- **Conjecture 6 results are floating point.** `sup F` values sit just below 0 and are margin-limited (`eps = 1e-4` for LP cells, `EPS = 1e-5` for the exhaustive MILPs). A publishable certificate needs rational arithmetic.
- **The exhaustive base case has two additional caveats.** (1) Big-M formulation (M = 4) with HiGHS floating point; a counterexample living entirely within margins < 1e-5 after normalization to the unit box would be missed. Cells are open cones, so margin normalization is legitimate, but the EPS threshold is a real (small) hole until the rational pass exists. (2) The margin scaling law rests on exhaustive k = 5 plus THREE spot-checked pattern triples each at k = 7 and k = 9. Treat the mixed-size form as a guess until (5,5,7) is run (predicts `-8·EPS`).

## Known dead ends

Do not repeat these. Each cost real time.

| Attempt | Failure | Lesson |
|---|---|---|
| Acceptance tolerance `1.5e-3·diam`, KD distance only | Accepted tiny "squares" (side 0.04) with 8% corner error | Acceptance must be scale-aware: cap at a fraction of SIDE, not just diameter |
| Nelder-Mead restarts with scale-aware filter, still KD point-cloud distance | Residuals stuck at ~2e-3 on big squares | KD nearest-sample distance has a resolution floor at sample spacing and cannot certify below it |
| Segment-exact distance polish inside the main NM loop, per-point Python loops | Timeout, over 42 s for 400 candidates | Vectorize the distance (batch einsum segment projection) and polish AFTER dedupe |
| `nohup ... &` background runs in the sandbox | Process killed at the call boundary, log never written | Sandbox bash calls are independent. Chunk long runs, checkpoint to JSON, use `python3 -u`, stay under ~42 s per call |
| Global max-slope `S = max|r'(θ)|/r(θ)` as the onset criterion | `S` at onset scatters 10.7 to 18.5 across seeds; at fixed `h` it does not separate presence from absence | Global steepness is the wrong statistic. If revisited, use LOCAL statistics such as co-occurrence of steep segments at ~90° relative orientation at matched scales |
| Random-start annealing for Conjecture 6 at k >= 7 | Zero valid configurations in ~170 restarts with 3500-step schedules | Axiom (iii) is combinatorially rigid at k = 7. Grow configurations with `inflate.py`, do not search for them |
| Perturbative counterexamples to Conjecture 6 | Killed to 6th order by Tao's conserved integral | Per Tao 2017 |
| Pure parity arguments in the periodic setting | Homological count is EVEN there, so parity certifies nothing | Per Tao 2017; this is why tiny squares always sit at the detection floor |

Also ruled out as an activity: attempting a proof of the general conjecture in-session.

## Open next steps

Geometry track:

1. Partner-match branch deaths to confirm fold / pair annihilation with corners converging.
2. Push the side floor down with an adaptive multiscale grid near the small-diagonal region, and estimate a scaling law `min_side(h, K)`.

Combinatorics track (higher value), in attack order:

3. Extract LP dual multipliers at the optimum for one cell: the dual is a nonnegative combination of margin constraints certifying `F ≤ -7·EPS`, i.e. a machine-written proof for that cell. Inspect its structure.
4. Check dual uniformity across the 120 pattern triples. A pattern-independent dual is a single human-readable inequality chain proving the margin law for all of k = 5.
5. Induction over k via arc insertion: does the certificate extend by `+3/2` per inserted pair? The uniform `+3·EPS` per k-step strongly suggests each new arc contributes exactly three margin constraints.
6. Test mixed sizes `(5,5,7)` to pin the general constant (predicts `-8·EPS`).
7. Rational-arithmetic re-verification of the k = 5 base (exact certificate; removes the EPS hole).
8. Characterize which insertions and deletions preserve validity, toward step 3 of the induction strategy; prime-configuration search at k = 7 (cells whose arc-deletions all break validity).

## Context: what is already known about the problem

Not to be re-derived.

- Toeplitz (1911) conjectured every Jordan curve inscribes a square. Still open for general continuous curves.
- Proven cases: convex (Toeplitz), analytic (Jerrard 1961), locally monotone (Stromquist 1989), smooth via every-aspect-ratio rectangles (Greene-Lobb 2020, symplectic).
- Rectangle existence holds for ALL Jordan curves via Vaughan's Möbius-band argument.
- The entire difficulty: on rough curves, inscribed squares in approximating curves can degenerate (side → 0) in the limit. Any proof must keep one square uniformly large.
- Tao's roadmap: prove combinatorial Conjecture 6 (sign patterns of `y1+y2+y3` over finite sets) → area formulation (Conjecture 4) → Toeplitz. See [tao_2017_notes.md](tao_2017_notes.md).

What this repo adds to that picture: the continuation experiment shows individual squares are mortal while the population survives, which independently confirms Tao's family-level conserved integral. Any proof must attach its invariant to the population (parity or degree), and must forbid the LAST large branch from dying in the rough limit.

On the combinatorial front, the repo settles the `k=(5,5,5)` base case of Tao's Conjecture 6 exhaustively (all 512 pattern triples, MILP, no counterexample) and observes the margin scaling law `sup F = -((3k-1)/2)·EPS`, exact at k = 5, 7, 9, a sharp quantitative strengthening of the conjecture with a concrete certificate-based path to a proof.
