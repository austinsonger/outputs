# Tao 2017 (arXiv:1611.07441) mapped onto our experiments

Sources: arXiv abstract + Tao's blog post (Nov 2016). Blog post covers the mechanism;
symplectic/Lagrangian framing is NOT in it (that is Greene-Lobb's later toolkit).

## The mechanism in four lines

1. Follow a one-parameter family of squares with vertices gamma_1..gamma_4(t). Then
   int_{g1} y dx - int_{g2} y dx + int_{g3} y dx - int_{g4} y dx = int d[(a^2 - b^2)/2],
   where (a,b) is the square's edge vector. The alternating sum is an EXACT differential:
   it vanishes whenever the family starts and ends degenerate (a=b=0).
2. Theorem: if the curve is two Lipschitz graphs f,g with constant < 1 agreeing at endpoints,
   a square exists. Proof: rotating g by 90 deg about gamma_1(t) hits f in a UNIQUE point
   (slope < 1), so the square family is well defined; the fourth-vertex trace gamma_3 is
   SIMPLE; if it avoided g, Jordan + Stokes + the conserved integral force the enclosed
   region to have zero area. Contradiction, so gamma_3 hits g, which IS an inscribed square.
3. Breaks at Lipschitz >= 1 in exactly two places: (a) the rotated graph can hit f in
   multiple points, so the moving square becomes multivalued; (b) gamma_3 can self-cross,
   Jordan fails, and signed areas cancel across crossings so "zero area" no longer means
   "no region". The integral identity itself stays true; the topology quits.
4. Blow-up analysis of the small-square regime gives a PERIODIC version on the cylinder
   (still open even for polygons). Crucially, the homological count there is EVEN, not odd:
   parity gives nothing. Reduction chain: periodic problem -> area formulation (Conj 4)
   -> three-curve special case (Conj 5) -> purely combinatorial sign-pattern statement
   (Conj 6) about triples y1+y2+y3 from finite real sets. Tao's roadmap: prove 6, then 4,
   then Toeplitz.

## Mapping to our numerical findings

| Tao | Our experiment |
|-----|----------------|
| Slope < 1: square family unique, single-valued | Smooth cases: exactly 1 square (blob), 1 (ellipse), low counts at h=1.0. Our rough curves have max slope ~ 9+ at h=0.5, far outside the regime |
| Rotated graph hits f in multiple points when slope >= 1 | Count explosion 4 -> 70 as h drops: each extra intersection is another square family |
| gamma_3 self-crossing kills the Jordan argument | The churn zone: small-square branches born/dying over narrow h windows. The obstruction is not that squares vanish; it is that no single-valued family survives to carry the invariant |
| Conserved integral is attached to FAMILIES, vanishes only against degenerate endpoints | Matches our continuation result exactly: individual branches are mortal (abrupt deaths, no warning shrink), only population-level invariants persist |
| Periodic (blow-up) problem has EVEN parity: homology cannot certify existence | Our small-square band IS the numerical shadow of the blow-up regime. Even parity = tiny squares appear and die in pairs, guaranteed nothing. Consistent with 44/47 branch deaths and with tiny squares always sitting at the detection floor |
| Homological odd count in the classical setting | Our fully resolved counts were all odd (1, 1, 5, 29) |

## Testable predictions we can check with existing tooling

1. Pair annihilation: small-band branch deaths should occur in PAIRS at (nearly) the same h,
   with the two squares' corners converging before death (fold bifurcation). Our tracker
   records deaths but does not match partners. Cheap upgrade.
2. Slope threshold: compute max effective slope of our Weierstrass curves vs h; the count
   explosion onset should track where local slope crosses ~1 (Tao's uniqueness boundary),
   seed by seed. This would explain WHY onset varies by seed (local slope is phase-dependent).
3. Conjecture 6 stress test: it is a finite combinatorial statement. Adam Wagner checked
   ~500 random instances at k=(7,7,7) in 2016. A structured search (constraint solver /
   MILP over sign patterns, or simulated annealing on the y values) at k=9-13 is well within
   our compute. Finding a counterexample would be a real event; failing to find one under
   adversarial search adds evidence and might suggest an induction pattern.

## Do not repeat / dead ends noted by Tao himself

- Perturbative counterexample constructions: he expanded to 6th order; the conserved
  integral kills all of them. Do not attempt perturbative counterexamples.
- Pure homology/parity in the periodic setting: parity is even; insufficient by itself.
- The blog post has no symplectic content; for that, read Greene-Lobb 2020 directly.
