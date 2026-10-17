# PROTOCOL — 12-model-small-x (frozen before any data; not edited afterwards)

Date: 2026-10-15. Campaign: pair-46 factor-base campaign, item 12.

## Question
Does choosing the "small coordinate" factor base in a non-Weierstrass curve model
(Montgomery u, twisted Edwards y and x, Hessian x and x+y) change the m = 2
decomposition behaviour of the factor base, compared with Weierstrass small-x and with
random factor bases of the same size, on the same curve and group?
Prediction (task): indistinguishable.

## Part A — pair 46 (exact)
Statement to check exactly: Montgomery, twisted Edwards (and Jacobi quartic / Huff)
models over F_p need a rational point of order 2 (Montgomery and Edwards need 4 | #E);
Hessian / twisted Hessian / generalized Hessian models need a rational point of order 3.
For pair 46 (crater and floor, order N; their quadratic twists, order 5 v^2) check
N mod 2, N mod 3, (5 v^2) mod 2, (5 v^2) mod 3, and that x^3 + a x + b has no root in
F_p (gcd(x^p - x, f) = 1) for both curves. If all hold, no such model exists over F_p for
pair 46 or its twists, and the experiment can only run on toy curves where the models
exist (Part B). This is an exact obstruction, reported as such.

## Part B — toy experiment (comparative)
### Populations (curves where the models exist)
- Pop M: random Montgomery curves B v^2 = u^3 + A u^2 + u over F_p (A, B uniform,
  B != 0, A^2 != 4), with the isomorphic short Weierstrass model X = (u + A/3)/B,
  Y = v/B, and the twisted Edwards model y_E = (u-1)/(u+1), x_E = u/v.
- Pop H: random Hessian curves x^3 + y^3 + 1 = 3 d x y (d uniform, d != 0, d^3 != 1).
  All affine points are enumerated natively (lines through (0,-1)); a short Weierstrass
  model is obtained from u0 = 1/(x+y+d) (double pole at O = (1:-1:0)) and
  v0 = (y-x)/(x+y+d), with the cubic relation v0^2 = c3 u0^3 + c2 u0^2 + c1 u0 + c0 fitted
  and then VERIFIED on every point; #H(F_p) must equal the Weierstrass point count.
  Any failed verification aborts the curve (recorded).
- Primes: largest prime below 2^b for b in {16, 20, 22}.
- Curves per (pop, size): 60 (b=16), 40 (b=20), 20 (b=22). Seeds: splitmix64 from
  base seed 46012 and (pop, b, curve index, arm). Factor-base size K = 2*round(sqrt(p)/2)
  (256, 1024, 2048).

### Arms (every set has exactly K points; all arithmetic done in the Weierstrass model)
Test arms ("take coordinate values c = 0, 1, 2, ... as integers in [0, p), include every
point with that coordinate, truncate the last fibre to reach exactly K"):
- Pop M: W-x (Weierstrass X), M-u (Montgomery u), E-y (Edwards y), E-x (Edwards x).
- Pop H: W-x, H-x (Hessian x), H-s (Hessian x+y).
Controls on every curve:
- R-sym: random negation-closed set (random X values with two points, both points).
- R-sym2: an independent second R-sym draw (null control: R-sym2 vs R-sym).
- R-asym: K distinct uniformly random points.
- PC1 (positive control, strong): {+-kG : 1 <= k <= K/2}, G a random point of order > K.
- PC2 (positive control, weak): {+-kG : 1 <= k <= K/4} union K/2 points of R-sym type.
Matched control for the primary comparison: negation-closed arms (W-x, M-u, E-y, H-s,
PC1, PC2, R-sym2) vs R-sym; the others (E-x, H-x) vs R-asym.

### Instrument (identical on every arm)
One Rust binary in this directory. For each set F: compute all K(K+1)/2 unordered sums
P_i + P_j (i <= j) in affine Weierstrass arithmetic and mark them in a bitmap of the
group (key 2X + [Y > p/2], identity separate).
- Primary metric: coverage C = |F + F| / N (fraction of the group decomposable as a sum of
  two factor-base points, i.e. m = 2 decomposition success probability for a uniform target).
- Secondary: collision excess = (K(K+1)/2 - |F+F|) / (K(K+1)/2).
- Semaev cross-check (on negation-closed arms, after mapping back to Weierstrass X):
  for 200 random non-identity targets R per (curve, arm), decide decomposability by
  solving S3(x1, x2, X(R)) = 0 for x2 over all x1 in X(F) and testing x2 in X(F); record
  the number of disagreements with the bitmap (expected 0).
- Exact identity (E1): the mean over all targets of the number of ordered representations
  equals K^2 / N for every set; so only the spread (coverage/collisions) can differ.

### Effect size and decision rule
Per curve: rel = C_arm / C_matched_control - 1. Report mean rel over curves with a
percentile bootstrap CI (20000 resamples, seed fixed), at 95% and at the Bonferroni level
1 - 0.05/21 (7 test arms x 3 sizes). Equivalence margin delta = 1% (relative).
- indistinguishable: Bonferroni CI inside [-1%, +1%];
- different: Bonferroni CI excludes 0 and its whole interval lies outside [-1%, +1%]
  on one side;
- otherwise: small/inconclusive (report the CI).
Detector validity: PC1 must be "different" at every size and PC2 must have Bonferroni CI
excluding 0 at every size; R-sym2 must be "indistinguishable" (else the instrument is not
trusted and that is the result).

### Power/resolution note
Expected C ~ 1 - exp(-1/2) ~ 0.39 for K ~ sqrt(p); binomial-scale relative sd of C on one
curve ~ 1/sqrt(0.39 N) ~ 0.6% (b=16), 0.15% (b=20); curve-to-curve means over 20-60
curves resolve well below 1%.

### Stop rule / budget
Single core, at most ~40 min wall-clock. Run order b=16, then 20, then 22. If the budget is
exceeded, stop and report the completed sizes; the cut is stated in RESULTS.md.
No change of metric, margin, arms or sizes after data; any deviation is a labelled
post-hoc analysis.
