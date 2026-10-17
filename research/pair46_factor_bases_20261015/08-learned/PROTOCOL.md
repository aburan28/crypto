# 08-learned — PROTOCOL (frozen before any data; not edited afterwards)

Item 08 of the pair-46 factor-base campaign: a *learned* factor-base selector
(family `learned_sparse`, IDEA-20260807-*), required to beat the
greedy-optimised base (H-FB3 family, null) and a matched random base on
**held-out** targets.

## Claim class

Comparative, empirical (method A beats B on held-out targets), with an exact
side statement. The experiment can support "learned beats greedy and random
on held-out targets" or fail to.

Exact statement checked numerically (not assumed): for a symmetric set S of
points (|S| = 2L, L x-values, both lifts) in a group of prime order N, the
number of ordered pairs (P1,P2) in S x S with P1+P2 = R, summed over all R,
is exactly |S|^2. So the mean number of m=2 decompositions of a uniform
target is 4L^2/N for every base of size L, and the success probability is at
most |S+S|/N <= L(2L+1)/N. Against a random base (success ~ 1-exp(-mu),
mu = 2L^2/N) no base can win by more than about mu/(1-exp(-mu)).

## Instrument (same on every arm)

Toy prime-order curves E/F_p, y^2 = x^3 + a x + b, N = #E(F_p) prime.
A generator G is fixed; the table k -> x(kG) is built once per curve
(this is a test oracle for exact counting; the selectors see only
x-coordinates and target points). A base is a set V of L liftable x-values;
S = lifts of V. For a target R = rG, the m=2 decomposition count is
c(r) = #{(s1,s2) in S^2 : s1+s2 = r}, computed exactly for all r by cyclic
FFT convolution of the indicator of S (rounded to integers).

- Primary metric: held-out success rate H = fraction of held-out targets
  (all r in 1..N-1 not in the training set) with c(r) >= 1.
- Secondary: held-out mean count; training-set success rate (in-sample).
- Effect size: per (curve, seed) log-ratio log(H_arm / H_random); reported
  as mean with 95% percentile bootstrap interval over (curve, seed) units
  (2000 resamples, seed 0).

## Arms (matched size L on every arm)

- A0 `random`: L uniformly random liftable x-values (null control).
- A0' `random2`: an independent second random base (null-vs-null control;
  its log-ratio must have CI covering 0).
- A1 `greedy`: lazy greedy max-coverage of the training targets over the
  pool of all liftable x (first element random by seed).
- A2 `learned_greedy`: L1-regularised logistic regression on x-features,
  label = membership in the greedy base vs. random pool negatives; base =
  top-L pool x by score.
- A3 `learned_success`: same model, label = x participates in at least the
  median number of training-target decompositions against a reference random
  base of size 4L (decomposition-success label); base = top-L by score.
- PC1 `known_log` (positive control for the detector): x(kG), k = 1..L.
  Predicted: held-out success far below random (|S+S| <= 4L).
- PC2 (positive control for the learner): label = [x mod 7 == 3] on the same
  features; the learner must reach test AUC >= 0.95. Also report the test AUC
  of A2/A3 labels on held-out pool x (predicted ~0.5).

Features of x (attacker-legal): the low `bits` binary digits of x, x mod q
one-hot for q in {3,5,7,11,13}, Legendre(x), Legendre(x+1), parity of the
canonical y = sqrt(f(x)) (smaller root), 1 bias. L1 strength 1e-3, 300
proximal-gradient steps, step 0.5, balanced classes.

## Sizes, seeds

- Bits of p: 14, 16, 18. Curves per bit size: 4 (first prime-order curves
  from a seeded search, seed = bits). Base seeds per curve: 3.
- Operating points: mu = 2L^2/N in {0.25, 1.0} (L = round(sqrt(mu N / 2))).
- Training targets: T = floor(N/16) uniform r (seeded), disjoint from
  held-out by construction.
- Total units: 3 bits x 4 curves x 3 seeds = 36 per mu.

## Decision rule

"Learned advances" iff for A2 or A3, at every bits and both mu, the held-out
log-ratio vs random AND vs greedy has bootstrap 95% lower bound > 0, and the
effect does not shrink with bits. Otherwise null. The controls are valid iff
PC1 log-ratio upper bound < 0, random2 CI covers 0 and PC2 AUC >= 0.95; if a
control fails, the run is reported as inconclusive.

## Stop rule / budget

Fixed design, run to completion once. Budget 1 core, ~40 min. If over
budget, cut in this order and report the cut: 18-bit arm seeds 3 -> 1,
then the 18-bit arm entirely.
