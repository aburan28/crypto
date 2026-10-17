# PROTOCOL — 02-deuring-lift (frozen before data)

Campaign: pair-46 factor-base campaign, item 02. Frozen 2026-10-15 before any
data file in this directory was generated. Not edited after.

## Question

Does lifting the pair-46 crater (j = -96^3, End = Z[(1+sqrt(-19))/2]) to its
canonical (Deuring) lift over Q — the CM(-19) curve E: y^2 + y = x^3 - 38x + 90
(Cremona 361a1) and its quadratic twists E_d with chi_d(p) = +1 — give a
Xedni-style factor base with *lifted* (characteristic-0) relations, where a
generic Xedni lift (EV-XEDN-001 lineage) does not?

Prediction (stated in the task): fails by rank obstruction.

## Claim classes

- **Exact (derivation + exact computation):** Mordell-Weil rank and generator
  of E(Q) and of the twists E_d(Q), |d| small (PARI `ellrank`, with
  `ellanalyticrank` as a cross-check where 2-descent bounds do not close);
  that E reduces mod p to a curve isomorphic to the crater; the structure of the
  reduction image; the rank-obstruction statement (number of Q-relations among
  lifted factor-base points = |F| - rank of the height-pairing Gram matrix,
  bounded by the Mordell-Weil rank).
- **Empirical (toy population):** the distribution of the shortest mod-p
  relation among reductions of independent canonical-lift generators,
  normalized by the Gaussian heuristic; coverage of E(F_p) by small
  combinations; the generic-lift (Xedni) dependence rate.

## Population and sizes

- Toy primes: for each bit band b in {16, 20, 24}, the first 8 primes p >= 2^b
  (scanning upward, deterministic, no RNG) with chi_{-19}(p) = +1, good
  ordinary reduction of 361a1 and #E(F_p) prime (mimics pair 46: prime order,
  crater = reduction of the canonical lift). 24 primes total.
- Pair 46 itself (p = 5947698030651343) for the exact parts only (ranks,
  reduction identification, image of E(Q), height growth of target lifts). No
  DLP is attempted at pair-46 size.
- Twist set D: fundamental discriminants d with 1 < |d| <= 200, d != -19, one
  representative of each {d, -19d/(square)} isogeny pair (keep the smaller |d|),
  ranks computed once over Q (rank is a property of d only).
- Factor-base size for the lattice metrics: R = 6 (first six d in D, ordered by
  |d|, with rank E_d(Q) >= 1 and chi_d(p) = +1, using one generator each).
  If fewer than 6 qualify for a prime, that prime is reported and skipped for
  the R = 6 metrics.
- RNG seed: 20261015 (Python `random.Random(seed ^ p)` per prime) for random
  points / targets.

## Arms (same instrument on every arm)

Instrument I-Q (lifted relations): for a set F of points on curves over a
number field, rank of the canonical-height Gram matrix (PARI
`ellheightmatrix`, relative tolerance 1e-20 on the smallest eigenvalue relative
to the largest); lifted relations = |F| - rank.

Instrument I-p (mod-p relations): for a set F of R points of E(F_p) (prime
order N), discrete logs via PARI `elllog`, kernel lattice
L = {c in Z^R : sum c_i F_i = 0} (det N), LLL-reduced (`qflll`); metric
rho = lambda_1(L) / GH, GH = sqrt(R/(2 pi e)) N^(1/R).

Instrument I-cov: fraction of E(F_p) equal to sum c_i F_i, |c_i| <= 2 (exact
enumeration via logs).

- **A1 canonical, single curve:** F = {G} (generator of E(Q)); exact image of
  E(Q) in E(F_p); naive-height size of a lift of a uniformly random target.
- **A2 canonical, twist family:** F = {G_d reduced mod p : 6 twists}. I-Q (over
  the multiquadratic field; Gram is block-diagonal across distinct d), I-p,
  I-cov.
- **A3 generic Xedni lift:** r in {2,3,4,5} random points of E(F_p), lifted to
  integer coordinates in [0, p), curve coefficients lifted a_i + p k_i with
  rational k solving the incidence equations (free k set to 0). I-Q on the r
  lifted points. One trial per (prime, r).
- **NC (negative control) for I-p / I-cov:** F = 6 uniformly random points of
  E(F_p) (no lift exists / no structure).
- **PC (positive control):** I-Q and I-p on F = {G, 2G, 3G, G_d1, G_d2, G_d3}
  (two planted Q-relations; small mod-p relation (2,-1,0,...) must be found,
  rho << 1). Xedni positive control: r = 3 lifted points including P and -P
  (one planted relation).

## Metrics and decision

- Primary: lifted (Q-)relation count per arm, against the |F| relations an
  index-calculus factor base needs.
- Secondary: rho distribution A2 vs NC (difference of medians, bootstrap 95%
  CI, 2000 resamples; Mann-Whitney U two-sided); PC rho must be < 0.5 on
  every prime for the detector to be validated.
- Coverage A2 vs NC: difference of means with bootstrap CI.
- Decision: "rank obstruction confirmed" if A1/A2 give 0 lifted relations
  beyond planted ones, A2 rho indistinguishable from NC (CI contains 0 or
  |effect| < 0.1), and PC detects. "Advances" only if A2 shows lifted
  relations or rho significantly below NC.

## Stop rule / budget

One core, at most ~40 min wall-clock. Order of execution: exact parts, twist
ranks, then bands 16, 20, 24. If budget runs out, stop at the last completed
prime and report what was cut. Per-twist ellrank capped at 60 s; a twist whose
rank is not determined within cap is recorded "undetermined" and excluded.
