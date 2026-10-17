# PROTOCOL — 01-endo-chain (frozen before any data; not edited afterwards)

Date: 2026-10-15. Campaign: pair-46 factor-base campaign, item 01.

## Question
Does a factor base closed under the degree-5 CM endomorphism psi = w,
w = (1+sqrt(-19))/2 (chains x, R(x), ..., R^{k-1}(x), one unknown per chain because
log(psi^i P) = lambda^i log P), lower ECDLP cost below matched rho, or below the
abscissa (small-x) factor base, on CM(-19) curves of the pair-46 family shape?
Prediction (stated in the task): constant at best.

## Claim classes
- Exact: (E1) psi's x-map R is a degree-5 rational function with x(lambda P) = R(x(P)),
  recovered by interpolation and verified on independent points; (E2) one unknown per chain;
  (E3) cost-model identity: a single known-head chain of length k with an x-only table
  is a baby-step/giant-step search whose table entries cost c_E group-op equivalents each,
  so its optimum is sqrt(2 c_E N) group-op equivalents (versus sqrt(2N) for ordinary
  +-BSGS); (E4) psi-orbit canonicalisation (the glv-orbit trick) needs the orbit length
  ord(lambda) to be small; recorded exactly.
- Comparative/statistical: per-arm mean cost S = G / sqrt(N) with bootstrap 95% CI;
  decomposition rate (m = 2) of endo-chain sets vs random sets of equal size.
- Heuristic: cost formulas in the README used as predictions.

## Family and sizes
Curves: p = 5t^2 + 19t + 19 prime (so t^2 - 4p = -19 (t+2)^2, same as pair 46),
crater curve j = -96^3 (y^2 = x^3 + 3k x + 2k, k = j/(1728-j), or its quadratic twist)
whose order N = p+1-t or p+1+t is prime. For each target size b in {16, 18, 20, 22}
take the smallest t >= 2 whose p has exactly b bits and at least one of the two orders is
prime. lambda = root of X^2 - X + 5 mod N whose x-map interpolates to degree (5,4).
Pair 46 itself is used for the psi-vs-addition cost measurement and the E1 check only.

## Instrument (identical code path on every arm)
One standalone Rust binary (`endo_chain`, this directory). Affine Weierstrass arithmetic
over F_p with u128 reduction and extended-Euclid inversion. Counters: A = affine
additions/doublings, E = psi x-evaluations, L = Legendre/sqrt membership tests for chain
heads. Group-op equivalents G = A + c_E * E + c_L * L, where c_E, c_L are the ratios of
median ns per operation to median ns per affine addition, measured once per curve in
phase 0 (10 repeats x 200000 ops) and frozen before phase 1. Secondary metric: S with
c_E = 1 (most favourable to endo-chain), and raw counts.
Every solve is verified by checking d*P == Q.

## Arms (phase 1: full DLP solves, 40 targets per arm per size)
Targets: Q = d P with d uniform from a seeded RNG (seed 20261015 + size), the same targets
for every arm. F = ceil(sqrt(N)).
1. rho (reference): Pollard rho, r = 20 Teske adding walk, full storage, no negation.
   Theory S = sqrt(pi/2) = 1.2533.
2. rho-neg (reference): same on +-classes (canonical y <= (p-1)/2), useless self-collision
   -> restart at a fresh random combination (restart cost charged). Theory 0.886.
3. bsgs-pm: x-only +-table {jP : 1 <= j <= m}, m = ceil(sqrt(N)/2), giant stride 2m+1.
4. abscissa-k1 (small-x arm): F points with smallest x (x = 0,1,2,... on the curve), logs
   unknown; walk R_j = Q + (a0 + j) P; solve when two hits land on the same unknown.
5. endo-chain-k{4,16,64}: C = ceil(F/k) chains, heads at smallest on-curve x, each extended
   by k-1 psi x-evaluations, one unknown per chain; same walk; solve on the first pair of
   hits on the same chain (u eliminated, exact linear equation in d).
6. endo-chain-kF: one chain of length F with unknown head (same decoder).
7. endo-chain-known: one chain from the base point P (head log = 1), length
   k* = ceil(sqrt(N / (2 c_E))); solve on the first hit. This is the best possible
   endo-chain arm under the cost model (E3).
8. positive control for the fold: endo-chain-k16-nofold: the k = 16 table, but the decoder
   treats every element as its own unknown (as abscissa). Must cost more than k16 if the
   fold is real and visible to the instrument.
Walk start a0 uniform per (arm,target) from the seeded RNG.

## Phase 2: decomposition structure (m = 2), size 20 only
For FB types of equal size F2 = ceil(sqrt(N)/4): random points (negative control / null),
small-x, endo-chain k = 16 and k = F2, and an arithmetic progression {jP : 1..F2}
(positive control: its pairwise sums +-P_i +-P_j cover only ~4 F2 points, so its rate must
be far below the null). Rate = fraction of T = 3000 uniform random points R that equal
+-P_i +-P_j for some i, j (both in FB). Reported with Wilson 95% intervals; the endo-chain
claim "decomposes like a random set" is supported only if its interval overlaps the null's
and the positive control is separated from the null.

## Effect sizes
For each arm and size: mean S and percentile bootstrap 95% CI (2000 resamples, seed 7);
ratio of mean S to rho-neg mean S with bootstrap CI.

## Outcomes and what they license
- advances: some endo-chain arm has CI upper bound of S below rho-neg's CI lower bound
  at every size and the ratio does not grow with size. Predicted not to happen.
- constant_only: endo-chain beats the abscissa arm by a constant (the fold) but not rho.
- null/negative: endo-chain no better than BSGS/rho; with E3 as the exact reason.

## Stop rule and budget
At most 1 core, 30 minutes of total compute for phases 0-2. Order: phase 0, then sizes
16, 18, 20, 22 for phase 1, then phase 2. If the budget is reached, stop, and report the
completed subset and exactly what was cut. No re-runs with altered parameters after seeing
data except to fix a crash, which would be reported.
