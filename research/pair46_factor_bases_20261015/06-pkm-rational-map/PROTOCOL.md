# 06-pkm-rational-map: PROTOCOL (frozen before data)

Pair-46 factor-base campaign, item 06. Written 2026-10-15 before any data was generated.
Not edited after the first data point; any change goes in a dated amendment at the bottom.

## Source being reproduced

C. Petit, M. Kosters, A. Messeng, "Algebraic approaches for the Elliptic Curve Discrete
Logarithm Problem over prime fields", PKC 2016, LNCS 9615, pp. 3-18,
doi:10.1007/978-3-662-49387-8_1 (author PDF: christophe.petit.web.ulb.be/files/16PKC_primeECDLP.pdf,
read in this session). Factor base F = {(x,y) in E(F_p) : L(x) = 0}, L a large-degree rational
map that is a composition L = L_n0 o ... o L_1 of small-degree maps, so x in V is replaced by the
low-degree chain x_{i,j+1} = L_j(x_{i,j}), 0 = L_n0(x_{i,n0}) (their System (4)), solved with
Groebner bases together with S_{m+1}(x_11,...,x_m1, X) = 0. Two instances: Sec. 3.3, V = subgroup
of F_p^* of smooth order (L_j = x^{p_j}, L_n0 = 1 - x^{p_n0}); Sec. 3.4, V = x-coordinates of a
coset of a smooth-order subgroup H of an auxiliary curve E'(F_p), L_j the x-parts of Velu
isogenies. Their experiments (Sec. 4, Tables 4-5): m = 2, |V| = 2^n1, p >= 2^(2 n1), Magma F4,
D_av^corr = 4,4,4,4,4,5,5,5,5,6,6 for n1 = 1..11 (Table 4), time growing about 8x per n1.

## Question

Does the PKM chained (low-degree) membership make the m = 2 relation search on prime-field
curves cheaper, in a way that changes the cost exponent relative to matched Pollard rho, compared
with the same set V given flat (Dickson images / subgroup x^d), small-x and random bases?

## Claim classes

- Exact (proved, checked by computation): X1 the chained and flat systems for the same V have
  the same solution set (eliminating auxiliary variables), so the relation yield of a PKM base
  equals that of the flat base with the same V, target by target; X2 every solution of the
  system lies in V x V subset F_p^2 (L splits over F_p), so nbsols is computed exactly by
  enumeration; X3 smoothness facts about p-1, p+1 for pair 46.
- Empirical: Groebner degrees and times, their growth exponent in p, relation-rate ratios,
  and the projected IC/rho cost ratio.
- Heuristic: extrapolation to p = 2^53 (pair 46) of fitted growth.

## Arms (m = 2, level n1 => |V| = d = 2^n1)

Prime for level n1: p = smallest prime >= 2^(2 n1) with 2^n1 | p - 1 (PKM Sec. 4.1 choice).
Target curve E: seeded random y^2 = x^3 + ax + b over that p with #E prime (naive count).
- A1 pkm_mult_chain: V = subgroup of order d of F_p^*; chain x_{j+1} = x_j^2 (n1-1 times),
  x_{n1}^2 - 1 = 0. 2*n1 variables.
- A2 subgroup_flat: same V, membership x^d - 1 = 0. 2 variables.
- A3 pkm_isogeny_chain: auxiliary E'/F_p (seeded search) with a point h of order d; V = x(R0 + <h>)
  with phi(R0) not 2-torsion; chain of n1 2-isogeny x-maps (Velu), last equation phi_n1(x) = chi.
- A4 dickson_chain (EV-PFDR image base, chained): its own prime p' = smallest prime >= 2^(2 n1)
  with 2^n1 | p'+1 and (p'+1)/2^n1 >= 3; V = {u + 1/u : u in u0*mu_d in the norm-1 torus},
  u0^2 not in mu_d, |V| = d; chain x_{j+1} = x_j^2 - 2, last x_{n1}^2 - 2 - lambda = 0
  (D_{2^n1}(x) = lambda). Its own prime-order target curve over p'.
- A5 dickson_flat: same V as A4, membership D_d(x) - lambda = 0 (flat).
- A6 small_x (Semaev): V = {0, ..., d-1}, membership prod (x - v) (flat).
- A7 random_x (negative control): seeded random d-subset of F_p, flat prod.
- A8 planted (positive control, relation detector only): V = x({jP : 1 <= j <= d}), targets
  R = aP + bP with 1 <= a, b <= d. No Groebner run.

## Metrics (same instrument on every arm)

- M1 relation yield per target (exact enumeration over V x V): valid decompositions
  R = +-P1 +- P2 with P1, P2 in E(F_p), x(Pi) in V (unordered); and nbsols = algebraic
  solutions {x1,x2} subset V with S_3(x1,x2,x_R) = 0. Ratio rho1 = mean valid / (|F|^2 * 2 / N)
  with |F| = #{x in V liftable} (one representative per +- pair); 95% Poisson interval on the
  pooled count. Prediction: rho1 ~ 1 for A1-A7; A8 >> 1.
- M2 Groebner cost of System (4) (arms A1-A7, same targets for arms sharing p): one
  pure-Python Buchberger (grevlex, sugar normal strategy, Gebauer-Moeller criteria), single core.
  Recorded: wall time T, D_max = max sugar degree of a processed pair, D_corr = max sugar
  degree of a pair whose reduction was nonzero (analogue of PKM's D_av / D_av^corr),
  #nonzero reductions. Correctness: GB checked against sympy's groebner on levels n1 <= 3,
  and against the enumerated solutions on every system (all must satisfy the GB; ideal must be
  zero-dimensional with solution count consistent).
- M3 matched rho: Pollard rho (r-adding walk, 16 partitions, Floyd) on the same E, 12 random
  DLPs per level, group operations and wall time in the same Python. Projected IC time per
  arm: (d/2 + 10) / (mean valid per target) * mean T (relation stage only; linear algebra
  ignored, which favours IC). Ratio IC/rho with a bootstrap (2000 resamples over targets) 95% CI.
- Growth exponents: OLS slope of log2(mean T) on log2 p over completed levels per arm, with
  bootstrap CI (resample targets within level); IC exponent = slope + 1/2 (relations needed
  ~ p / d ~ p^(1/2) attempts) vs rho exponent 1/2.

## Positive / negative controls

- Positive control (degree instrument): A1 must reproduce PKM Table 4 D_corr within +-1 at
  every completed n1 (paper: 4 for n1 <= 5, 5 for n1 = 6..8). Positive control (relation
  detector): A8 rho1 >> 1 (at least 10x).
- Negative control: A7 random-x, flat; expected rho1 ~ 1, D_max ~ d.

## Sizes, seeds, stop rule

- Levels n1 = 2, 3, 4, 5, 6, 7 ascending; 12 targets per level (seed = sha256("pkm06|arm-p|n1|t")).
- Per-system timeout 60 s (censored, recorded). An arm stops after the first level whose
  mean T > 15 s or that has a timeout. Total compute cap 35 min wall on 1 core; on hitting it,
  finish the current system and report the completed subset, naming what was cut.
- Pair 46 (p = 5947698030651343): only exact facts (X3) and heuristic extrapolation; no
  Groebner run at that size (out of budget by many orders of magnitude; stated, not hidden).

## Outcome -> licensed statement

- If the IC exponent (slope + 1/2) of the best chained arm has CI above 1/2 and IC/rho > 1 at
  every level: "at toy sizes with this solver, PKM chained bases do not beat matched rho and
  the gap widens with p" (null/negative for pair 46, empirical).
- If a chained arm's IC exponent CI lies below 1/2: report as an advance candidate needing an
  F4/msolve replication (not a claim).
- Otherwise inconclusive.

## Amendment 1 (2026-10-15, before the recorded run; after a discarded smoke test at n1 <= 3)

A code smoke test (n1 = 2, 3; output discarded, not in results.json) showed that the sugar degree
(D_corr, as preregistered) is not the quantity Magma reports: PKM's D_av is the largest F4 "step
degree", i.e. the total degree of the lcm of the critical pairs selected in a step. The instrument
therefore additionally records Dl_max / Dl_corr = max total degree of lcm(LM f, LM g) over processed /
nonzero-reducing pairs. The preregistered positive control (D_corr within +-1 of PKM Table 4) is
still evaluated and reported as stated; the lcm-degree comparison is reported separately and
labelled as the amendment metric. Nothing else changes.
