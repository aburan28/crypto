# RESULTS — 03-floor-pullback (pair-46 factor-base campaign, item 03)

Protocol: `PROTOCOL.md` (frozen before data), `PROTOCOL_EXT.md` (population
extension E1, frozen before E1 data; also records deviation D1).
Raw records: `results.json` (main run, `run.log`), `results_ext.json`
(extension, `run_ext.log`; first attempt crashed in the coset arm before
writing v=1549, log kept as `run_ext_attempt1_crash.log`; the fix removed a
p < 1e7 guard on the divisor enumeration, and the rerun reproduced the three
instances that had finished, value for value). Summary: `summary.json`. Graph:
`fig_yield.svg` / `fig_yield.png`. Report: `report.pdf`.

Host: shared 4-CPU container, one core (`taskset -c 3`), Python 3 + numpy.
Wall clock: main run 13 s, extension about 3 min. The budget was 40 min.

## Verdict

**Null (an exact obstruction, confirmed empirically).** A factor base pushed
from the crater through the descending v-isogeny gives the same decomposition
yield on the floor as the crater base, exactly. That yield is the random-set
yield: |z| <= 1.40 on all 8 toy instances for both pushed bases. Membership in
the pushed base costs an extra dual-isogeny evaluation and a scalar
multiplication per test, so on pair 46 this direction is strictly worse than
working on the crater.

## Exact statements (derivation; verified numerically on every toy)

1. N = #E(F_p) is prime and deg phi = v != N, so phi: E(F_p) -> E'(F_p) is a
   group isomorphism. For any S, phi(S)+phi(S) = phi(S+S), so
   Y_floor(phi(S)) = Y_crater(S). Every decomposition R' = phi(P1)+phi(P2) on
   the floor is the crater decomposition phi^-1(R') = P1+P2. The floor adds no
   relation that the crater lacks. Check X0: the crater yield and the floor
   yield, computed with two independent log tables, were bit-identical for
   A1 and A2 on all 8 instances.
2. Membership: Q in phi(S) iff [v^-1 mod N] * phi_hat(Q) is in S. One
   evaluation of the degree-v dual isogeny plus about log2 N doublings. At
   pair 46 that is about 6e3 field ops (sqrt-elu, figure from the task) plus
   about 53 doublings and 27 additions, roughly 1e3 field ops (derived). The
   crater small-x test is one integer comparison. In summation-polynomial
   form, x(phi(P)) is a rational function of degree v in x(P), so the
   pushed small-x set's defining relation has degree v = 3.4e7 per variable
   (derived). That pushes the membership degree up, the wrong way.
3. Pair 46 structure (exact for this family shape): pi + 1 = v*w and
   End = O_K, so the crater twist group is O_K/(v*w), which is Z/5v x Z/v.
   Numerical checks on the toys: on all 8, [5]P produced two independent
   order-v points, so the v-torsion is fully rational. At v=359, [5v]P = O
   for every sampled point. This confirms the task's "derived, unverified"
   structure for this family shape. The pair-46 curve itself was not
   checked. The floor twist is cyclic
   (some [5v]P != O), and all v+1 v-isogenies from the crater are F_p-rational
   (pi acts as -1 on E[v]).

## Empirical results (toy population)

Toys: p = 5v^2 - v + 1, N = 5v^2 - 2v + 4, both prime, as for pair 46. The
preregistered range v in [29, 919] holds only 3 such v. E1 adds all 5 in
[920, 2000]. m = 2*round(sqrt(N)/2). Metric: Y = |F+F|/N on the floor, exact.

| v | p | m | mu0 (C0, 10 seeds) | sd0 | z A1 pull-smallx | z A2 pull-coset | z C1 floor-smallx | P+ ratio |
|---|---|---|---|---|---|---|---|---|
| 359 | 644047 | 802 | 0.39365 | 0.00067 | -0.09 | -0.21 | -0.61 | 0.0063 |
| 409 | 835997 | 914 | 0.39371 | 0.00070 | +0.82 | +1.28 | +0.25 | 0.0056 |
| 443 | 980803 | 990 | 0.39365 | 0.00062 | -1.40 | +0.39 | +0.87 | 0.0051 |
| 1213 | 7355633 | 2712 | 0.39364 | 0.00018 | +0.80 | +0.43 | +0.33 | 0.0019 |
| 1283 | 8229163 | 2868 | 0.39357 | 0.00021 | -1.27 | -0.40 | +0.23 | 0.0018 |
| 1409 | 9924997 | 3150 | 0.39351 | 0.00026 | +0.66 | +0.73 | -0.63 | 0.0016 |
| 1549 | 11995457 | 3464 | 0.39358 | 0.00018 | +1.23 | +0.75 | +1.41 | 0.0015 |
| 1613 | 13007233 | 3606 | 0.39355 | 0.00014 | +0.47 | -0.66 | +0.25 | 0.0014 |

- Effect sizes: the yield ratio Y(arm)/mu0 lies in [0.9978, 1.0023] for A1
  and A2 on every instance (|effect| <= 0.23%). The 95% intervals in
  `summary.json` cover the uncertainty in mu0 only, not the single-draw spread
  of the arm, so the z-score is the test. Mean z over 8 instances: A1 +0.15,
  A2 +0.29, C1 +0.26.
- Positive control P+ (log-space interval) has |z| >= 562 on every instance
  (ratio <= 0.0063). The instrument detects structured bases.
- The random-set heuristic 1 - exp(-m^2/(2N)) is 0.3932-0.3936, in line with
  mu0.
- Secondary metric: the fraction of pushed-base points with floor x < p/16 is
  0.055-0.076 (A1) and 0.057-0.072 (A2). The C0 mean is 0.058-0.065 and
  uniform is 0.0625. The pushed small-x set does not stay small-x on the
  floor; C1 is 1.0 by construction.
- Decision rule (preregistered): transport-null is supported, with |z| < 3 for
  A1 and A2 everywhere and |z| > 3 for P+ everywhere.

## What was cut or changed

- Pair-46-size isogeny evaluation (v = 3.4e7) was not run. The yield identity
  is exact, so it would add nothing. The cost there is derived, not measured.
- The preregistered main population had only 3 instances, not 6, because few
  v in range make both p and N prime. Extension E1 added 5 under a separate
  frozen protocol.
- D1: the order-v kernel points come from [5]*P, not [5v]*P as PROTOCOL.md
  said. The twist exponent is 5v. Arms, metric and decision are unchanged.
- Only m = 2 decompositions and one base size m ~ sqrt(N) were run. No
  relation-collection or linear-algebra phase was run.

## Relation to prior work

This extends the isogeny-image family (EV-ICI-0bdaaa) and the small-x and
coset nulls. Moving to an isogenous curve through a rational isogeny cannot
change any group-theoretic yield statistic of a transported base. A floor
base is only new if it is defined natively on the floor. C1 is the native
floor small-x base, and it is also random-like (|z| <= 1.41).

## Reproduce

```
cd research/pair46_factor_bases_20261015/03-floor-pullback
taskset -c 0 python3 run_floor_pullback.py results.json
taskset -c 0 python3 run_floor_pullback.py results_ext.json --ext
python3 analyze.py && python3 make_report.py
```
