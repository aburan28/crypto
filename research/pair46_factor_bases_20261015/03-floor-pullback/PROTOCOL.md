# PROTOCOL — 03-floor-pullback (frozen before data)

Campaign: pair-46 factor-base campaign, item 03. Frozen 2026-10-15 before any
data file in this directory was generated. Not edited after.

## Question

Pair 46 has a crater curve E (j = -96^3, End = O_K, K = Q(sqrt(-19))) and a
floor curve E' (End = Z[pi], conductor v = t+2 = 34489703) joined by a
descending isogeny phi: E -> E' of prime degree v. Does pushing a structured
crater factor base S (small-x, multiplicative-coset x) through phi give a
floor factor base phi(S) whose decomposition yield on E'(F_p) differs from a
random-x floor base of the same size? And what does membership in phi(S) cost
on the floor?

At pair-46 size phi has degree 3.4e7 and cannot be evaluated here at scale, so
the empirical part runs on toy analogues of the same arithmetic shape.

## Toy population (deterministic, no RNG)

Pair 46 satisfies pi + 1 = v*w with w = (1+sqrt(-19))/2 of norm 5, i.e.
p = 5v^2 - v + 1, N = #E(F_p) = 5v^2 - 2v + 4, twist order 5v^2.
Toy family: primes v in [29, 919] with p = 5v^2 - v + 1 and N = 5v^2 - 2v + 4
both prime (exactly pair 46's shape, so the v-torsion of the crater twist is
fully rational and the v+1 v-isogenies from the crater are F_p-rational).
From the qualifying list (sorted by v, length L), take the 6 instances at
indices round(k*(L-1)/5), k = 0..5 (duplicates dropped). p ranges up to ~4.2e6.

Per instance:
- crater E: y^2 = x^3 + 3j(1728-j)x + 2j(1728-j)^2, j = -884736, or its
  quadratic twist, whichever has order N (checked by [N]P = O on 8 points).
- kernel: from E's twist (order 5v^2), [5v] * random points give two
  independent order-v points; the v+1 subgroups are tried in order and the
  first whose Velu codomain has j != -884736 is taken as the descending phi.
  Floor is certified by: #E'(F_p) = N (on 8 points), j(E') != -96^3, and the
  twist of E' has a point of order divisible by v^2 (cyclic twist, as stated
  for pair 46's floor).
- phi is evaluated by Velu's formulas (x via sum over half-kernel; y = y*X'(x)).
  Homomorphism check: phi(P+Q) = phi(P)+phi(Q) on 20 pairs.

## Metric (same instrument on every arm)

For a factor base F (a set of m points of E'(F_p), closed under negation,
m even), the 2-term decomposition yield
  Y(F) = |{R in E'(F_p) : R = P1 + P2, P1, P2 in F}| / N
computed exactly: discrete logs of all points of F via a full log table of the
floor built by iterating a generator G' (independent of phi), then the set of
pairwise sums mod N. Identical code on every arm, every instance.

Base size: m = 2 * round(sqrt(N) / 2).

## Arms (all are point sets on the floor E'(F_p))

- A1 pull-smallx: phi(S), S = crater points with smallest x (x = 0,1,2,...,
  both signs), first m points.
- A2 pull-coset: phi(S), S = crater points with x in g*mu_d, d = smallest
  divisor of p-1 for which the coset holds >= m points, g = smallest
  non-residue... precisely g = smallest integer >= 2 that is a primitive root
  mod p; points taken in order of k (x = g * zeta^k, zeta = g^((p-1)/d)) both
  signs, first m.
- C1 floor-native small-x: floor points with smallest x, first m (the prior
  campaign's small-x base, placed on the floor).
- C0 random-x control (null): x drawn uniformly (seeded) until m points, both
  signs; seeds 0..9 (10 draws per instance).
- P+ positive control: {+-iG' : i = 1..m/2} (an interval in log space). Known
  exact yield ~ 2m/N, far from random. Instrument validity requires |z| > 3.
- Identity check X0: Y_crater(S) computed with the crater's own log table
  (generator G, phi(G) = G') must equal Y(phi(S)) exactly for A1, A2.

Secondary metric (same on A1, A2, C1, C0): fraction of base points whose
floor x-coordinate is < p/16 (x-uniformity of the image).

## Statistics / decision

Per instance: mu0, sd0 = mean, sd of Y over the 10 C0 seeds.
z(arm) = (Y(arm) - mu0)/sd0; effect ratio Y(arm)/mu0 with a 95% interval from
the C0 t-interval on mu0. Pooled: mean z over instances.
- Transport-null supported if |z| < 3 for A1 and A2 on every instance and P+
  has |z| > 3 on every instance.
- Any |z| >= 3 for A1/A2 is reported as an effect (sign stated).

## Exact part (derivation, checked numerically on the toys)

phi restricted to E(F_p) is a group isomorphism onto E'(F_p) (N prime, v != N).
Hence Y(phi(S)) = Y(S) exactly, and decomposition on the floor with base
phi(S) is the crater decomposition after phi^-1. Membership of Q in phi(S)
needs phi^-1(Q) = [v^-1 mod N] * phi_hat(Q): one degree-v isogeny evaluation
plus one scalar multiplication. Costs are derived (formula), not measured.

## Seeds, budget, stop rule

Seeds: random points for model/kernel selection use random.Random(1000+v);
C0 seed s uses random.Random(10*v + s). Budget: 1 core, <= 40 min wall clock.
Stop rule: run instances in increasing v; if cumulative wall clock exceeds
30 min, stop and report the completed instances and what was cut.
