# Protocol: Torelli rarity heuristic on the n = 5 and n = 7 toy ladder (item 10)

Frozen 2026-10-14, before any curve data in this directory was collected.
Design source: `reports/Index calculus leads next steps.md`, "### 10.".
Background: `research/notes/ecc2k130/RESEARCH_ECC2K130_CURVE_CONSTRUCTION.md` §2
(toy ladder), `research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md` §6, §8.

## Check-first record

- `research/`, `benchmarks/`, git log: no Torelli/toy-ladder run beyond the
  2026-10-06 curve-construction ladder (genus 2-4 complete, genus-5
  hyperelliptic only). Non-hyperelliptic genus 5 and genus >= 6 were never
  enumerated. `conductor check` failed (no control plane on localhost:8080).

## Objects

`E_a : y^2 + xy = x^3 + a x^2 + 1`, trace t = -1 (a=0) or +1 (a=1).
`W_n(T) = T^{2n} - s_n T^n + 2^n` (Weil restriction from F_{2^n}),
`P_{A_n} = W_n / (T^2 - tT + 2)` (trace-zero part, dim n-1).

Pre-data facts computed with PARI/GP (structural, not curve data):

- `P_{A_3}`, `P_{A_5}` irreducible (dim 2, dim 4). Class numbers h_K = h_K+ = 1.
- `P_{A_7}` is **not** simple: it splits into two degree-6 factors
  `B_7, B'_7` (dim 3 each) for each sign, because Q(sqrt -7) ⊂ Q(zeta_7) and
  ord_7(2) = 3 (2 is not primitive mod 7). The n = 7 rung is therefore a
  structurally unfaithful analogue of n = 131 (where 2 is primitive and A is
  simple). It is run anyway, as the design asks, and reported separately.

## Hypothesis under test (the heuristic)

H (Torelli codimension heuristic): Jacobians are a "generic" subset of
principally polarised abelian varieties, so the mass of genus-G curves over
F_2 whose Jacobian has a fixed d-dimensional isogeny factor A is

  pred(G, d) = h_A * M_G * 2^{-c(G,d)},   c(G,d) = dG - d(d-1)/2,

where c is the codimension in A_G of {A x B : B in A_{G-d}}, M_G is the mass
(sum of 1/|Aut|) of the curve family considered, and h_A (number of ppav in
the isogeny class) is O(1) here (class numbers 1; we use a band
h_A in [1, 8]). The same formula is applied per family (uniformity assumption).

## Arms

Each arm is a complete (or uniformly sampled) family of smooth curves over
F_2; L-polynomials from point counts over F_{2^k}, k <= g.

| arm | genus | family | size |
|---|---|---|---|
| C3h | 3 | hyperelliptic y^2+hy=f, all models | 2^14 models |
| C3p | 3 | plane quartics, all | 2^15 - 1 |
| G5h | 5 | hyperelliptic, all (reduced mod y -> y + e) | equiv. 2^20 models |
| G5t | 5 | trigonal: plane quintics with a double point (delta = 1) at [0:0:1], all | 2^18 |
| G5c | 5 | canonical complete intersections of 3 quadrics in P^4, uniform random nets | 2^20 nets (cheap filter); 2^13 nets full L-poly |
| G6p | 6 | smooth plane quintics, all | 2^21 |
| G6h, G7h, G8h | 6, 7, 8 | hyperelliptic, all (reduced mod y -> y + e) | complete |

Hyperelliptic + trigonal + canonical is a partition of M_5, so genus 5 at
n = 5 is complete up to the sampling of G5c. Genus 6-8 is partial (plane
quintics + hyperelliptic only); its coverage is reported as mass.

Targets (both signs where defined):
- n = 5: `A_5` (dim 4) at genus 5 and 6; `W_5` (dim 5) at genus 5 and 6.
- n = 7: `A_7` (dim 6) at genus 6, 7, 8; its factors `B_7`, `B'_7` (dim 3) at
  every genus enumerated (3, 5, 6, 7, 8).
- Positive control (known to occur): `A_3` (dim 2) for both signs, in every
  arm. Instrument control: C3h and C3p must reproduce the 2026-10-06 counts
  (C3h: 6,144 smooth, 96/96 A_3 hits; C3p: 56/224 A_3, 56 W_3 = Klein).

## Metric

Observed carrier mass per arm (models carrying / |group acting on models|) and
the ratio observed / pred(G, d). For G5c the mass is estimated from the
sample (binomial / Poisson 95% intervals). Positive control: observed A_3 mass
vs pred(G, 2) in the same arm (calibrates the heuristic's constant). Smoothness
is checked by absence of singular points over F_{2^k}, k <= 6 (plane curves,
hyperelliptic exact gcd criterion); every target hit is re-checked to k <= 8.

## Kill criterion (from the report)

- Zero carriers of A_5 through genus 6 and of A_7 through genus 8:
  "confirms rarity" (verdict NEGATIVE for reopening the route).
- More than one carrier at genus <= dim A + 2 (genus <= 6 for n = 5,
  <= 8 for n = 7): reopens the route (verdict POSITIVE).
- Reported alongside: whether the heuristic's predicted count was large
  enough for zero to be informative (power). If pred << 1 at every rung, a
  zero is consistent with the heuristic but does not test it; that limits the
  verdict and will be said.
