# Protocol: combined-system solution-count census for algebraic factor bases

Item 2 of the IC next-10 campaign (`reports/Index calculus leads next steps.md`, "### 2.").
Frozen 2026-10-14 before any census data was collected.

## Prior-run check (done before this protocol)

- `grep` over `research/`, `benchmarks/`, `reports/` for "solution count", "census",
  "combined-system", "ideal degree": the only ideal-degree hits are the binary/cover fall-degree
  tables (`research/fall_degree_bounds_20261008/`, `RESEARCH_COVER_DECOMPOSITION_LEDGER.md`),
  which are extension-field systems. No prime-field `{S_{m+1}=0, L(x_i)=0}` count exists.
- `git log --grep census` finds only the GHS cover census (#1821), unrelated.
- Closest existing data: `research/ic_exponent_falsifiers_20261013/` enumerated decompositions on
  one Kummer base (p = 946177, t = 8) to test exponent sums. It did not compare counts with any
  control.
- Conclusion: the census has not been run.

## The question and an a-priori observation (stated before data)

Program §2.3 asks for the number of solutions over F̄_p of

    I_full = ( S_{m+1}(x_1, …, x_m, x_R),  L(x_1), …, L(x_m) )

for a base V = roots(L), |V| = B, and compares it with B^{m−1}.

Every family in scope (Kummer coset, Dickson split and torus, subgroup μ_B, small-x interval,
random subset) has L split over F_p with distinct roots. Then every F̄_p zero of I_full has
x_i ∈ V ⊂ F_p, so the count is the number of grid points of V^m on the hypersurface S_{m+1} = 0,
i.e. the number of ordered x-tuple decompositions of the target. Summing over all targets
counts signed m-tuples from the base, so for **any** base the mean over uniform targets is the
counting identity

    E_R[N(R)] ≈ (2·B_E)^m / n   (B_E = number of abscissae of V on E(F_p), n = #E(F_p)),

which is ≈ B^m / p when half of V lies on the curve (the "product law"). This is far below
B^{m−1} for every B < p, so the comparison with B^{m−1} cannot by itself separate structured
from random bases. B^{m−1}·2^{m−1} is instead the degree of the *partial* system that leaves one
membership equation out, which is what an enumerating solver walks through.

This observation is a prediction. The census tests it, and it also tests whether structure
moves anything that the mean does not fix: B_E/B, the per-target spread, and the F̄_p count
for a non-split L.

## Hypotheses

- **H0 (null, product law).** For each structured family, the mean count N̄ per target is within
  a factor of 2 of the random-subset control on the same curve, B and targets, and both follow
  (2B_E)^m/n.
- **H1 (live lead, as written in the report).** A family's N̄ is below B^{m−1}/16 at two
  consecutive t, *and* that is a property of the family, not shared by the random control.

## Arms

Factor bases, all of size B = 2^t, on the same curve and the same targets:

| arm | V | role |
|---|---|---|
| `kummer` | g·μ_B, g a seeded element of F_p^* outside μ_B | structured (PKM) |
| `subgroup` | μ_B (x^B = 1) | structured |
| `dickson` | {w + w^{-1} : w ∈ u·μ_B}, u ∈ F_p^*, u² ∉ μ_B (split case) | structured |
| `dickson-torus` | same with u in the norm-1 torus of F_{p²}, B ∣ p+1 (own primes) | structured |
| `smallx` | {1, …, B} | structured (§2.1) |
| `random` | uniform random B-subset of F_p (seeded) | null control |
| `onlycurve` | random B-subset of curve abscissae only (B_E = B) | positive control for the mean detector: predicts N̄ × 2^m |
| `dloginterval` | x({kG : 1 ≤ k ≤ B}) | positive control for the spread detector: same mean as `onlycurve`, counts concentrated on few targets |

Algebraic (msolve) validation arm, small sizes only: the ideal's F̄_p solution count from the
msolve rational parametrization (degree of the eliminating polynomial), for every arm above,
plus `random-nonsplit` (a random monic degree-B L over F_p, generally not split, roots in
F̄_p) and the partial system `(S_{m+1}, L(x_1), …, L(x_{m−1}))` whose predicted count is
B^{m−1}·2^{m−1}.

## Sizes

- p of 20, 24 and 28 bits; three prime-order curves y² = x³ + ax + b per size, found by a seeded
  search with PARI `ellcard`; primes p ≡ 1 (mod 2^11) for the split arms and p ≡ −1
  (mod 2^11) for `dickson-torus` (three curves each).
- m = 2 and m = 3; t = 2 … 10 (B = 4 … 1024).
- Targets per cell: 2000 for m = 2, 500 for m = 3, uniform random points R with x_R ∉ V,
  seeded, shared by all arms on a curve. (The report says 20; more targets are cheap here and
  only tighten intervals. Recorded as a deviation.)
- msolve validation: one 20-bit curve, m = 2 at t = 2 … 6, m = 3 at t = 2 … 4, 5 targets per cell,
  under `prlimit --as=700MB`, `-t 1`.

## Instrument

`examples/ic_next10_02_solution_count.rs` (Rust, single-threaded): for m = 2, for each signed base
point P, test x(R − P) ∈ V; for m = 3, a table of x(P_a + P_b) over ordered signed pairs, matched
against x(R − P_c). Distinct ordered x-tuples are counted. This enumerates E(F_p) decompositions.
Twist-type abscissae (x ∈ V not on E) can only contribute degenerate solutions requiring a pair
to cancel and x_R ∈ V, which the target rule excludes; the msolve arm checks that the
enumeration equals the F̄_p count on every small cell.

## Metrics

- N̄: mean over targets of N(R), with a 95% bootstrap interval (2000 resamples over targets).
- Ratios N̄/N̄_random, N̄/((2B_E)^m/n), N̄/(B^m/p), N̄/B^{m−1}.
- Spread: index of dispersion Var(N)/N̄ (Poisson ≈ 1) and the fraction of targets with N ≥ 1.

## Kill criterion (from the report, verbatim)

"If every structured family's count is within a factor of 2 of B^{m−1}, or of the
random-membership control, close all algebraic bases on prime fields and drop item 3. A family
whose count is below B^{m−1}/16 at two consecutive t is a live lead."

Reading fixed now: the kill fires if every structured family's N̄ is within a factor of 2 of the
random control in every cell with N̄_random ≥ 0.5 (smaller cells are too noisy at these target
counts and are reported but not scored). The live-lead clause is scored too. If the random
control is itself below B^{m−1}/16, the live-lead clause does not distinguish the family from
the null and is reported as non-discriminating rather than as a lead.

Positive controls must read as predicted (`onlycurve` N̄ ≈ 2^m × random; `dloginterval`
dispersion ≫ 1) or the run is INCONCLUSIVE.
