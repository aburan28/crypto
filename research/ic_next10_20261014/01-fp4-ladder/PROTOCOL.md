# Protocol: the E(F_{p⁴}) ladder with plain, D₄ and double-large-prime arms

Item 1 of `reports/Index calculus leads next steps.md`. Frozen 2026-10-14,
before any data was collected.

## Prior work, checked first

- `git log --all`, `research/`, `benchmarks/`: no end-to-end k = 4 run
  exists. The pieces that do exist are:
  - §11.17 of `RESEARCH_RESIDUAL_WALKS.md`: the S₄ solve, C₄ = 1.213·10¹².
  - §11.19: the linear algebra (LA) against rho, r∞ = 0.518 against the
    unfolded walk.
  - E18 (`RESEARCH_GLV_INVARIANT_FACTOR_BASES.md` §9.3): the τ_T column
    fold on 2-torsion curves. Every curve there had h = 4.
  - `gaudry_quartic_d4_20261013`: C₄(D₄) = 2.54·10¹⁰, 48× below S₄.
- No double-large-prime run at k = 4 is recorded.

## Hypotheses

- H1 (constant-factor lead, plain arm). With the relation phase priced at
  the solver's per-trial cost, S_IC / S_rho falls as about #E^{−1/4}, and
  the LA constant stays below 0.9 of rho-negation.
- H2 (D₄). The D₄ arm is below the plain arm at every size, by a constant
  factor. Its exponent equals the plain arm's.
- H3 (asymptotic lead, double large primes). The double-large-prime (DLP)
  arm's total cost exponent in #E is below 0.45. Theory predicts 3/8; rho is
  1/2.
- H4 (D₄ census). The construction that puts T in reach of the Y-line
  (f'(x₀) = c², c ∈ F_p^*) forces 4 | #E, so h ≥ 4. The report's "h ≤ 2" is
  then unattainable for the D₄ arm.

## Curves and sizes

- p ∈ {101, 181, 269, 397, 521, 769, 1033}, all ≡ 1 (mod 4), so
  log₂ #E ≈ 26.6 to 40.0. Three curves per p (seeds 1, 2, 3).
- Each curve comes from `quartic_folds::generate_translation_instance(p,
  seed, 16)`. This is the D₄ module's construction: a = c² − 3x₀²,
  b = −x₀³ − a·x₀. It gives #E = h·r with r prime and the smallest h it
  finds. The h reached is recorded per curve; that record is the census.
- Every arm runs on the same curve, the same subgroup ⟨G⟩ of order r, and
  the same planted target Q = [d]G.
- Stretch size: p = 1301 (log₂ #E ≈ 41.4), if the 2 GB memory cap allows.

## Arms

| arm | base and columns | relations | linear algebra |
|:--|:--|:--|:--|
| `plain` | x-line base {P : x(P) ∈ F_p}, one column per abscissa (about p/2) | every 4-term decomposition over the base | filter to a square core, then Wiedemann (§11.7) |
| `d4` | Y-line base {P : Y(P) ∈ F_p}, one column per {±P, ±(P+T)} (about p/4) | same | same |
| `dlp` | x-line base; small base F' = the first ⌈\|F\|^{3/4}⌉ points (Gaudry's \|F\|^{1−1/k}) | decompositions with at most 2 points outside F', cancelled by `gaudry_cubic::LargePrimeEliminator` (uncapped) | same, over \|F'\| + 1 unknowns |
| `dlp-null` (null control) | the DLP code path with F' = F | must reproduce `plain` exactly | same |
| `rho` | the unfolded r-adding walk, `quartic_folds::rho4_folded(None)`, 16 walks per curve | | |
| `rho-neg` (reference) | the negation walk, `rho4_folded(Negation)`, 16 walks per curve | | |

Shared rules for the IC arms:

- Each target contributes one decomposition, drawn uniformly from its
  verified decompositions that the arm admits (E18 rule 2).
- Every decomposition is re-added in the group before it is used.
- A logarithm counts only if [d]G = Q.
- The three x-line arms see the same trial stream: trial t's (a, b) is
  derived from (seed, t). The D₄ arm sees the same stream.

## Oracle and charge (the one modelling step, stated in advance)

Relations come from an exact enumeration oracle, not from the S₅ or D₄
solve:

1. Build a table of x(±P_i ± P_j ± P_k) once per base.
2. For each trial, probe it with R − s·P_l.
3. Recover each hit through the pair table.

The oracle returns every decomposition over the base, the set the
algebraic solvers return. The D₄ README checked this on 19 of 19 targets,
and §11.17 on 4 of 4.

Running the solver for every trial does not fit the budget. One S₄ solve
takes 25 min and 3.7 GB; one D₄ solve takes about 15 s. So each trial is
charged the solver's measured per-trial cost:

- C₄(S₄) = 1.213·10¹² F_p multiplications (§11.17) for the plain and DLP
  arms;
- C₄(D₄), measured in this run on the positive-control solves, for the D₄
  arm. The committed value, 2.54·10¹⁰, is the fallback.

C₄ is independent of p: the Macaulay matrix has fixed shape. Its only
p-dependence is root finding, about 10⁻⁶ of it.

### Metric

S is measured in ecbench's unit: group-addition equivalents (GAE) over √r.
One E(F_{p⁴}) affine addition is 97 F_p multiplications (measured, §11.17).
The S of an IC arm is the sum of these terms, divided by √r:

- the factor base;
- the trial scalar multiplications;
- trials × C₄ / 97;
- every relation re-add;
- LA multiplications mod r (Wiedemann, filtering and large-prime merges) at
  16/97 GAE each, §11.16's charge;
- the final [d]G check.

This is S_charged, the primary statistic. Also reported, never used for a
claim:

- S_executed: the oracle's own counted work in place of trials × C₄;
- every phase separately.

Rho's S is its group operations over √r, set-up included. The walk's
counted F_p multiplications are also reported.

### Fits

- Per arm and per phase, an OLS fit of log₂(cost) against log₂ #E over all
  curves, with a 95% bootstrap interval (2,000 resamples of curves).
- The slope of log₂(S_IC / S_rho-neg) against log₂ #E.
- The LA constant r_LA = LA GAE / (S_rho-neg · √r), pooled.

## Positive controls

- **D₄ solver against the oracle.** On one curve each at p = 269 and the
  largest p reached, run the real `gaudry_quartic_d4::solve_d4` (Macaulay
  degree 14) on 3 planted targets (sums of 4 base points) and 3 uniform
  ones. Its decomposition set must equal the oracle's on every target.
  These solves also give C₄(D₄) in situ.
- **Fit pipeline.** The relation-phase trial count of the plain arm must fit
  an exponent near 1/4 in #E: the trials needed scale with |F| ∝ p.

## Null control

- `dlp-null` must give the same trials, relations and LA counts as `plain`
  on every curve.

## Kill criteria (from the report, unchanged)

- **Plain.** Kill the constant-factor lead if either holds:
  - the measured slope of S_IC / S_rho against #E is shallower than −0.12
    (that is, greater than −0.12);
  - the in-unit LA constant against rho-negation exceeds 0.9 (no crossover
    at any size).
- **DLP.** Kill the asymptotic lead if the DLP arm's fitted total cost
  exponent in #E is ≥ 0.45.

## What would be inadmissible

- Changing C₄, the 97 or 16 charges, the |F'| rule, the sizes or the seeds
  after seeing data.
- Dropping a curve or an unverified run.
- Quoting an extrapolated crossover as measured.
