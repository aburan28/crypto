# Protocol: item 02, the G2 equivariant pass as a Koblitz decomposition oracle

Frozen 2026-10-14, before any main-arm data. Rules: [`../RULES.md`](../RULES.md).
Prior work: [`../../ic_next10_20261014/06-equivariant-f4-m4/README.md`](../../ic_next10_20261014/06-equivariant-f4-m4/README.md) §3,
[`RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md`](../../../RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md).

## 1. The incumbent, identified in the code

The Koblitz IC pipeline decomposes a target in
`koblitz_index_calculus::groebner_decompose_with_geometry`
(`src/cryptanalysis/koblitz_index_calculus.rs`). It:

1. builds the **target-fixed** Semaev system with
   `polynomial_reuse::build_decomposition_system_reusing` (the target's `x`
   is a constant of the system);
2. solves it with `koblitz_groebner::solve_boolean_system_filtered`, options
   `SolveOptions { engine, split_rule: split_rule_default(), .. }`, where
   `engine` is `SolverEngine::default()` = `InheritedF4 { max_degree: 3 }`;
3. inside the solver every node's algebra is `reduce_system` /
   the inherited-basis path: one Macaulay matrix at degree 3 (systems with
   ≥ 24 variables go straight to the top degree), reduced by the crate's
   `F₂` kernel, consuming only a contradiction (`1`) or forced/linear rows,
   then splitting.

**No system the pipeline hands to `reduce_system` is Frobenius-symmetric.**
The target-fixed system is not invariant under `σ` (squaring maps the
target `x_R` to `x_R²`). The Frobenius-symmetric system exists only in the
*barrel* formulation of the G2 scope note, where the target carries its own
one-hot selector `s₀` (`X_R = Σ_k s₀[k] R^{2^k}`). Fixing `s₀ = e₀` (the
gauge) gives back a target-fixed system: that **gauge-fixed barrel system
is what the incumbent oracle would solve**, so it is the incumbent here.

- **Incumbent oracle (I).** `solve_boolean_system` with the pipeline's
  default options (`SolverEngine::default()`, `split_rule_default()`, no
  node budget cap) on the gauge-fixed barrel system, written as
  `F2BoolPoly` exactly as the pipeline's engine takes it.
- **Incumbent step (I_step).** The Macaulay step that engine runs at the
  root: the degree-3 Macaulay matrix of the gauge-fixed system, columns in
  descending degree, reduced by the crate kernel `gf2_elim::rref_counted`.

Structural fact recorded before data: the gauge-fixed system is
**quadratic** (`x_R` is a constant, so `x₁x₂x_R` is quadratic) while the
symmetric system is **cubic** (`X_R` is linear in `s₀`). At the same
Macaulay degree 3, the specialisation `s₀ = e₀` of the symmetric span lies
inside the gauge-fixed span and can be strictly smaller (it lacks
`S₃|_{gauge} × variable`). The step-level gate below tests this.

## 2. The candidate

The equivariant pass of `examples/ic_next10_06_g2_full_pass.rs`
(bit-identical to full `F₂` RREF at `n = 7, 11, 13, 17`), substituted at
the only node where the system is `σ`-symmetric: the root, before the gauge
is fixed.

- **Candidate step (C_step).** Symmetric barrel system at degree 3: row and
  column orbits, trivial block over `F₂` (reduced with the same crate
  kernel `gf2_elim::rref_counted`), character blocks over `F_{2^d}`,
  inverse transform of the **degree-fall rows**, specialisation `s₀ = e₀`,
  then RREF of the specialised rows in the gauge ring to read off what a
  solver consumes (contradiction, linear rows). Assembly is charged in full:
  orbit enumeration, block assembly, inverse transform, specialisation.
- **Candidate oracle (C).** C_step, then the incumbent solver on the
  gauge-fixed generators plus the reduced specialised fall rows. Charged:
  C_step + that solve. (An iterated equivariant F4 is out of reach: a split
  on any variable breaks `σ`.)

No `src/` file is changed: no node of the production pipeline receives a
`σ`-symmetric system, so a `reduce_system` hook would be dead code. Both
arms are driven from `examples/ic_improve_02_equivariant_oracle.rs`, which
calls the public pipeline solver for the incumbent.

## 3. Paired design

- Curve: `y² + xy = x³ + 1` (Koblitz, `b = 1`) over `F_{2^n}`; field
  polynomial `find_irreducible(n)`; window `W = span(1, z, …, z^{ℓ−1})`,
  `ℓ = 3`; `m = 2` summands; normal basis from a seeded RNG.
- Targets: `x_R` uniform in `F_{2^n}^*` from `StdRng::seed_from_u64(20261014 ⊕ n·1000 + t)`.
  Each target is run by every arm in the same process, arms interleaved
  per target, single thread, `nice`, `prlimit --as=4e9`.
- **Step sizes:** `n ∈ {7, 11, 13, 17, 19, 23}` (these have every
  nontrivial character orbit of size `d = ord_n(2)`, which the trace
  descent needs). Targets: 8 at `n ≤ 13`, 6 at 17, 4 at 19, 3 at 23
  (fewer if the 60-min job cap binds; reported).
- **Oracle sizes:** `n ∈ {5, 7, 11, 13}` (the pipeline's 64-bit monomials
  hold the gauge-fixed system's `4n + 2ℓ` unknowns up to `n = 14`).
  Targets: 12 per size.

## 4. Metric

- Primary: per-target **thread CPU seconds** (`CLOCK_THREAD_CPUTIME_ID`),
  each arm's whole cost including assembly. Field-operation counts are not
  commensurable across arms (`F₂` word XORs vs `F_{2^d}` table products), so
  they are recorded as secondary columns: incumbent kernel word XORs;
  candidate trivial-block word XORs plus character-block `F_{2^d}` products;
  solver `f4_word_ops_thread`.
- Paired ratio per target `ρ = cost(I) / cost(C)` (> 1 means the candidate
  is cheaper). Per size: geometric mean of `ρ`, 95 % percentile bootstrap
  CI (10,000 resamples of targets, fixed seed).
- Scaling: OLS of `log ρ` on `log n` over all targets; slope `b` with a
  95 % bootstrap CI (targets resampled within each size).

## 5. Gates (correctness)

- **Step gate:** on every target, I_step and C_step agree on the
  refutation bit and on the space of linear consequences (identical RREF of
  the degree-≤1 rows in the gauge ring). Containment
  `span(spec(C rows)) ⊆ span(I rows)` is checked (rank test) and reported.
- **Oracle gate:** on every target, I and C return the identical set of
  roots, equal to a brute-force enumeration of `(k₁, c₁, k₂, c₂)` with
  `S₃(x₁, x₂, x_R) = 0`.
- An arm that fails its gate on any target fails at that level; no speed
  verdict is given for it there.

## 6. Decision rule

At each level that passes its gate:

- **ADVANCES** — ratio CI lower bound > 1 at every size of that level.
- **REGRESSES** — ratio CI upper bound < 1 at every size.
- **MATCHES** — CI covers 1 at every size.
- **TRADES** — mixed across sizes, or the step advances while the oracle
  regresses (or the reverse).
- Exponent-level statement only from the slope CI over ≥ 4 sizes.

The overall verdict is the **oracle** level's (the item's question); the
step level is reported as the primitive.

## 7. Controls

- **A/A (null):** the incumbent step and the incumbent oracle are each run
  twice per target; the ratio CI must cover 1.
- **Positive control:** at `n ∈ {7, 11, 13}` the step instrument must see
  the known win of the equivariant pass over the crate kernel's full `F₂`
  RREF of the *symmetric* matrix (prior: 44–73× against a naive kernel),
  i.e. ratio CI lower bound > 1.

## 8. Budget

4 GB address-space cap per job, 60 min per job, single-threaded.

## Amendment A1 (2026-10-14, after a one-target pilot at n = 11, before main data)

The pilot (n = 11, target 0; kept as `runs/pilot.jsonl`) showed the
candidate oracle cut the incumbent's splits from 4,265 to 29. Appending
*any* root degree-fall rows to the generators could do that, equivariant or
not. One arm is added to the oracle level, as an attribution control:

- **F (fall-row augmentation, no equivariance).** The incumbent's own
  degree-3 Macaulay step on the gauge-fixed system (crate kernel), its
  degree-fall rows (pivot degree < 3) appended to the generators, the same
  solver. Charged: the step plus the solve. Gate as for C.

Reported: `I/C` (the item's verdict), `I/F`, and `F/C` (what the
equivariant pass adds over plain fall-row augmentation), each with its CI
under the same decision rule. Nothing else changes.
