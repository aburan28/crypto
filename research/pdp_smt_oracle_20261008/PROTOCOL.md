# Protocol: off-the-shelf SMT solvers as the point-decomposition oracle

Frozen 2026-10-08, before the confirmatory run.  **Stage diagnostic.**  `S`,
end-to-end cost and speedup are **unset**.  Nothing here touches the
`vs_rho` verdict of [`../../docs/ic/BOUNDARY_TARGETS.md`](../../docs/ic/BOUNDARY_TARGETS.md).

## Question

Index calculus on a binary Koblitz curve spends its oracle calls deciding
whether a target abscissa is a sum of `m` factor-base abscissae: the
Semaev system `S_{m+1}`, Weil-descended to `n` Boolean equations in
`m·ℓ` unknowns.  The repository solves that system with its own CDCL
solver with native XOR rows
([`semaev_sat`](../../src/cryptanalysis/semaev_sat.rs),
[`sat`](../../src/cryptanalysis/sat.rs)) and with Trimoska's WDSat
([`wdsat_oracle`](../../src/cryptanalysis/wdsat_oracle.rs)).  Neither is
an SMT solver.  The question is whether a general SMT solver, handed the
identical system, is a competitive oracle, and if not, where the gap is.

**Hypothesis.**  It is not, and the gap is parity reasoning: an SMT
solver's propositional core has no Gauss–Jordan step, so each dense
`xor` row of the descent is Tseitin-split and the search rediscovers
linear algebra by resolution, as the CNF control arm of
[`semaev_sat`](../../src/cryptanalysis/semaev_sat.rs) already shows for
the native solver.  Routing the rows through a bit-vector theory instead
(bit-blasting to the solver's SAT back end) does not restore it.

**Derivation of the floor.**  The oracle is not the method: the product
law of
[`../notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md`](../notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md)
§5 bounds any oracle that enumerates sub-tuples, and an SMT search over
the descended system is such an oracle.  A solver that merely searches
faster moves a constant; the exponent is fixed by the system.  So the
most this protocol can establish is an engineering class, never an
advance, and its negative outcome is a boundary on the solver route.

## Instrument (Rust)

- [`src/cryptanalysis/smt_oracle.rs`](../../src/cryptanalysis/smt_oracle.rs):
  emits the descended system as SMT-LIB 2 in two encodings, `bool`
  (`QF_UF`: `xor`/`and` over Booleans) and `bv1` (`QF_BV`: `bvxor`/`bvand`
  over `(_ BitVec 1)`), runs cvc5 or Z3 as a child process, enumerates
  models with blocking assertions up to `max_models`, checks every model
  against the original equations, lifts it through the group, and treats
  a final `unsat` as a refutation.  The rows are the native path's rows:
  the descended system, the Kosters–Yeo trace row, and the degree-2
  Macaulay rows (`macaulay_degree = Some(2)`), so the three arms differ
  only in the solver.
- [`examples/pdp_smt_oracle_bench.rs`](../../examples/pdp_smt_oracle_bench.rs):
  per cell, draws a target as a random multiple of the subgroup generator,
  decides it by exhaustive enumeration (`enumerate_decompose`), then runs
  the native arm (`sat_decompose_with`, default options) and every
  (solver, encoding) arm with the same `max_models = 64`, Macaulay degree
  2 and trace row.  Every verdict is checked against the enumeration and
  every returned decomposition is re-added in the group.

## Frozen inputs

| item | value |
|:--|:--|
| curves | `icv1-f2m13-t181-515ee569` (`K_0`, `n = 13`), `icv1-f2m19-t797-b6cf2467` (`K_0`, `n = 19`), `icv1-f2m23-t5197-69e76b73` (`K_0`, `n = 23`) |
| model | `y² + xy = x³ + 1`, the repository's irreducible for each `n` (`KoblitzCurve::new(0, n)`) |
| summands | `m = 3` (`S_4`) |
| factor base | standard polynomial-basis subspace `⟨1, z, …, z^{ℓ−1}⟩`, `ℓ = ⌈n/3⌉`: `ℓ = 5, 7, 8` |
| targets | `n = 13`: 20; `n = 19`: 12; `n = 23`: 8; random scalars from `StdRng::seed_from_u64(20261010)` |
| arms | `native-cdcl-xor`, `cvc5-bool`, `cvc5-bv1`, `z3-bool`, `z3-bv1` |
| solvers | cvc5 1.4.2 (static Linux release, CaDiCaL back end), Z3 5.1.0 (`z3-solver` 5.1.0 wheel binary) |
| budget | 120 s wall per solver call; `max_models = 64` |
| host | 2-vCPU Linux container, no isolation (`tools/isolated_bench.py` not run) |

The smoke test at `n = 13` with 3 targets of seed 20261008 (disclosed in
the PR, not cited) was used to confirm the instrument runs and to pick the
budget; the confirmatory seed differs.

## Unit

Verdicts and counts are the evidence; wall time is a practicality note.

- The native arm reports its own CDCL conflicts (cumulative over model
  enumeration), the repository's usual machine-independent effort count.
- Z3 reports `:conflicts` and `:decisions` (`-st`); zero counters are
  omitted and read as 0.
- cvc5 1.4 exposes no SAT conflict count.  On `bool` its
  `DecisionStep` resource counter is recorded as decisions; on `bv1` the
  search runs inside CaDiCaL and nothing is exposed, so its counters are
  `null`.
- **Counts are compared only within an arm across `n`**, never across
  solvers as if they were one unit.  Cross-arm comparison is on decided
  cells within budget and on wall time, read as practicality on this host.

## Predictions (pass/fail)

- **P1 (correctness).**  Every decided cell of every arm agrees with the
  exhaustive enumeration, every returned decomposition re-adds to its
  target, and no arm produces a spurious model.
- **P2 (parity cost, Z3).**  On refuting cells, Z3 in either encoding
  reports at least 5× the native arm's conflicts, at every `n`.
- **P3 (completeness).**  At every `n`, no SMT arm decides more cells
  within the 120 s budget than the native arm.
- **P4 (gap does not close).**  For every SMT arm that decides at least
  one refuting cell at both `n = 13` and `n = 23`, the ratio of its
  median refuting-cell wall time to the native arm's is at least as large
  at `n = 23` as at `n = 13`.

## Decision rule (registered)

- P1 fails: a bug, not a result.  The run halts for review and nothing is
  classed.
- P1–P4 pass: **boundary**.  Off-the-shelf SMT search, with or without
  the bit-vector route, does not replace native parity reasoning as the
  decomposition oracle; the solver route is closed at this level and any
  SMT-side continuation must add XOR reasoning (a Gauss–Jordan theory or
  an XOR-aware back end), which is a different protocol.
- P1 passes and P3 or P4 fails for some arm: that arm is an
  **engineering** candidate for the oracle constant only, to be carried
  to the §8a `m = 83` gate under its own protocol before any claim; the
  exponent is untouched by the floor above.
- P2 fails while P3 and P4 pass: the parity explanation is wrong but the
  verdict stands; record the counts and leave the mechanism open.

## Stop condition and inadmissible moves

Bounded: three `n`, 40 cells, five arms, one run.  It stops when P1–P4
are scored and the table is in `README.md`.

Inadmissible: changing `ℓ`, `m`, the Macaulay degree, the trace row or
`max_models` between arms; counting a timeout, `unknown` or an exhausted
model budget as a refutation; comparing cvc5 decisions or Z3 conflicts
with native conflicts as one unit; reporting a wall-time ratio from this
unisolated host as a speed; reading any outcome at `n ≤ 23` as a
statement about `n = 83` or `n = 131`.
