# Protocol: Hamming-weight normal-basis factor base under off-the-shelf SAT

Frozen 2026-10-11 before the confirmatory run.  Item 3 of
[`RESEARCH_PDP_SAT_SMT_NAGAO_SURVEY.md`](../../notes/index-calculus/RESEARCH_PDP_SAT_SMT_NAGAO_SURVEY.md) §6.
**Stage diagnostic.**  `S`, end-to-end cost and speedup are **unset**.

## Question

The one factor base on which a SAT solver has a structural edge over
Gröbner methods is the Hamming-weight base in a type-II optimal normal basis,
`F_w = {P : wt_ONB(x(P)) ≤ w}` on `K_0 : y² + xy = x³ + 1` over `F_{2^n}`:
Frobenius-stable at every `n`, membership by popcount, no storage, and a
cardinality constraint the solver takes natively while the Hamming ideal
costs Gröbner `|F_w|^{1.9–2.5}` calls
([`RESEARCH_HAMMING_IDEAL_PDP.md`](../../notes/index-calculus/RESEARCH_HAMMING_IDEAL_PDP.md)).
The repository's only measurement of it is the Python-era ONB Karatsuba
circuit under CryptoMiniSat (`ecc2k130/README.md`): `n = 9`: 0.11 s,
`n = 11`: 2.64 s, `n = 23`: nothing in 25 min (`k = 3` points, `w = 3`),
"the encoding has not been tuned".  This round replaces that path natively
(`AGENTS.md`, "Implementation language") and re-measures.

**Hypothesis.**  The compact descended encoding decides `n = 23` instances
the Python circuit could not, but the solver still enumerates: on refuted
targets its conflict count is at or above the `m − 1`-summand enumeration
floor, so the gain is a constant and the class is engineering.

## Instrument (Rust)

[`examples/hamming_nb_sat_panel.rs`](../../../examples/hamming_nb_sat_panel.rs).
Symbolic ONB arithmetic over Boolean affine forms (squaring = coordinate
permutation, product = `≤ 2n − 1` AND terms per coordinate using the fold
rule of `Type2Onb::mul`), `S₃` for `m = 2` and the chained
`S₃(x₁, x₂, u) = S₃(u, x₃, x_R) = 0` for `m = 3`, Sinz sequential counter
`wt(x_i) ≤ w`, lexicographic `x₁ ≤ x₂ (≤ x₃)`.  DIMACS with native `x` rows
for CryptoMiniSat; Tseitin-split 3-ary parity for CaDiCaL.  Every model is
lifted and its signed sums checked against `R`; a model that is not a point
relation is blocked and the solver re-run (`max_models = 64`), so `found` is
a verified relation and `unsat` refutes the algebraic system.  An exact
pair-table oracle decides every target first (`m = 2`: `|F_w|` subtractions;
`m = 3`: `|F_w|` subtractions into a `C(|F_w|, 2)` pair table).

## Frozen inputs

| item | value |
|:--|:--|
| curve | `K_0` in the type-II ONB (`Type2Onb`), `n ∈ {9, 11, 23}` |
| base | `w = 3` at `n = 9, 11`; `w = 2` and `w = 3` at `n = 23` |
| summands | `m = 2` and `m = 3` |
| targets | 6 uniform (random point of the curve, SplitMix64 seed `20261011`) + 2 planted (random `m`-sums of base points) per cell |
| solvers | CryptoMiniSat 5.11.15 (Ubuntu 24.04 package, Gauss–Jordan default, native XOR rows), CaDiCaL 1.7.3 (Ubuntu package, pure CNF) |
| budget | 300 s wall per (target, solver) at `n ≤ 11`; 180 s at `n = 23`; `max_models = 64` |
| host | shared 4-vCPU container, no isolation; other runs concurrent |

**Arm not available.**  The PB-XOR solver the survey named (LinPB, Yang–Meel
CP 2021) and CryptoMiniSat's BNN input are unreachable here: github.com
answers 403 under the environment's network policy and the packaged
CryptoMiniSat has no BNN parser.  Cardinality is therefore the sequential
counter in both arms; the PB-XOR arm stays pending.

## Unit

Verdicts and solver conflict counts (each solver's own counter, summed over
re-runs) against the enumeration floor `T_enum = |F_w|` (`m = 2`) or
`C(|F_w|, 2)` (`m = 3`) per target.  Wall time is a practicality note.  The
smoke test at `n = 9` with 3 + 1 targets of seed 1, used to confirm the
instrument, is disclosed and not cited.

## Predictions

- **P1 (correctness).**  Every decided cell agrees with the exact oracle;
  every reported relation re-adds to its target; zero disagreements.
- **P2 (reach).**  At `n = 23, w = 3, m = 3`, CryptoMiniSat decides at least
  one uniform target within 180 s (the Python circuit decided none in 25 min).
- **P3 (enumeration).**  On every refuted (`unsat`) uniform cell at
  `n = 23`, each solver's conflicts are at least `0.5 · T_enum`.
- **P4 (scaling).**  The median wall of found cells at `m = 3, w = 3` grows
  by at least `2^{1.5}` per field bit from `n = 11` to `n = 23`.

## Decision rule

- P1 fails: bug, nothing classed.
- P1–P4 pass: **engineering** (the native encoding extends reach; the solver
  remains an enumerator; the product law is untouched).
- P3 fails on two or more cells: the §1 falsifier of the survey is met for
  this base and a dedicated protocol follows.
- P2 fails: the base is as hostile to SAT as to Gröbner at this size;
  **boundary** for the solver route on this base.

Inadmissible: changing `w`, `m`, symmetry breaking or `max_models` between
arms; counting a timeout or a model-budget stop as a refutation; reading
`n = 23` as a statement about `n = 131`.
