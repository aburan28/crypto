# Off-the-shelf SMT solvers as the point-decomposition oracle

**Stage diagnostic.**  `S`, end-to-end cost and speedup are **unset**.
Protocol: [`PROTOCOL.md`](PROTOCOL.md), frozen 2026-10-08 before the run,
with two amendments (both budget-only, both made before any cited `n = 19`
or `n = 23` cell was written).  Instrument:
[`src/cryptanalysis/smt_oracle.rs`](../../src/cryptanalysis/smt_oracle.rs),
[`examples/pdp_smt_oracle_bench.rs`](../../examples/pdp_smt_oracle_bench.rs),
scored by [`examples/pdp_smt_oracle_score.rs`](../../examples/pdp_smt_oracle_score.rs)
from the frozen cell lines in [`results/`](results/).

## What was asked

Can a general SMT solver (cvc5 1.4.2, Z3 5.1.0), handed the same
Weil-descended Semaev `S_4` system the native CDCL+XOR solver and WDSat
solve, serve as the decomposition oracle — and if not, is the gap parity
reasoning?  Two SMT-LIB encodings were run: `bool` (`QF_UF`, `xor`/`and`)
and `bv1` (`QF_BV`, one-bit vectors, bit-blasted to the solver's SAT back
end).  Every arm received the identical rows (descended system, Kosters–Yeo
trace row, degree-2 Macaulay rows), enumerated models with blocking, checked
every model against the rows and lifted it through the group, and treated a
final `unsat` as a refutation.  Every verdict was checked against exhaustive
enumeration.

Curves: `icv1-f2m13-t181-515ee569`, `icv1-f2m19-t797-b6cf2467`,
`icv1-f2m23-t5197-69e76b73` (`K_0`, `y² + xy = x³ + 1`), `m = 3`,
standard subspace factor base with `ℓ = ⌈n/3⌉`.  Seed 20261010.  Budget
120 s wall per target per arm (amendment 2), `max_models = 64`.

## Results

Counts are compared only within an arm.  "Common refuted" cells are those
both the arm and the native arm refuted; the native-conflict and own-counter
sums and the wall medians are over exactly those cells.  Wall time is from a
2-vCPU container with no isolation and is a practicality note, not a speed.
cvc5 exposes no SAT conflict count; its `bool` decisions are the
`DecisionStep` resource counter and its `bv1` search runs inside CaDiCaL
with nothing exposed.

### n = 13

| arm | cells | found | refuted | exhausted | disagree | spurious | common refuted | native conflicts (common) | own conflicts (common) | own decisions (common) | median wall refuting, ms | native median wall (common), ms | wall ratio | median wall found, ms |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| native-cdcl-xor | 20 | 13 | 7 | 0 | 0 | 0 | 7 | 1151371 | 1151371 | n/a | 5589 | 5589 | 1.00 | 266 |
| cvc5-bool | 20 | 13 | 7 | 0 | 0 | 0 | 7 | 1151371 | n/a | 223431 | 2445 | 5589 | 0.44 | 967 |
| cvc5-bv1 | 20 | 13 | 7 | 0 | 0 | 0 | 7 | 1151371 | n/a | n/a | 7482 | 5589 | 1.34 | 500 |
| z3-bool | 20 | 13 | 7 | 0 | 0 | 0 | 7 | 1151371 | 878802 | 936412 | 12333 | 5589 | 2.21 | 3694 |
| z3-bv1 | 20 | 13 | 7 | 0 | 0 | 0 | 7 | 1151371 | 2098876 | 2188226 | 14679 | 5589 | 2.63 | 4184 |

`n = 19` (`ℓ = 7`): **partial, not cited.**  `results/n19-m3-ell7.jsonl`
holds the first 23 of 60 cells before the session was interrupted; it is
kept as a labelled partial run.  In those cells every arm, the native one
included, reached the 120 s cap on each refuting target, and the SMT arms
reached it on found targets too (native found one in 34.5 s).  `n = 23`
was not run.

## Scoring

Scored on `n = 13` only (`pdp_smt_oracle_score results/n13-m3-ell5.jsonl`):

- **P1 (correctness): PASS** — 100 cells, every decided verdict agrees with
  enumeration, every decomposition re-adds, zero spurious models.
- **P2 (Z3 parity cost ≥ 5× native): FAIL** — on the seven common refuted
  cells Z3 `bool` used 878,802 conflicts against the native arm's
  1,151,371 (0.76×); Z3 `bv1` used 2,098,876 (1.82×).
- **P3 (no SMT arm decides more than native): PASS** — every arm decided
  all 20.
- **P4 (gap does not close with `n`): not scored** — needs `n = 23`.

## Reading

At `ℓ = 5` the parity-reasoning explanation is not supported: a general
CDCL core with Tseitin-split `xor` rows refutes with fewer conflicts than
the native solver's Gauss–Jordan path, and cvc5's Boolean route refutes in
0.44× the native wall time on this host.  The degree-2 Macaulay rows every
arm receives may already carry most of the linear structure at this size,
leaving search heuristics to decide.  None of this is a speed: counts are
within-arm, wall time is unisolated, and `n ≤ 23` says nothing about
`n = 83` or `n = 131`.

## Class and what follows

**Not classed.**  The decision rule needs P4, which needs `n = 23`; the
protocol's bounded run was not completed.  The `n = 13` evidence rules out
the registered "boundary" verdict's mechanism (P2) without establishing an
engineering candidate (P4 unscored).  Completing the ladder is a follow-on:
`pdp_smt_oracle_bench --n 19 --targets 12` and `--n 23 --targets 8`, seed
20261010, `--native-wall-budget-s 120`, on an isolated host, then
`pdp_smt_oracle_score` over all three files.

## Completed ladder (2026-10-11, amendment 3)

`n = 19` (`ℓ = 7`, 12 targets, 139 factor-base points) and `n = 23`
(`ℓ = 8`, 8 targets, 275 points) were run with the frozen seed, budget and
rows; cvc5 is 1.1.2 (see amendment 3), Z3 5.1.0 as frozen.  Scored with
`pdp_smt_oracle_score` over the three files
([`research/pdp_followups_20261010/smt_ladder/score_n13_n19_n23.txt`](../pdp_followups_20261010/smt_ladder/score_n13_n19_n23.txt)).

### n = 19

| arm | cells | found | refuted | exhausted (120 s) | disagree | spurious | median wall found, ms |
|:--|--:|--:|--:|--:|--:|--:|--:|
| native-cdcl-xor | 12 | 4 | 0 | 8 | 0 | 0 | 17,541 |
| cvc5-bool (1.1.2) | 12 | 0 | 0 | 12 | 0 | 0 | — |
| cvc5-bv1 (1.1.2) | 12 | 5 | 0 | 7 | 0 | 0 | 60,452 |
| z3-bool | 12 | 1 | 0 | 11 | 0 | 0 | 25,536 |
| z3-bv1 | 12 | 0 | 0 | 12 | 0 | 0 | — |

The native arm's found cells took 4.8 s, 11.6 s, 17.5 s and 32.8 s
(0.14–0.59 M conflicts); its exhausted cells stop at 1.37–1.46 M conflicts.
cvc5-bv1 found targets 3, 6, 8, 9, 10 in 62.3, 11.6, 6.6, 61.3 and 60.5 s;
the native arm found targets 1, 6, 7, 10.  Two of cvc5-bv1's five are
targets the native arm did not decide, and one of the native arm's four
(target 1) cvc5-bv1 did not.  No arm refuted anything at `n = 19`.

### n = 23

Every arm exhausted its 120 s budget on all 8 targets (native: 1.09–1.23 M
conflicts per cell).  No cell was decided; nothing is comparable.

### Scoring

- **P1 (correctness): PASS** — 100 + 60 + 40 cells, every decided verdict
  agrees with enumeration, zero spurious models.
- **P2: FAIL** (unchanged from `n = 13`; no refuting cells exist above it).
- **P3: FAIL** — at `n = 19` cvc5-bv1 decided 5 cells to the native arm's 4.
- **P4: not scorable** — no arm refutes a cell at `n = 23` (none at `n = 19`
  either); the prediction's precondition is never met at the frozen budget.

### Class

By the registered decision rule, "P1 passes and P3 fails for some arm":
**cvc5-bv1 is an engineering candidate for the oracle constant only**, on
the strength of one extra decided cell out of twelve at `n = 19`, with a
different cvc5 version from the `n = 13` run.  The mechanism the protocol
hypothesised (parity reasoning) is refuted at `n = 13` and untestable above
it.  The exponent is untouched: at `n = 23` no solver decides a single
uniform target in 120 s, while the pair table decides each in microseconds
(`RESEARCH_IC_BOUNDARY_LEDGER.md` §5).  Carrying cvc5-bv1 to the §8a
`m = 83` gate would need a protocol of its own; nothing here proposes one,
because the §3 argument of the survey note says what it would find.

## Follow-ons (not started)

- (Done 2026-10-11) complete the ladder; see above.
- Wire `DecompositionStrategy::Smt` into the `ic` pipeline beside `wdsat`
  so the oracle can be charged inside relation collection.
- An XOR-aware SMT route (a Gauss–Jordan parity theory, or cvc5/Z3 over
  an XOR-native back end) is a different protocol; nothing here speaks to it.
- cvc5's `ff` theory is prime-field only; a `GF(2^n)` finite-field theory
  would let the un-descended `S_4` be posed directly and is also a
  different protocol.

## Disclosed, not cited

The 3-target smoke test at `n = 13`, seed 20261008, used to confirm the
instrument and pick the budget; the unbudgeted and the conflict-capped
`n = 19` attempts stopped under amendments 1 and 2 before any cell was
written.
