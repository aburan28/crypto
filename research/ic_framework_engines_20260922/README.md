# Matched engine suite: Gröbner, crossbred and exhaustive engines on Weil-descent systems

This directory freezes the matched baseline/candidate suite that
`AGENTS.md` §8 asks every index-calculus performance iteration to run,
for the polynomial-system engines behind the benchmarking framework's
`SystemSolver` plug point (`docs/ic/FRAMEWORK.md`).  The protocol, the
boundaries it is read against and the targets it was run to test were
declared in `research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md`
§17.1–17.2 and committed before the frozen run; §17.3 onward reads it.

## Why not the WDSat regression suite

The frozen suite at `research/index_calculus_baseline_20260914/regression/`
pins the WDSat SAT protocol: its inputs, its counters and its
algebraic-only rejection control are SAT-solver quantities, and its own
README says it is not an adapter for other engines and that broader
solver experiments use the parent accounting contract.  The engines here
— F4, three hybrid F4/F5 engines, crossbred, fast exhaustive search — do
not have a conflict counter and are not SAT solvers, so this is the
equivalent matched suite §8 provides for, frozen under
`research/index_calculus_baseline_20260914/ec_index_calculus_contract.json`.

What it keeps from that contract:

| contract item | here |
|:--|:--|
| identical targets for every comparison | every engine solves every target of a cell in the same process; each system is fingerprinted with blake3 and `manifest.json` freezes the fingerprints |
| timeout is never UNSAT | a budget verdict is `budget`, never `unsatisfiable`, is counted in every row, and is not repeated on that target |
| every counter names its unit; missing is null | each engine's `op_unit` is carried; an engine without a degree reports none, never zero |
| first-solution and enumeration leaderboards apart | `sat-cdcl` is marked `first` and checked only for membership in the reference's solution set |
| a boolean root is not yet a group relation | the stage parts stop at boolean roots (status `SAT_NOT_YET_POINT_VERIFIED`); the bench parts lift, verify on the curve, and recover the logarithm |
| fixed setup charged once; cold cost first | the bench parts are cold whole runs, setup included |
| full-rank recovery and scalar verification for full-DLP claims | every bench row carries `verified`; `compare.py` fails on an unverified row |
| ratios to the floor and the matched reference | bench rows carry `S`, the counted-rho reference on the same instance, and the pair-table row on the same base |

Verdicts map to the contract's statuses as: `solved` →
`SAT_NOT_YET_POINT_VERIFIED` at the stage (verified on the equations,
not yet on the curve); `unsatisfiable` → `SOLVER_REPORTED_UNSAT`,
checked against the exhaustive reference on every target; `budget` →
`TIMEOUT`; `declined` → `NOT_RUN`.

## Layout

| path | what |
|:--|:--|
| `suite.json` | the frozen protocol: parts, cells, seeds, engines, repetitions, budget, the analysis and the declared targets |
| `sweeps/` | the whole-pipeline sweep files the bench parts run |
| `run.py` | runs parts into `results/<run-id>/`, append-only, with a provenance record per part |
| `compare.py` | reads a run, checks it, evaluates the declared targets, writes `compare.json` and `compare.md` |
| `manifest.json` | the frozen run's input fingerprints and verdict digests, per cell |
| `confirm_20260923.json`, `strength_20260923.json`, `hybrids_20260923.json`, `hybrids_34_20260923.json` | the checks declared after parts of the frozen run were read (ledger §17.3–§17.5), each its own frozen protocol |
| `results/baseline_v1/` | the frozen run |
| `results/confirm_v1/`, `strength_v1/`, `hybrids_v1/` (a failed first attempt, kept), `hybrids_v2/`, `hybrids_v3/` | the checks |

## The protocol, in brief

**Stage parts** (`ic descent --solver …`): the Weil descent of the
Semaev polynomial over the seeded targets of the ledger's §14/§16
table, eight targets a cell, families `K` (Koblitz) and `R` (random
curve), three repetitions interleaved per target with the engine order
rotated each repetition, a 120-second budget per call, one thread.
The reference is `fes-f2` on the quadratic (two-summand) cells and
`exhaustive` on the three-summand ones — and `fes-f2-wide`, the vector
form of the fast search, wherever a run lists it (the checks from ledger
§17.4 on); the baseline is `buchberger-f2`.  Every engine's answer on every target is compared
with the reference's.

**Bench parts** (`ic bench --sweep …`): whole runs — factor base,
walk, descent oracle with each engine, relation matrix, logarithm,
verification — on `E(F_{2^13})`, `E(F_{2^15})` (three curves each) and
`E(F_{2^17})`, two planted targets each, with the pair table as the
oracle reference row and sixteen counted rho runs per instance.

**Wall time is the common unit and a practicality note.**  Each engine
counts in its own unit (word XORs of an elimination, monomial
operations, conflicts, monomial tests), and none of the candidates'
counts covers its whole run, so the framework prices them at measured
wall time (ledger §15.4) and this suite compares wall times paired on
the same target in the same process.  The ratios are host-dependent;
the native counts stay in every report as each engine's own regression
measure.

## Running it

```
cargo build --release --bin ic
python3 research/ic_framework_engines_20260922/run.py --run-id baseline_v1
python3 research/ic_framework_engines_20260922/compare.py --run baseline_v1
```

`run.py --parts A1-m2-to18,W13-s20260922` runs a subset; several
invocations with disjoint parts may share a run directory, which is how
the frozen run spread its parts over the host's four cores.  The
engine environment overrides (`KIC_F4_*`, `SOLVER_SPLIT_RULE`,
`IC_REDUCTION_CACHE`) must be unset and `run.py` refuses to start
otherwise.

## Comparing a new engine against the frozen baseline

1. Implement `SystemSolver` and register it (`docs/ic/FRAMEWORK.md`).
2. Run the suite with it added to every part, into a new directory:
   `run.py --run-id candidate_<name> --add-engine <name>[:k=v,...]`.
   The frozen engines run again beside it, so the comparison is paired
   on this host rather than against another host's numbers.
3. `compare.py --run candidate_<name> --manifest research/ic_framework_engines_20260922/manifest.json`
   fails if any system fingerprint or any cell's verdict digest differs
   from the frozen run's, and otherwise reports the candidate's paired
   ratios to the baseline and to the reference on every cell, and its
   whole-pipeline `S`.
4. Classify the result by `AGENTS.md` §3 and write it into the ledger
   and the scoreboard in the same change.  A stage ratio is a stage
   diagnostic; only the bench parts' `S` is a speed.

## Results

The frozen run is `results/baseline_v1/`; its `compare.md` grades the
targets of ledger §17.2.  `confirm_v1`, `strength_v1`, `hybrids_v2` and
`hybrids_v3` are the checks declared in §17.3–§17.5, each run from its
own suite file.  The ledger reads them in §17.6–§17.14:

- targets 1 (reach) and 2 (the contract's gate) are met by every
  F4-family engine;
- target 3 (a stage crossing) is met by crossbred against the scalar
  `fes-f2` and confirmed on a holdout, and does not survive the vector
  `fes-f2-wide` (`2.45`–`3.0×`);
- the inherited-F4 hybrid is `1.20×` the vector reference at 34
  unknowns, with a crossing extrapolated near 35;
- nothing measured end to end is below rho.
