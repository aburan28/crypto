# WDSat child-status and empty-output correctness gate

Status: pre-change protocol, 2026-09-25, based on `main` ec990d85b61eacccefa9dee97f2a835e7a20b785. This is a triage repair, not a new solver-performance experiment. The frozen `research/wdsat_rr_20260920/` Riemann–Roch panel uses the separate Python `ecc2k130/codegen/wdsat.py` exporter/runner, not this Rust Semaev adapter. Its published outcomes and prior limitations are unaffected by this code change; no historical result without stored child status is retrospectively reclassified. No new index-calculus boundary is claimed.

## Defect and upstream contract

`src/cryptanalysis/wdsat_oracle.rs::run_wdsat` currently drops `Output.status`. Then `wdsat_decompose` treats an empty stdout as a completed UNSAT result. A child that exits nonzero before printing anything can therefore be counted as a refutation. Upstream WDSat revision [`61c6ff3f49445af2729274819edb27b85a89efc1`](https://github.com/mtrimoska/WDSat/tree/61c6ff3f49445af2729274819edb27b85a89efc1) prints `UNSAT` for exhausted search and `UNSAT on XORGAUSS init` for an initialization refutation in [`src/wdsat.c`](https://github.com/mtrimoska/WDSat/blob/61c6ff3f49445af2729274819edb27b85a89efc1/src/wdsat.c#L462-L561). Its [`src/main.c`](https://github.com/mtrimoska/WDSat/blob/61c6ff3f49445af2729274819edb27b85a89efc1/src/main.c#L113-L145) normally exits successfully after invoking the solver. Neither source establishes a successful *empty* output as an UNSAT certificate.

## Frozen classification rule

| Child result | Model output | Result |
|:--|:--|:--|
| Any nonzero exit or signal | Any | unknown/exhausted; never refuted, including if stdout contains an UNSAT word |
| Zero exit | Parseable model that passes ANF, original-equation and group-law checks | SAT relation |
| Zero exit | Exact upstream `UNSAT` or `UNSAT on XORGAUSS init` stdout line, with no model | refuted |
| Zero exit | Empty stdout or no recognized model/status line | unknown/exhausted |
| Timeout, capacity, spawn or output-collection failure | Any | unknown/exhausted |

A model that fails any independent check remains spurious/exhausted. The adapter must not infer UNSAT from stderr diagnostics or an incidental substring. Preserve `solver_calls = 1` for a spawned child and existing capacity/timeout semantics.

## Bounded verification before merge

1. Add a mock-child regression for zero and nonzero empty stdout, explicit successful UNSAT, nonzero exit with misleading UNSAT text, and successful model output. Check both the process boundary and oracle verdict; use a tiny frozen Koblitz fixture if the full oracle is needed.
2. Re-run the existing WDSat ANF/parse/validation unit tests and the pinned planted n=7 WDSat/native-SAT witness comparison where the existing WDSat binary is available. If the binary is unavailable, report that test as unavailable rather than substitute a mock for actual solver agreement.
3. Replay the existing Rust ANF/parser fixtures and inspect archived WDSat corpus receipts without changing them; check hashes when stored. Distinguish the separate Riemann–Roch Python panel from this Rust Semaev path. Do not infer a changed historical verdict where child status was not stored, and do not rewrite timing rows or the canonical scoreboard from a parser fix alone.
4. Keep the change to child status and verdict triage plus focused tests/docs. Any subsequent WDSat comparison needs a new frozen corpus, raw stdout/stderr/exit receipts, and independent witness checks.
