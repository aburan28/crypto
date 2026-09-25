# Generic query accounting controls

Classification: **accounting**. The generic worker now retains every attempted
collection/descent query and its frontend result, including unsuccessful terminal
descents. No performance comparison or improvement round was run. One bounded
improvement round has completed without a qualifying winner; two remain.

The observed APIs run the same query sequence once, retain native counters, and
sort parallel collection records by trial. Existing callers may keep using
`collect` and `solve`; `solve_report` additionally returns unsuccessful terminal
reports, and `solve_observed` includes their complete query history. A witness
remains subject to independent group verification. Pair-table/window misses are
unresolved, encoding rejections are unsupported, and incomplete searches are never
promoted to proved negative answers. The worker exports the persistent LA report
and counts actual matrix solve attempts, excluding checks with too few rows.

Local release controls passed on macOS arm64 / Rust 1.93.1 with the repository's
autolab dependency lock and one test thread:

| Control | Result | Evidence |
|---|---|---|
| New observed-query Rust controls | 4/4 passed | [raw log](rust-controls.log) |
| Descent-filter regressions (includes one new query test) | 63/63 passed | [raw log](existing-regressions.log) |
| Unsupported-input and failed-LA regressions | 5/5 passed | [raw log](existing-regressions.log) |
| Collected-relation/precompute equivalence | 1/1 passed | [raw log](equivalence.log) |
| Python harness tests | 126 passed, 2 platform skips | [raw log](python-tests.log) |
| Actual generic worker controls | 28 complete solves and 7 intentional incomplete runs passed | [summary](worker-summary.json), [all raw reports](worker-raw.jsonl) |

The Rust commands executed 73 passing tests, representing 72 distinct tests.
Local Clippy also completed successfully, with existing repository/toolchain
warnings retained in [its log](clippy.log).

The worker controls cover all seven built-in backends (pair table, enumeration,
F4, F5, inherited F4, native-XOR SAT and CNF SAT), dense/sparse final LA, and both
Koblitz degree-9 curve coefficients. The independent Python group checker verifies
each witness, scalar recovery and factor-base census. Its bounded two-summand
controls also enumerate all pair sums to reject false negative answers. Mutation
controls reject missing/reordered queries, forged coefficients/witnesses, lost
terminal descent history, false recovered scalars and incomplete UNSAT claims.
The query-law RNG and general matrix/stage auditor are not certified by this
checker. External WDSat/MQ-FES and Weil-chart performance are not qualified here.

One first regression command used a filter matching zero tests. Its output is
retained; it is not counted as validation. The corrected exact-name command
executed the equivalence test successfully. No worker or solver control failed.
The Linux integration workflow now repeats the query controls, replays native and
profiled records, and retains seven expected incomplete worker runs alongside its
34 existing integration cases. Its CI result is recorded on the PR.

Reproduce the local worker controls from the repository root, in a new directory:

```sh
cp research/ic_candidate_tournament_20260915/ci/Cargo.lock Cargo.lock
cargo test --locked --release --lib query_accounting_ -- --test-threads=1
cargo build --locked --release --example ic_tournament_worker
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/generic-query-accounting/run_controls.py --worker target/release/examples/ic_tournament_worker --out /tmp/ic-query-controls-fresh
python3.12 -m unittest discover -s research/ic_candidate_tournament_20260915 -p 'test_*.py' -v
```

The [source/receipt manifest](controls.json) binds the implementation and retained
local evidence; [inputs](worker-inputs.json) fix every job. Local elapsed values
are raw correctness-run diagnostics. Online wall time, complete cold cost, Ir,
S, rho ratio and speedup remain **unknown**, with `promotion_eligible=false`.
Public-point input, exclusive scientific phases and full generic admission remain
required before these pipelines enter a comparative tournament.
