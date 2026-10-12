# Validation record

The two measurement sessions pass their native audits: 192/192 executions
verified, including warm-ups; all 96 measured executions replayed with identical
outcomes and deterministic counts. No measured failure or A/A mismatch exists.
The native analyzer's two statistical controls pass in both source versions.

The full command `cargo test --locked --release --lib -j 3` completed with
3,940 passed, nine failed and 112 ignored. It is **not** reported as a green
whole-library run. Six failures were unavailable sparse-checkout fixtures:
the archived source-42 base and five CryptoPro chain tests. Their committed
fixtures were restored; all six targeted retries pass. Two JV-cover F4 tests
failed their valid-system-count thresholds in the full run (21 < 30 and
15 < 25); both pass in a serial retry. The initial process suspension is
retained in PROCEDURE_DEVIATIONS.md rather than removed from the record.

One failure remains reproducible individually:
`cryptanalysis::koblitz_index_calculus::tests::n83_relation_queries_use_full_width_zero_inclusive_additive_draw`,
at `src/cryptanalysis/koblitz_index_calculus.rs:15515`. The first recorded
coefficient is 899991076855513638720967; the test expects
1842027667052925756179639. This branch changes no library source, Cargo.toml
or Cargo.lock relative to base e5a8c7b8cb66cc1d3c332f6ac80c7492482fde41.
The N83 issue remains open; no N83 measurement or correctness conclusion is
made here. This diagnostic PR must not be described as passing every library
test or as fixing that sampler contract.

Logs: `library-tests-missing-fixtures.log`, `library-tests-full.log`,
`library-chain-fixtures-retry.log`, `library-orbit-fixture-retry.log`,
`library-jv-serial-retry.log`, `library-n83-baseline-retry.log`, and
`analysis-tests-v1.log` / `analysis-tests-v2.log`. Graphs and all PDF pages were
visually inspected; static publication checks are recorded before the final
commit.
