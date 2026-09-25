# Relation-LA accounting controls

Both new release controls passed, in dense and sparse modes, and the existing
collected-relation/precompute equivalence test passed. The failure control keeps
two singular-matrix attempts in the stored report; the success control keeps
the returned and stored report consistent and verifies every column log in the
group. A check with too few rows remains zero attempts.

The change fixes lost bookkeeping in `FactorBaseLogSolver::try_solve`. It does
not change the matrix solver, relation generation or success criterion. Failed
attempts now update the persistent counters, elapsed LA time and sparse report
before returning `None`; successful attempts persist their full report as well.

These [protocol-defined](PROTOCOL.md) tests used macOS arm64, Rust 1.93.1,
the checked-in autolab lockfile and one test thread. Exact source/log hashes and
null performance fields are in [controls.json](controls.json), with raw output
in [accounting.log](accounting.log) and [equivalence.log](equivalence.log).
The homogeneous matrix is explicitly a unit fixture, not collected IC evidence.
The success control uses actual verified toy relations from 600 ordinary queries
per LA mode. No time, operation count or natural-yield estimate from these
controls is used to rank candidates. Linux integration is the PR acceptance gate.

```sh
cp research/ic_candidate_tournament_20260915/ci/Cargo.lock Cargo.lock
cargo test --locked --release --lib la_accounting_ -- --test-threads=1
cargo test --locked --release --lib logs_from_collected_relations_match_the_single_process_precompute -- --test-threads=1
```

This follow-on depends on PR 789's PDP outcome correction. Both are admission
repairs for the generic path; the measured optimized round-one source and its
negative result remain unchanged. Classification: **accounting/correctness**.
Complete IC online/cold costs, Ir, S and speedup are unknown for this change.
The goal remains active with two improvement rounds available. Generic per-query
outcomes, terminal failed-descent reports and exclusive phase instrumentation
remain required before scoring its F4/F5/SAT pipelines.
