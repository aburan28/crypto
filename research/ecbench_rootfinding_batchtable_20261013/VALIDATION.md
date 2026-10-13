# Validation receipt before the current-main merge

The frozen spec was accepted by `target/release/ecbench plan` and expanded
to 128 executions on eight registered public workloads. The session was
created from clean source commit `3904a4f9d4ad0def62c0a4db9bf847543e770eaa`
with release binary SHA-256
`4c55e01badbba18a43c7b096ff42c4e701996cc5311ce708a295b3235b4db391`.
It completed with 128/128 verified executions. The [pair
checker](pair-check.json) reports 32/32 exact paired outputs and counted
phases. The [full replay](replay-all.json) reports 96/96 identical
measured replays and zero problems, SHA-256
`fafd4865ea91fbf2e50a46a470a8fbf5108513fa2186342cd4bf7344322ab0db`.

The focused source test
`cargo test --release --lib prime_batch_table_matches_scalar_table_and_reduces_inversions -- --nocapture`
passed 1/1; its [raw log](focused-test.log) has SHA-256
`9ff70a17fc2bd303788b2a76e22544d9abfdb8e4ffa9f3a98707208633159ada`.
The same test and the current-source
`rootfind_subtract_follows_subtract_query_stream_on_a_prime_curve` test
also passed in the full suite.

The required `cargo test --release --lib` exited 101 after 774.86 seconds:
**4,045 passed, 9 failed, 114 ignored**. All nine failures are
`cryptanalysis::pollard_collab` socket/listener tests whose raw errors are
`PermissionDenied: Operation not permitted` under the default sandbox.
The [full log](full-lib-test.log), SHA-256
`626ac667ca76c57c08fa6cb292b878916061c5022abd70558611d3a450c55e89`,
is retained. The focused loopback-enabled command
`cargo test --release --lib cryptanalysis::pollard_collab:: -- --test-threads=1`
passed **37/37** with zero failures; its [raw log](pollard-loopback-test.log)
has SHA-256
`bba13396827b5ae14907b6478a058eae5f378f5ce0ce38398f589e92e496b70b`.
This classifies the nine sandbox failures; it is not a clean whole-suite
pass receipt.

Both result and run PDFs compile to one page and were visually inspected.
The graph PNG is a 150 dpi render of the same SVG source; the benchmark
catalog PDF/SVG/PNG copies were byte-compared with the study originals.
`git diff --check` passed. All measured wall records are L0; no wall-time
claim was admitted.

The branch will merge `origin/main` before PR delivery. That integration
requires a separate check on the merged revision; this receipt binds the
pre-merge source and is not silently reused as a post-merge test result.
