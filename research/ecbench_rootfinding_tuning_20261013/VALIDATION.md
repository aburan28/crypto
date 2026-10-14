# Validation receipt

No Rust source was changed in this iteration. The frozen spec was accepted by
`target/release/ecbench plan`, expanding to 128 executions on eight public
workloads, and the run completed with 128/128 verified outcomes. The saved
[`replay-all.json`](replay-all.json) reports 96/96 identical measured replays,
zero problems, and an audit receipt SHA-256 of
`cdd1a0bcde5fca4939b912688f20ab662f53b3557197c742f3f01788b8135c58`.
The new cohort's eight target values have zero exact overlap with the prior
2026-10-12 prime cohort's eight target values.

The required `cargo test --release --lib` command ran against the unchanged
library and exited 101 after 612.86 seconds: **3,999 passed, 11 failed, 113
ignored**. The exact output is retained in
[`cargo-test-release-lib.log`](cargo-test-release-lib.log), SHA-256
`cf984987c0a6ac62d2736c56ee97f4e8d6720ac7291598c03a8fca4c22e88383`.
Two `jv_cover` solver tests failed in the heavily concurrent suite. Nine
`pollard_collab` tests failed when sandbox policy denied local socket or
listener operations (`Operation not permitted`). This is not a clean
whole-suite pass.

Focused reruns, without source changes:

| Command | Result | Raw log |
| --- | --- | --- |
| `cargo test --release --lib cryptanalysis::jv_cover::tests:: -- --test-threads=1` | 14 passed, 1 ignored, 0 failed | [`cargo-test-jv-cover.log`](cargo-test-jv-cover.log), SHA-256 `1ba81b04340d5adfc1bf1bf80c78f4fa5d0e668f9d23a82ae392d7a93c80b5ec` |
| `cargo test --release --lib cryptanalysis::pollard_collab:: -- --test-threads=1` with local loopback permission | 37 passed, 0 failed | [`cargo-test-pollard-collab.log`](cargo-test-pollard-collab.log), SHA-256 `6c88d77d7377bbb01ef2c6333b7c0d866dfd6805894a070e7f94ffe953ed43d5` |

Both previously failing `jv_cover` tests and all nine previously failing
`pollard_collab` tests passed in those focused reruns. The benchmark session
and all reported numbers were produced before these reruns and were not
altered by them. The results and run PDFs each compile to one page and were
visually inspected. The benchmark catalog PDF/SVG copies were byte-compared
against the study originals.
