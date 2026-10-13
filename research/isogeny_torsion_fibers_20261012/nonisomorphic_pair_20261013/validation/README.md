# Local validation: non-isomorphic F617 pair

Host: shared Darwin 25.6.0 arm64, L0; Rust/Cargo 1.93.1, Typst 0.15.1.
Run dependency revision: `fa1a6d4dee12fec5e5f2d7b3b211520b9e2b2bf3` plus
the source snapshot named below. Accepted release-runner source SHA-256:
`72b25721734346d82e4268c54622251312b050aa2f23e573ec52616e165d60b2`;
release binary SHA-256:
`cc28892bcaa8830ddc6d1c643088c715c7ad4e0c6d097a5c40d618952c277cb6`;
v2 raw TSV SHA-256:
`f17c47549eb02abfaf3fec1cf30d4181d34a491d274bfcdf06eca64f5561e68e`.
The output was produced after the [protocol](../PROTOCOL.md) and
[wrapper remediation](../REMEDIATION_V2.md) were committed.

| check | exact command or receipt | outcome |
|:--|:--|:--|
| native release build | `cargo build --release --example isogeny_fiber_nonisomorphic_20261013` | passed; existing library warnings recorded by Cargo, no new example warning |
| structural preflight and substantive v1 child | `/usr/bin/time -l target/release/examples/isogeny_fiber_nonisomorphic_20261013 > .../raw-v1.tsv 2> .../run-v1.log` | child emitted 24 verified rows and final summary; outer wrapper exit 1 on denied `sysctl kern.clockrate`; no admitted wall time |
| clean v2 direct execution | `target/release/examples/isogeny_fiber_nonisomorphic_20261013 > .../raw-v2.tsv 2> .../run-v2.log` | exit 0, 24/24 verified arm records |
| deterministic replay | same binary to `replay-v2.tsv` and `replay-v2.log`; `cmp` each with v2; `cmp raw-v1.tsv raw-v2.tsv` | all three byte comparisons passed |
| independent arithmetic and relation replay | `rustc --edition 2021 -O -D warnings verify.rs -o /private/tmp/isogeny_noniso_verifier`; run verifier on `raw-v2.tsv` with output `fibers-v2.tsv` | passed; 24 scalar records, 16 independently reconstructed relation solves, 8 exact [target-fiber triples](../fibers-v2.tsv), 381,924 homomorphism pairs, 3-point kernel, 206-point image, exact phase/map counters; [receipt](../independent_v2.log) |
| exact model identities | compile [identity.rs](../identity.rs) against the release `crypto_lib` and matching `num_bigint` rlibs; compare fresh output with [identities.tsv](../identities.tsv) | passed byte for byte |
| derived tables and graph | `rustc --edition 2021 -O -D warnings summarize.rs -o /private/tmp/isogeny_noniso_summarize`; regenerate v2 summary, phase, paired-target table, and SVG in `/private/tmp`; `cmp` with study | four byte comparisons passed |
| Rust formatting | `rustfmt --check` on runner, verifier, summarizer and identity probe | passed; [log](rustfmt.log) |
| Lean finite certificate | `lake build Nonisomorphic617`, then `lake build` in the pinned [Lean project](../../lean/lakefile.toml) | both passed; [theorem receipt](../../lean/nonisomorphic617_build.log) prints no axioms for the sixth-power and coprimality theorems, and standard `propext` for the scalar-kernel implication; [full project receipt](../../lean/nonisomorphic617_full_build.log) passed; `lean4checker --fresh` unavailable |
| touched example test target | `cargo test --release --example isogeny_fiber_nonisomorphic_20261013` | passed; zero `#[test]` functions, so substantive checks are the native solve and independent verifier; [log](example_test.log) |
| first release library suite | `RUST_TEST_THREADS=2 cargo test --release --lib`, with local loopback access | exit 101 after 3,889.89 s: 4,007 passed, 2 failed, 113 ignored; the unchanged `jv_cover` Nagao test hit its 120-second F4 budget, and the replayed F4 test completed 29 eligible cases against its threshold of 30; [full failure log](release_lib.log) retained |
| isolated checks of the two failures | `RUST_TEST_THREADS=1 cargo test --release --lib <exact test name> -- --exact --nocapture` for each `jv_cover` test | both passed under lighter host load: [Nagao](jv_nagao_isolated.log) 1/1 in 5.57 s; [replayed F4](jv_replayed_isolated.log) 1/1 in 34.89 s, with 40 eligible cases; this supports, but does not prove, host contention as the earlier cause |
| required release library rerun | `RUST_TEST_THREADS=2 cargo test --release --lib`, with local loopback access | exit 0: 4,009 passed, 0 failed, 113 ignored in 897.27 s; [clean full-suite log](release_lib_rerun.log) |
| visual render | Typst compiled result and v1/v2 run PDFs, plus graph/diagram PDFs; Poppler rendered every page and the result and run pages were inspected | passed; result PDF 2 pages; v1/v2 run PDFs 1 page each; parent addendum PDF 3 pages |

The native runner counts every group addition and doubling in the target
solves; its curve-specific discovery and exhaustive map preflight are not
charged per target. The map, factor-base square roots, row operations,
lookups and memory have separate counters or documented gaps, so these
checks do not admit a complete native-cost speed claim. The replay and
arithmetic checker used the same host; independent-host validation remains
open. No Python cryptographic research or performance path was introduced.
