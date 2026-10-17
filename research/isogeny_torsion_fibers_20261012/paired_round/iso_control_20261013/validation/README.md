# Local validation for isomorphism and six-class controls

Host: shared Darwin 25.6.0 arm64, L0; Rust/Cargo 1.93.1, Typst 0.15.1.
The source and binary of the accepted v3 runner have SHA-256
`8b37c4953ec9f65195ad7564eaedfb4fa2c8b1c0d9746b6bb811f7c4b48e7a0a`
and
`2cdd8d72d52e51ab38f998a590ab30599ed32bf1edddbaa3d8799da6cd55e21f`.

| check | command or receipt | result |
|:--|:--|:--|
| release example build and substantive v3 run | `cargo build --release --example isogeny_fiber_iso_control_20261013`; run the release binary to `raw-v3.tsv` and `run-v3.log` | passed; 80/80 verified arm records, three charged six-class restarts |
| deterministic replay | run the same binary to `replay-v3.tsv` and `replay-v3.log`; `cmp` each with its first-run counterpart | both passed byte for byte; raw/replay SHA-256 `8456a90d17d1b67c6817943e5a4e5139541c33d88a7ab7626c06d4cd638f834c` |
| independent arithmetic and accounting | `rustc --edition 2021 -O -D warnings .../verify.rs -o /private/tmp/isogeny-iso-control-verifier`; run on v3 raw, frozen v4 raw, and v3 log | passed; 80 records, 16 source targets, five arms, three restarts, 89,764 cubic-action pairs, 44,882 scaling pairs, 1,840 repeated six-class point checks; [receipt](../independent_v3.log) |
| historical failed-run replay | rebuild reconstructed v1/v2 native sources with the release library artifact; compare stdout/stderr with original files | passed byte for byte with expected exit 2; v1 source hash exactly matches pre-run receipt; v2 original source hash unavailable; [details](historical_replay.md) |
| derived table and graph replay | rebuild `plot.rs` with `rustc -O -D warnings`, regenerate into `/private/tmp`, and `cmp` summary, phase summary, relation SVG, rho SVG | all four `cmp` checks passed byte for byte |
| Rust formatting | `rustfmt --check` on the example, checker, plotters, and retained runner sources | passed; [log](rustfmt.log) |
| touched example test target | `cargo test --release --example isogeny_fiber_iso_control_20261013` | passed; zero `#[test]` functions, so substantive checks are the run and independent verifier; [log](example_v3.log) |
| required release library suite | `RUST_TEST_THREADS=2 cargo test --release --lib` with loopback access outside the restricted sandbox | passed: 4,009 passed, zero failed, 113 ignored; test execution 1,046.86 seconds; [full log](release_lib_v3.log) |
| Lean conditional group-action component | `lake build` in `research/isogeny_torsion_fibers_20261012/lean` | passed; [receipt](../../../lean/six_class_build.log) includes `#print axioms`; `lean4checker --fresh` unavailable |
| PDF render and inspection | Typst 0.15.1 compiled result, transport diagram and v1/v2/v3 run reports; Poppler rendered every page | passed; result PDF has three pages including the explicit transport diagram, run PDFs have one, one, and two pages |

The full release library suite was run after the v3 example source was
finalized. No library or v3 example source changed during it. The
historical v1/v2 sources were reconstructed afterward without changing
the accepted v3 executable. No Python cryptographic research path or
relevant Python suite was changed by this follow-on.

The verified v3 runner's `mapped_qx/y` field remains arm-specific, as the
[report](../REPORT.md) explains. The independent checker validates its
meaning for each arm and checks every answer against the source target.
