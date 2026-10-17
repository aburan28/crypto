# Validation and failure receipts for the order-229 batch

The [V2 prefix correction](PROTOCOL_V2_PREFIX_CORRECTION.md), final
[source/binary freeze](FREEZE_SPEC_V1.json), [public target
freeze](RUN_SPEC_V1.json), and their commits preceded the 21:14:48 UTC
measurement. The frozen run produced 84/84 successful child processes;
the separate native [checker](REPLAY_V1.json) replayed 1,785/1,785
scalars and 13,141 relation witnesses. The checker did not reuse the
runner's point arithmetic or relation-elimination functions. The
solver had no access to the checker scalar input.

| Check | Exact command or action | Exit and receipt |
| --- | --- | --- |
| Controlled release library | `TMPDIR=target/research-tmp cargo test --release --lib -- --test-threads=2` in a network-enabled context | Exit 0; 4,083 passed, zero failed, 117 ignored in 1,026.37 s. [stdout](../different_j_order_screen_20261015/VALIDATION_LIB_V2.stdout) SHA-256 `93e8fd2e0d3e5f77a3d95a43773ec8d33c4d7139ef0508012795cd2082cdc397`; [stderr](../different_j_order_screen_20261015/VALIDATION_LIB_V2.stderr) SHA-256 `0c73c36a0c996bcb10fe1c94315f9ecaaa06a57ebe458a945903a1d5efac2b9d`. |
| Four touched release example harnesses | `CARGO_BUILD_JOBS=2 CARGO_TARGET_DIR=target/order-batch-check cargo test --release --example different_j_order_batch --example different_j_order_batch_check --example different_j_order_process_meter --example different_j_order_batch_summary` | Exit 0; all four compiled and reported zero built-in tests. [stdout](VALIDATION_EXAMPLES_V1.stdout) SHA-256 `52e64529d3fe7bfb8771e2beee8782ae3171424375e92b1995b3d506ddb93a52`; [stderr](VALIDATION_EXAMPLES_V1.stderr) SHA-256 `e11d1f3672e24ede7e67788da5249c0b0133633bda0303716a047d213048c5a6`. Their executable correctness paths are checked by the full native run/replay. |
| Touched-example Rust Clippy | `TMPDIR=target/research-tmp rustup run 1.98 cargo clippy --release --example different_j_order_batch --example different_j_order_batch_check --example different_j_order_process_meter --example different_j_order_batch_summary` | Exit 0. [stderr](VALIDATION_CLIPPY_V5.stderr) SHA-256 `a60fd07cb68cc53400aed2f7e179899ae37d14498a6df87f039ed596ab9533c0` has only the inherited `clippy::chunks_exact_to_as_chunks` unknown-lint warning. |
| Format and shell syntax | `rustfmt --check` on the four new examples; `zsh -n` on `freeze_v1.zsh` and `run_v1.zsh` | Both exit 0. |
| Run audit | `zsh run_v1.zsh` after its hash checks, then the independent Rust checker; `different_j_order_batch_summary` on the frozen raw/checker files | All exit 0. [Run](RUN_V1.jsonl) SHA-256 `d0e4616236aeb8978f5ac1549dde373be294bdf52503c954a6f964904f68faf8`; [meters](METERS_V1.jsonl) `8cbb7e5d5737792231fa703d1bdf35535230132bff24062746abced60daaf776`; [replay](REPLAY_V1.json) `c4572cc1dda31ddbfaec37190433eac60e8862769f77333de0d6ea6940d92037`; [summary](SUMMARY_V1.json) `ad4b969aec7e6cd5aebf7f7d2ca0dfd851d22a96acfaef492838fc194d4e8af5`. |
| Visual render and mirror | `typst compile` for the run and results PDFs; `typst compile --format png --ppi 150` for the SVG wrapper; `pdfinfo`, page render and visual inspection; `cmp -s` of both PDFs and SVG/PNG against `benchmarks/20261015-different-j-order-batch/` | All exit 0. Both PDFs have two A4 pages; the full 1000×535 SVG graph, tables, IDs and limits are legible. The four mirror files are byte-identical. |

The first library attempt ran inside a sandbox that denied loopback
sockets and recorded **4,074 passed, 9 failed, 117 ignored** in
[V1 stdout](../different_j_order_screen_20261015/VALIDATION_LIB_V1.stdout)
and [V1 stderr](../different_j_order_screen_20261015/VALIDATION_LIB_V1.stderr).
The network-enabled V2 run passed all 4,083 active tests; the V1
failures were environment permission failures, not a mathematical
result. A compile-check V1 failed on a missing `CountedGroup` trait
import and remains in [CHECK_V1.stderr](CHECK_V1.stderr); the
corrected [CHECK_V2](CHECK_V2.stderr) and final release build passed.
The first Clippy launcher used unsupported `cargo +1.98`; the second
strict attempt stopped on inherited `unusual_byte_groupings`; the
third strict attempt stopped on inherited unrelated all-target
warnings. Their versioned logs remain in the [screen
directory](../different_j_order_screen_20261015/) and were not
overwritten by the successful default Clippy V4/V5 runs.

The final source and binary hashes for the runner, checker, meter and
summary are in [RUN_SPEC_V1.json](RUN_SPEC_V1.json). The run was L0
and has only same-host arithmetic replay; a second-host performance
receipt and a portable ecbench primitive conversion remain open.
