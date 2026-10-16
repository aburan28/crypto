# Validation receipt for the inverse-Vélu V2 result

**Workspace:** Darwin arm64 T6041, Rust 1.93.1, release profile;
branch `codex/isogeny-portable-trend-20261015`. The frozen native run
was built at commit `fc1e639d0ac3302b62684419347117faf3e87f64`;
its exact runner, checker, meter, public-input and binary hashes are in
[`RUN_SPEC_V2.json`](RUN_SPEC_V2.json) and
[`RUN_V2/build-sha256.txt`](RUN_V2/build-sha256.txt).

| check | command or record | outcome |
| --- | --- | --- |
| frozen functional run | [`RUN_V2/`](RUN_V2/), [independent replay](RUN_V2/replay.json) | 60/60 children, 1,275/1,275 scalar outputs and 6,126 relation witnesses checked; all 255 inverse/full relation rows matched by the native analyzer |
| raw-file integrity | `(cd RUN_V2 && shasum -a 256 -c ARTIFACT_FILES.sha256)` | 368/368 files matched |
| polynomial certificate | native `different_j_inverse_velu_identity` issuer and `identity_certificate::check_with_extra` | exit 0; [`IDC1h9003f11809fd7a8f`](INVERSE_CUBIC.identity.certificate.json) checks 32 issue plus 32 independent-stream evaluations modulo `2^61-1` |
| Lean proof | `(cd lean && lake build InverseVelu)` | exit 0; `#print axioms` for all three theorems: `propext`, `Quot.sound`; [`LEAN_BUILD_V2.stdout`](LEAN_BUILD_V2.stdout). `lean4checker --fresh` was unavailable. The finite-field curve and fiber-completeness bridge remains open. |
| release library, sandboxed | `cargo test --release --lib -- --test-threads=2` | exit 101: 4,074 passed, 9 failed, 117 ignored, 597.10 s. All nine failures were unrelated `pollard_collab::{cairn,net}` local-socket tests at `TcpListener::bind`, with OS error 1 `Operation not permitted`. This is an environmental suite failure, not a green full-suite receipt. |
| denied socket groups with loopback access | `cargo test --release --lib cryptanalysis::pollard_collab::cairn::tests:: -- --test-threads=2 --quiet`; corresponding `net::tests::` command | exit 0; Cairn 13/13 and net 3/3 passed, including all nine sandbox failures. The rest of the library ran only in the sandboxed full suite above. |
| touched examples | `cargo test --release --example different_j_inverse_velu_oracle --example different_j_inverse_velu_check --example different_j_inverse_velu_analyze --example different_j_inverse_velu_identity -- --quiet` | exit 0; all four compiled and expose zero unit tests. The frozen functional run and independent replay are their behavior checks. |
| source format and lint | `rustfmt --check` on the four touched examples; `cargo clippy --release` on those examples with `-- -D warnings` | both exit 0. Clippy reported the existing repository command-line `clippy::chunks_exact_to_as_chunks` unknown-lint warning in the library and examples; no example-specific diagnostic caused failure. |
| visual reports | `typst compile` on V1 run, V2 run and V2 results; `pdftoppm` preview | all PDFs open, each one page, and diagrams, graph labels and limits were visually checked; benchmark PDF/SVG copies match their study originals byte for byte |

The [V1 split-cap failure](reports/run-20261016T0556Z-v1/REPORT.md)
remains separate from the successful V2 run. The V2 host receipt lacked
a CPU-model string because `sysctl machdep.cpu.brand_string` was denied;
the OS/hardware architecture and toolchain remain recorded. Process CPU
values are shared-host L0 descriptions, and complete portable `S` is
still unset.
