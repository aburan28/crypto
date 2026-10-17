# Validation and retained failures

Date: 2026-10-17. The factor-base construction measurements were made
from the frozen protocol/source commit
`2d3a7f3723925311e74a29ff1cdc7c74a0a97f50`. The later test-only
primary-adapter fixture correction does not change the measured exporter
binary; its SHA-256 is recorded in [RUN_ENV.md](RUN_ENV.md).

| Command or check | Result | Receipt |
| --- | --- | --- |
| `rustc -O support_check.rs -o /private/tmp/n83-support-check && /private/tmp/n83-support-check` | byte-equal to committed `support_check.csv` | [source](support_check.rs), [output](support_check.csv) |
| Native panel manifests, historical hash listing and generic replay | three 54/54 complete manifests; all three hash comparisons exit 0 and all three recovery replays exit 0 | [`runs/`](runs/) |
| `cargo test --release --example koblitz_n83_factor_base_export` before fixture correction | 18 pass, one failure at the zero-trial primary adapter `unwrap`: default pair table cannot pack degree-83 keys | [first example log](example-test.log) |
| `cargo test --release --example koblitz_n83_factor_base_export` after fixture correction | 19 pass, zero fail | [fixed example log](example-test-fixed.log) |
| `cargo test --release --lib` in the workspace sandbox | stopped with `StorageFull` while listing tests after the system volume fell below 300 MiB; it has no final test summary | [first library log](release-lib-test-first.log) |
| `cargo test --release --lib -- --test-threads=1` with SSD temporary directory | interrupted to use four bounded threads; no pass claim | [serial log](release-lib-test-serial.log) |
| `cargo test --release --lib -- --test-threads=4` with SSD temporary directory in the workspace sandbox | 4,207 pass, nine `pollard_collab` local-socket tests fail with sandbox `PermissionDenied`; 117 ignored | [sandbox library log](release-lib-test-four.log) |
| `cargo test --release --lib cryptanalysis::pollard_collab::` with local-socket access | 37 pass, zero fail | [socket test log](pollard-collab-test-escalated.log) |
| Complete release library rerun with local-socket access and four test threads | 4,213 pass, three long `jv_cover` tests fail; 117 ignored. One failure reports a solver wall timeout; the other two have too few uncapped cases under the same long-running workload | [four-thread unrestricted log](release-lib-test-final.log) |
| Three `jv_cover` failures rerun individually with one test thread | all three pass: 18.20 s, 76.62 s and 70.83 s respectively | [status](jv-rerun-status.txt), [Nagao](jv-nagao_finds_planted_six_sums_and_agrees_with_the_oracle.log), [stopped F4](jv-stopped_f4_agrees_with_the_full_solver_and_the_oracle.log), [replayed F4](jv-replayed_f4_agrees_with_the_full_solver_and_the_oracle.log) |
| Complete release library rerun with local-socket access and two test threads | PASS: 4,216 passed, zero failed, 117 ignored; 3,679.41 s | [two-thread library log](release-lib-test-two.log) |

The adapter fixture now explicitly selects `DecompositionStrategy::Enumerate`
for its zero-trial preflight. Its default `PairTable` path returns `None`
for the 83-bit field because its compact key representation is `u64`.
No relation-generation behavior was changed by this test correction.

The PDFs were compiled from the committed Typst sources and visually
inspected. The study and `benchmarks/` PDF, SVG and PNG copies were checked
byte-for-byte. `git diff --cached --check` is the final whitespace gate.
