# Validation receipts for the 60-target reuse result

The native V3 and V4 runs preceded these repository checks. Their
answer-free solver children and separate replay processes are described
in the [run report](reports/run-20261015T1955Z/REPORT.md); each version
checked 1,260/1,260 recovered scalars and 10,960 relation witnesses.

| Check | Exact command or action | Exit and receipt |
| --- | --- | --- |
| Controlled library | `TMPDIR=target/research-tmp cargo test --release --lib -- --test-threads=2` | Exit 0; 4,083 passed, zero failed, 117 ignored in 1,017.42 s. [stdout](VALIDATION_LIB.stdout) SHA-256 `e2727e1b6e3ec6a8c14ada08298362ec26c27e60ddcfeca91153cc3a6b906444`; [stderr](VALIDATION_LIB.stderr) SHA-256 `8080e57368cb131b80e02068fa774fa064230ba1982df6dff44acc353bed2440`. |
| Six touched example harnesses | `TMPDIR=target/research-tmp cargo test --release --example nonisomorphic_extension_isogeny --example different_j_process_meter --example different_j_reuse_summary --example different_j_process_summary --example chi_polynomial_decomposition --example torsion_pairing_batch -- --test-threads=2` | Exit 0; six harnesses compiled and each reported zero built-in tests. [stdout](VALIDATION_EXAMPLES.stdout) SHA-256 `88605fcde34d35afd53808b1d5e2253ccbdcbbf80f59515b05e2b9c5226303c7`; [stderr](VALIDATION_EXAMPLES.stderr) SHA-256 `dfc908ef835e11b581e4464d4270f80c0018c3ae3e96c3eb3d474b7c4c5da4d7`. Their executable research paths were checked by the separate native run/replay receipts, not by these empty harnesses. |
| Rust format | `rustfmt --check examples/nonisomorphic_extension_isogeny.rs examples/different_j_process_meter.rs examples/different_j_reuse_summary.rs examples/different_j_process_summary.rs examples/chi_polynomial_decomposition.rs examples/torsion_pairing_batch.rs` | Exit 0. |
| Rust 1.98 touched-example Clippy | `TMPDIR=/Volumes/SSD990-2/.codex/worktrees/d92d/crypto/target/research-tmp rustup run 1.98 cargo clippy --example nonisomorphic_extension_isogeny --example different_j_process_meter --example different_j_reuse_summary --example different_j_process_summary --example chi_polynomial_decomposition --example torsion_pairing_batch -- -D warnings` | Exit 0. [stdout](VALIDATION_CLIPPY_V3.stdout) was empty; [stderr](VALIDATION_CLIPPY_V3.stderr) SHA-256 `d09a8c8f0cff8ae0ab9714f6d8c9838fe28d728233175b9759fd8fe050a591fc` includes an inherited unknown-lint warning from the repository's Clippy flag and the successful finish line. |
| Visual mirror | `cmp -s` for each source-study PDF, SVG and PNG against the six matching files in `benchmarks/20261015-different-j-reuse/` | All six byte comparisons exited zero. The report PDF has two A4 pages; the results PDF has one A4 page. |

Two Clippy launcher attempts are retained as infrastructure failures.
The first used Homebrew `cargo +1.98`, which exited 101 because that
`cargo` does not understand a rustup toolchain directive; see
[stderr](VALIDATION_CLIPPY.stderr) SHA-256
`f3fd728f94473137d3b8d97975046f6146a30b0230a9a035641cacad7942e8eb`.
The second used `rustup run 1.98 cargo` but a relative `TMPDIR` that
dependencies resolved below their own source directories; it exited
101 with missing temporary directories. Its [stderr](VALIDATION_CLIPPY_V2.stderr)
has SHA-256 `ee84fbbef1742dc338e7ec65f9cba85218f762542074d9cc75f5fed5bdb6f380`.
The V3 command fixed both launcher conditions; neither earlier exit
is a cryptographic or source-code test failure.

The six example harnesses were checked after the lint-only edits to
the χ and pairing examples. The research measurements still refer to
their historical hash-pinned binaries, as recorded in the
[post-run source note](../POST_RUN_LINT_FIX_20261015.md).
