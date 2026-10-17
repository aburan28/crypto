# Pairing-transfer validation

The cryptographic execution path, accounting, and scalar replay are native
Rust. The [V1](reports/run-20261015T0546Z/RUN_REPORT.md) and
[V2](reports/run-20261015T0555Z/RUN_REPORT.md) reports give the exact run
commands, frozen hashes, per-command status, and raw records. V1 and V2 each
completed 30/30 commands with zero exit status and empty command stderr.
The separate verifier checked all six recovered source scalars; the extra
pairing checks passed 18/18; both rho controls passed 12/12. The six entries
in [`validation/v1_v2_equivalence.json`](validation/v1_v2_equivalence.json)
confirm that V1 and V2 retained identical answers, base-field operation
counts, extension multiplication counts, and rho group-addition equivalents.
This equivalence concerns deterministic values, not timing.

## Local gates

| Gate | Exact command or receipt | Outcome |
|:--|:--|:--|
| Release library | `TMPDIR=<study>/validation/tmp RUST_TEST_THREADS=1 cargo test --release --lib --offline` | Exit 0: **4,083 passed, 0 failed, 117 ignored**, 1,741.89 s. [`release_lib_workspace_serial.log`](validation/release_lib_workspace_serial.log). The temporary directory was on the workspace volume. |
| Touched release examples | `TMPDIR=<study>/validation/tmp cargo test --release --example torsion_pairing_transfer --example torsion_pairing_accounting --offline` | Exit 0; both examples compiled and their test targets ran. They contain zero embedded unit tests; the frozen command runs and native accountant assertions above are the functional checks. [`touched_examples_final.log`](validation/touched_examples_final.log). |
| V2 native accountant replay | `cargo run --release --example torsion_pairing_accounting --offline -- <study> runs_v2` | Exit 0; its fresh [`accountant_replay.json`](validation/accountant_replay.json) is byte-identical to [`SUMMARY_V2.json`](SUMMARY_V2.json) by `cmp` exit 0. Compiler and command stderr: [`accountant_replay.stderr`](validation/accountant_replay.stderr). |
| Raw record hashes | `shasum -a 256 -c <manifest>` on [`SHA256SUMS`](SHA256SUMS) and [`SHA256SUMS_V2`](SHA256SUMS_V2) from the repository root | Both exit 0, 62 checked paths each. [`V1 log`](validation/v1_sha256_check.log), [`V2 log`](validation/v2_sha256_check.log). |
| Formatting and staged diff | `rustfmt --check --edition 2021 examples/torsion_pairing_transfer.rs examples/torsion_pairing_accounting.rs`; `git diff --cached --check` | Both exit 0. |
| Report renders | Typst compilation plus `pdfinfo`, `pdftoppm`, and visual inspection of rendered pages | Results PDF has three A4 pages; each run report PDF has one page. Every SVG and 150 dpi PNG benchmark copy matches its study source byte for byte. |

Create `<study>/validation/tmp` before either command that sets `TMPDIR`;
the suite's routine generated scratch file was removed after the completed
run. The exact study path is this file's parent directory.

The release build and runs used Rust 1.93.1 on Darwin 25.6.0 arm64. The
pairing measurement itself was shared-host L0; the release test is a
correctness gate and does not upgrade measurement isolation. No Python
cryptographic path or Python test suite was changed.

## Earlier attempts retained

The first release-library attempt was intentionally interrupted with exit
130 while V2 instrumentation was still being added; its partial output is
[`release_lib_sandbox_interrupted.log`](validation/release_lib_sandbox_interrupted.log).
A subsequent unsandboxed two-thread run exited 101 before reaching a
terminal test summary. Its log ends with the incomplete text `error`, after
one `jv_cover` test reported a failure; the system Data volume had only
106 MiB available when inspected. The exact termination cause is not
established from that truncated output. The unmodified
`cryptanalysis::jv_cover::tests::nagao_finds_planted_six_sums_and_agrees_with_the_oracle`
test then passed when run alone with workspace-backed scratch
([`jv_cover_standalone.log`](validation/jv_cover_standalone.log)) and passed
again in the complete serial release run. The incomplete two-thread record
is preserved as [`release_lib_tmp_full_failure.log`](validation/release_lib_tmp_full_failure.log),
not counted as a passing full-suite run.

Neither the interrupted nor the two-thread attempt changed the frozen V1
or V2 measurement files. The complete serial run supplies the release gate
used for this PR. Study inputs, run records, reports, and validation logs
written before manifest verification are sealed by
[`ARTIFACTS.sha256`](ARTIFACTS.sha256). The manifest verification log is a
derived receipt and is intentionally outside that manifest to avoid a hash
cycle.
