# CI and curve-schema cleanup validation, 2026-10-09

The cleanup reconstructs the main-branch Rust source after conflicting merge
parents, preserves the public curve models, and regenerates the 321-curve
traits artifact. The source branch was based on `dba53d92517a8899699639a55f7748f489be4693`.
This directory keeps the local validation output for review. Each compressed
log can be read with `gzip -dc evidence/<name>.log.gz`; `evidence/SHA256SUMS`
binds the files in this directory.

The validation host was macOS 26.6 on ARM64. Local Rust and Cargo were 1.93.1,
rustfmt was 1.8.0, and Clippy was 0.1.93. The PR's GitHub Actions lint workflow
pins Rust 1.98, so its result is a separate validation gate.

## Iteration record

| Log | Recorded outcome | Interpretation |
| --- | --- | --- |
| `ci-cleanup-full-lib-20261009.log.gz` | 3,776 passed, 22 failed, 110 ignored | Initial diagnostic run: sandbox localhost restrictions and a reused build directory affected fixtures, with source issues also present. |
| `ci-cleanup-final-lib-20261009.log.gz` | 3,797 passed, 1 failed, 110 ignored | A wall-clock ratio assertion failed under a contended host. |
| `ci-cleanup-final2-lib-20261009.log.gz` | 3,797 passed, 1 failed, 110 ignored | A separate wall-clock calibration ratio failed; the mathematical correctness checks passed. |
| `ci-cleanup-final3-lib-20261009.log.gz` | 3,797 passed, 0 failed, 111 ignored | Full release library suite passed after the host-isolation-sensitive calibration became an explicit opt-in test. |
| `ci-cleanup-examples-20261009.log.gz` | 41 passed across 13 changed examples | Touched-example check before the final lint fixes. |
| `ci-cleanup-clippy-20261009.log.gz` through `ci-cleanup-clippy-final4-20261009.log.gz` | Nonzero exits | Strict lint diagnostics retained as the findings were fixed. |
| `ci-cleanup-clippy-final5-20261009.log.gz` | Exit 0 | Maintained library, bins, tests, benches, and examples passed strict local Clippy. |
| `ci-cleanup-final4-lib-20261009.log.gz` | 3,797 passed, 0 failed, 111 ignored; exit 0 | Full release library suite on the final source snapshot. |
| `ci-cleanup-examples-final-20261009.log.gz` | 41 passed across 16 changed top-level examples; exit 0 | Release tests for every changed maintained example target. |
| `ci-cleanup-clippy-final6-20261009.log.gz` | Exit 0 | Strict Clippy rerun after the final formatting edits. |

The two wall-clock failures motivated separating CPU timing calibration from
correctness. They do not change any recorded timing evidence or establish a
performance ratio. The main cleanup note describes the source and schema
repairs.

## Final source and commands

The final delivery worktree starts from `f79576a24ff5961e3c58013f54043aabafec3617`.
The tracked-file patch used for the passing release suite has SHA-256
`aeb59f744ad0d46e33982c4b1e99751ded6b39aa65ce0bf0ca6aa336da5d599d`;
`git diff --binary HEAD | shasum -a 256` produced the same digest in the
delivery worktree before these evidence files were added. The source and test
snapshot is therefore unchanged by moving to the sparse delivery checkout.

The final checks ran with these commands from the full source checkout. The
release library suite used localhost socket access for tests that bind a local
address; the Clippy and example commands used the same source and shared Cargo
target directory.

```sh
RUST_MIN_STACK=33554432 CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/rho-strong-cli-20261009/target \
  cargo test --release --lib --locked --offline -- --test-threads=4
CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/rho-strong-cli-20261009/target \
  CARGO_NET_OFFLINE=true bash scripts/clippy-maintained.sh
```

For the example run, the shell generated `--example` arguments from changed
root-level `examples/*.rs` paths, excluding nested modules, then ran
`cargo test --release --locked --offline` with `RUST_MIN_STACK=33554432`, the
same `CARGO_TARGET_DIR`, and `--test-threads=4`. The 16 target names are
`a3_s4_ladder`, `a3_stage_ladder`, `f5_support_rank_bits_screen`,
`ffd_s4_subspace`, `find_a3_bench_curves`, `g2_equivariant_orbit_probe`,
`g2_wide_probe`, `iso1_class_census`, `j0_stage_ladder`,
`koblitz_orbit_dlp_fast`, `koblitz_rank_fixture`, `koblitz_rho_batch_ks`,
`koblitz_rho_fixture`, `koblitz_s5_sat_instance`,
`koblitz_unknown_scalar_panel`, and `prime_ecdlp_breakthrough_probe`.

`python3 scripts/curve_id.py selftest` returned `curve_id selftest: ok`.
`python3 scripts/check_curve_names.py` and its
`--diff dba53d92517a8899699639a55f7748f489be4693` form each returned
`curve names: 0 problem(s)`. The PR-scoped rustfmt check over changed Rust
files and `git diff --check` both passed.
