# Different-j full-image V2 validation

All commands ran in `/Volumes/SSD990-2/.codex/worktrees/d92d/crypto`
on 2026-10-15 unless a working directory is stated. The release library
suite ran after the native measurement, so its CPU load did not overlap
the measured V2 comparison.

| Check | Command and outcome | Receipt |
| --- | --- | --- |
| Controlled release library suite | `cargo test --release --lib -- --test-threads=2`: exit 0; 4,083 passed, 0 failed, 117 ignored; 957.05 s. | [stdout](validation_release_lib_20261015.stdout) SHA-256 `bf12cfd7033739e3f0a1632293633cb64670a40dd1ca6e37184d11a1d722245d`; [stderr](validation_release_lib_20261015.stderr) SHA-256 `b5185875c17aaa07245ef6d16f2bea1bc5ca0d0e7dcb41bf615dacd379eeaca8`. |
| Touched example harnesses | `cargo test --release --example nonisomorphic_extension_isogeny --example different_j_full_image_summary --example different_j_support_census`: exit 0; three compiled harnesses each have zero built-in tests. The native run/replay below are their substantive correctness check. | [stdout](validation_touched_examples_20261015.stdout), [stderr](validation_touched_examples_20261015.stderr). |
| Formatting | `rustfmt --check examples/nonisomorphic_extension_isogeny.rs examples/different_j_full_image_summary.rs examples/different_j_support_census.rs`: exit 0. | Checked on the final source. |
| Native measurement/replay | `zsh research/isogeny_torsion_fibers_20261012/different_j_full_image_20261015/run_v2.zsh`: exit 0; 18/18 scalar rows, eight base witnesses and 178 relation witnesses checked by a separate Rust process. | [run](RUN_V2.jsonl), [replay](REPLAY_V2.jsonl), 54 individual files in [raw-v2](raw-v2/). |
| Lean finite support | `lake build FullImageSupport` from [lean](lean/): exit 0; six `#print axioms` declarations report no axioms. | [stdout](LEAN_BUILD_V1.stdout), [stderr](LEAN_BUILD_V1.stderr), [proof map](PROOF_MAP.md). |

`lean4checker --fresh` is not installed here. The proof map records the
unformalized curve-level bridge. No Python research path was introduced
or run. The touched examples have no Python test suite.
