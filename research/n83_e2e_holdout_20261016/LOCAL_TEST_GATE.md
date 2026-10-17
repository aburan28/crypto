# Local release library gate: non-green receipt

The N83 holdout changes research Python, workflow, and result files only.
The branch and the then-current `origin/main` had identical files under
`src/` and `examples/` at this check. In particular,
`src/cryptanalysis/jv_cover.rs` had SHA-256
`a6d69fc93e10329b9608b0af6898fbd67867ace64d149e716d6830475424ca42`
on both refs. The focused N83 IC and rho example tests passed in
[screen run 38940108434](https://github.com/aburan28/crypto/actions/runs/38940108434);
the repaired holdout repeats them on its frozen cloud commit.

The repository-wide command
`cargo test --release --lib --locked --offline` compiled the release test
binary in 19m 11s and began 4,279 tests on the shared local macOS host.
It reported these four `jv_cover` failures before the agent interrupted the
already non-green run with exit 130:

- `the_solver_loses_no_rational_solution_at_small_p`
- `stopped_f4_agrees_with_the_full_solver_and_the_oracle`
- `replayed_f4_agrees_with_the_full_solver_and_the_oracle`
- `nagao_finds_planted_six_sums_and_agrees_with_the_oracle`

The first case passed alone in 10.89 s with
`cargo test --release --lib cryptanalysis::jv_cover::tests::the_solver_loses_no_rational_solution_at_small_p -- --exact --nocapture --test-threads=1`.
The second passed alone in 108.49 s with the same command shape and its
exact test name. The third failed alone after 124.63 s at
`src/cryptanalysis/jv_cover.rs:3473` with `n = 27` against its `n >= 30`
assertion; it was invoked with
`cargo test --release --lib cryptanalysis::jv_cover::tests::replayed_f4_agrees_with_the_full_solver_and_the_oracle -- --exact --nocapture --test-threads=1`.
The fourth was not rerun alone. The full suite and the third isolated test
are **non-green**; the first two isolated passes do not change that status.

The `n` assertion counts completed nondegenerate solver cases after skipped
timed-out attempts. Local host contention is a possible contributor, but
the recorded failure does not establish its cause. The N83 producer checks,
frozen fixture replay, independent grid replay, and the cloud holdout gate
remain separate from this library-suite status. A future full-suite green
result requires a completed run on a suitable host, or a justified fix to
the unchanged `jv_cover` code and tests; neither is claimed here.
