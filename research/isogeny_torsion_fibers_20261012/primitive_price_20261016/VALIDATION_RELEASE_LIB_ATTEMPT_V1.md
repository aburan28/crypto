# Release-library validation attempt V1

Command from the research worktree, started before the result and CI
commits of this PR:

```sh
cargo test --release --lib --offline -- --test-threads=2
```

This shared Darwin arm64 run recorded 1,595 passing tests, 51 ignored
tests and one failure before it was interrupted with SIGINT (shell exit
130). The failing name in the retained log is
`cryptanalysis::jv_cover::tests::nagao_finds_planted_six_sums_and_agrees_with_the_oracle`.
The harness had not printed its failure detail or final suite summary
before interruption, so this receipt does not establish the failure
cause or full-suite outcome. The older PR-head broad attempt and its
focused `jv_cover` reruns are described in PR #1939; this new receipt
must not replace their exact statuses.

The [raw log](VALIDATION_RELEASE_LIB_ATTEMPT_V1.log) is 152,241 bytes,
SHA-256 `379d1546ee3fbfee3b8a44b22cd213c07294aa1a2480d859862c0987c9ba88d2`.
It was stopped after the failure was known, while two shared-host
research tests were running, to free capacity for the changed examples'
release test and Rust 1.98 Clippy gate. There is no green full-suite
claim from this attempt. The focused failing test must be rerun with
full output before assigning a cause.
