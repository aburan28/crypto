# Transfer assessment validator

This dependency-light Rust tool checks the fail-closed, machine-readable core
of [`docs/transfer-assessment/METHODOLOGY.md`](../../docs/transfer-assessment/METHODOLOGY.md).
It validates the four-state vocabulary, evidence bindings, cumulative verdict
gates, complete-cost promotion gate, and diagnostic-only status boundary.

```sh
cargo test --release --locked --manifest-path tools/transfer-assessment/Cargo.toml
cargo run --release --locked --manifest-path tools/transfer-assessment/Cargo.toml -- \
  --output /tmp/summary.json research/transfer-methodology-validation-20261013/assessments/*.json
```

Exit status `0` means every assessment passed the contract. Status `2` means
invalid JSON, an unknown flag, or at least one semantic validation error;
status `1` is an I/O failure. Validation establishes internal consistency of
an assessment. It does not independently prove cited mathematics or promote a
transfer, speed, novelty, or security claim.
