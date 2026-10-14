# Transfer assessment validator

This dependency-light Rust tool checks the fail-closed, machine-readable core
of [`docs/transfer-assessment/METHODOLOGY.md`](../../docs/transfer-assessment/METHODOLOGY.md).
It also implements the solver and net-cost gates in
[`MEASUREMENT.md`](../../docs/transfer-assessment/MEASUREMENT.md). It validates
the four-state vocabulary, evidence bindings, cumulative verdict gates, named
solver applicability, three-arm cost records, nine-phase totals, net-ratio
arithmetic, measured-evidence promotion gates, and diagnostic-only status
boundary.

The v3 record deliberately separates:

- `solver_assessment`: the named solver, exact image-subgroup access, search
  or relation compatibility, and individual-log/recovery interface; and
- `advantage_measurement`: source reference, destination generic control, and
  complete destination candidate, in one unit and one cold one-target scope.

Start a new assessment from
[`assessment-template.json`](../../docs/transfer-assessment/assessment-template.json)
and replace every placeholder before citing it as evidence.

`modeled_advantage=supported` is rejected for stage-only or amortized records,
unpriced work, missing arms, incomplete phase sums, a mismatched ratio, or a
ratio no greater than one. `measured_advantage=supported` additionally needs
measured arms, at least five verified pairs, four sizes, a replay receipt, and
a 95% net-speedup interval whose lower endpoint exceeds one. These are
consistency gates; the validator does not authenticate a cited receipt.

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
