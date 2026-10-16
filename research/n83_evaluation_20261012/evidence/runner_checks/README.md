# Native runner and export validation

Performed after all 22 timed children completed. These checks do not retime either solver or modify the frozen panel, protocol, primary statistics or original compiled source.

- `runner_delivery_tests.log`: 7 tests passed, 0 failed. Delivery source uses portable temporary fixture paths.
- `descriptive_delivery_tests.log`: 4 tests passed, 0 failed. Tests cover reciprocal bounds, absence of a one-pair confidence interval, missing/invalid values, and separation of failed attempt charges from successful-solve costs without zero imputation.
- `descriptive_delivery_build.log`: release build passed.
- `descriptive_delivery_tests_attempt1.log`: retained failed compile caused by a delivery-only test function declaration typo, corrected before the passing test run. The original compiled exporter and generated statistics were unaffected.
- `secondary_integrity_checks.json` and `secondary_row_checks.json`: 15 Boolean checks passed. All 30 referenced formal-row raw-file hashes match; counts, values, descriptive means/min/max and primary reciprocal presentations agree with retained sources. Floating-point arithmetic comparisons use relative tolerance `1e-12`.
- `smoke_timeout/`: original tiny process timeout check, preserved from before the campaign. It observed `TIMEOUT`, SIGTERM 15, 0.305179865 seconds parent wall, and 524288 bytes Linux child peak RSS. This is runner validation, excluded from cryptanalytic measurements.

The delivery package commands, from the repository root, were equivalent to:

```sh
CARGO_BUILD_JOBS=2 CARGO_TARGET_DIR=/workspace/n83_evaluation/runner_delivery_test_build cargo test --release --locked --offline --manifest-path research/n83_evaluation_20261012/runner/Cargo.toml
CARGO_BUILD_JOBS=2 CARGO_TARGET_DIR=/workspace/n83_evaluation/descriptive_delivery_test_build cargo test --release --locked --offline --manifest-path research/n83_evaluation_20261012/descriptive_summary/Cargo.toml
CARGO_BUILD_JOBS=2 CARGO_TARGET_DIR=/workspace/n83_evaluation/descriptive_delivery_test_build cargo build --release --locked --offline --manifest-path research/n83_evaluation_20261012/descriptive_summary/Cargo.toml
```

The cached session toolchain supplied `PATH`, `CARGO_HOME` and `RUSTUP_HOME`; package targets were separate from the full repository's release test target. No large executable was copied into this evidence directory. Exact original compiled runner/exporter source and package manifests/locks are retained under `../source/`; delivery formatting and test changes have separate source identities.
