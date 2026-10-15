# Exact native workflow

Measurement source commit c716e462934128c75e99dd458748e08623a2708c,
Rust1.99.0, locked dependencies, two build jobs:

    CRYPTO_BUILD_GIT_COMMIT=c716e462934128c75e99dd458748e08623a2708c CRYPTO_BUILD_GIT_DIRTY=false cargo build --release --locked --bin ecbench --bin hyperelliptic-cover --example ecbench_open_questions

Initial protocols/specs were committed before measurement. Source baseline,
release baseline and measured binary digests are recorded separately.

    ecbench_open_questions freeze NEW_STUDY_DIR
    ecbench plan --spec SPEC --json
    ecbench_open_questions_run MEASURED_ECBENCH STUDY_DIR
    ecbench_open_questions analyze STUDY_DIR
    ecbench_open_questions controls-extended STUDY_DIR
    ecbench_open_questions_replay_failures MEASURED_ECBENCH STUDY_DIR

The dispatcher preserves all exact run/audit/table/resources/compare arguments
and exit codes in execution-manifest.json. Sessions are never reused. The
supplemental replay uses the exact measured binary, not a silently rebuilt one.
The first failed supplement and corrected second supplement both remain.

Publication uses native-publish/Cargo.toml and the same Rust publisher source.
The standalone crate has a tracked lockfile and matching serde_json features.

    cargo run --release --locked --manifest-path research/ecbench_open_questions_20261014/native-publish/Cargo.toml
    cargo run --release --locked --manifest-path research/ecbench_open_questions_20261014/native-publish/Cargo.toml -- --check

Checks: cargo test --release --locked --lib, --test ecbench, --bin
hyperelliptic-cover, and the touched native examples. The first whole-library
run had nine local socket permission failures; targeted net/cairn retries with
network access pass. Preserve those failures. Clippy uses the repository CI's
pinned Rust1.98; an additional Rust1.99 attempt exposed a pre-existing PMAC
bit-width lint and is preserved, not suppressed or called a passing check.
Rustfmt and the native metadata checker must pass.

Rendering is the only Python step: it reads render-data.json and draws/layouts
recorded values without deriving costs, statistics, success counts or verdicts.
Every PDF is inspected; SVG/PDF/150dpi PNG copies go to benchmarks.

Independent-host/L2/GPU submissions have not run. No background research job is
claimed without an actual job identifier and receipt.
