# Integration validation commands

Worktree: /workspace/crypto. Integration source:
7b0654c768f5b6775918d759446f400a6aebc03d, merging exact main
a80b15bbc5c1ea4f2db2b294a9bcddc31f06fdd7. Subsequent edits correct only the
leads example's redundant borrow, decimal-reference formatting and validation
documentation. Measurement remains source0cbc with its original binary.

    export PATH=/workspace/toolchains/cargo/bin:$PATH
    export RUSTUP_HOME=/workspace/toolchains/rustup
    export CARGO_HOME=/workspace/toolchains/cargo
    export CARGO_BUILD_JOBS=2
    cargo +1.98.0 test --release --lib
    cargo +1.98.0 test --release --bin ecbench --test ecbench --test ecbench_large_prime --example ecbench_open_questions --example ecbench_open_questions_run --example ecbench_open_questions_replay_failures --example ecbench_open_questions_publish --example ecbench_leads --example ecbench_leads_findings
    cargo +1.98.0 test --release --example ecbench_leads
    cargo +1.98.0 build --release --bin ecbench
    cargo +1.98.0 clippy --release --bin ecbench --example ecbench_open_questions --example ecbench_open_questions_run --example ecbench_open_questions_replay_failures --example ecbench_open_questions_publish --example ecbench_leads --example ecbench_leads_findings -- -D warnings
    cargo +1.98.0 test --release --bin curve_cover_check --bin curve_standards
    cargo +1.98.0 run --release --bin curve_standards -- --check
    cargo +1.98.0 run --release --bin curve_cover_check -- --check

The first clippy invocation fails on the redundant borrow. The corrected retry,
clippy-integrated-v2.log, passes without a lint override. Both attempts remain.
The full library log is a single successful invocation, not a sum of retries.

For each of the six new and nine parent session directories:

    target/release/ecbench verify --dir SESSION --replay-all --exit-code --out NEW_RECEIPT

All receipts are new files under integration-audits/. No receipt or session is
reused or overwritten. The verifier digest, build-source commit and receipt
index sit in that directory. Verification with the integrated binary is not an
independent-host measurement.

    cargo +1.98.0 test --release --locked --manifest-path research/ecbench_open_questions_20261014/native-publish/Cargo.toml
    cargo +1.98.0 clippy --release --locked --manifest-path research/ecbench_open_questions_20261014/native-publish/Cargo.toml -- -D warnings
    research/ecbench_open_questions_20261014/native-publish/target/release/ecbench-ten-question-publish --historical-inputs
    research/ecbench_open_questions_20261014/native-publish/target/release/ecbench-ten-question-publish --check
    python3 scripts/update_curve_standards.py --check
    python3 scripts/check_curve_names.py
    python3 scripts/test_cover_catalog.py
    python3 scripts/test_check_ic_leaderboard_paths.py
    python3 scripts/check_ic_leaderboard_paths.py
    python3 -m unittest discover -s research/sat_factor_base_review_20260908/autolab -p test_boundary_autolab.py
    node tools/render_ic_dashboard.mjs --check

Python here checks existing metadata/legacy control-plane contracts or renders
already recorded values; it supplies no new experiment, arithmetic or research
analysis. The first legacy-suite attempt lacks a tracked sparse input; the
retry follows scripts/sparse-checkout.sh add --path for oracle_u128.py.

In detached /workspace/historical-metadata-check at exact a80:

    python3 scripts/build_ic_leaderboard.py --check
    python3 scripts/build_lab_browser.py --check
    cargo +1.98.0 fmt --all -- --check
    python3 research/notes/ecc2k130/leaf_m10_support_outcome_20260930/verify_archive.py

The first two pass. The final two retain the baseline formatting and archive
failures. Current-tree cargo fmt produces the identical normalized baseline
diff. A separate rustfmt1.98 check of every changed Rust file passes after the
reference-crate layout correction.

The file-URL browser attempt fails under the runtime's managed navigation
policy. The passing retry uses the same document with all assertions retained:

    IC_DASHBOARD_DOCUMENT_MODE=inline IC_DASHBOARD_CHROME=/tmp/ecbench-dashboard-chromium IC_DASHBOARD_SCREENSHOTS=research/ecbench_leads_20261015/presentation node tools/check_ic_dashboard_browser.mjs

The temporary executable invokes /usr/bin/chromium --no-sandbox with its
original arguments, inside the managed execution sandbox. No policy is edited.
The receipt identifies inline mode and the exact source-page hash. Screenshots,
both failed file-navigation attempts and both successful inline runs remain.

    python3 research/ecbench_leads_20261015/render_reports.py
    sha256sum -c research/ecbench_leads_20261015/pdf-artifacts.sha256

PDFs were rasterized with pdftoppm and inspected; benchmark PDF/SVG/PNG copies
were compared byte for byte. Frozen inputs were hashed from git show at exact
measured commit0cbc, not from the changed integration tree.
