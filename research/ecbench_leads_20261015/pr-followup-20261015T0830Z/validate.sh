#!/usr/bin/env bash
# Integration validation only; repository tests use temporary smoke fixtures.
set -uo pipefail
cd /workspace/crypto
export PATH=/workspace/toolchains/cargo/bin:$PATH
export RUSTUP_HOME=/workspace/toolchains/rustup
export CARGO_HOME=/workspace/toolchains/cargo
export CARGO_BUILD_JOBS=2
evidence=research/ecbench_leads_20261015/pr-followup-20261015T0830Z
run_check() {
    local name=$1
    shift
    printf '%q ' "$@" > "$evidence/$name.command"
    printf '\n' >> "$evidence/$name.command"
    "$@" > "$evidence/$name.log" 2>&1
    local result=$?
    printf '%s\n' "$result" > "$evidence/$name.exit"
    gzip -n "$evidence/$name.log"
    printf '%s: exit%s\n' "$name" "$result"
    return "$result"
}
run_check library cargo +1.98.0 test --release --lib || exit $?
run_check cli-examples cargo +1.98.0 test --release --bin ecbench --test ecbench --test ecbench_large_prime --example ecbench_open_questions --example ecbench_open_questions_run --example ecbench_open_questions_replay_failures --example ecbench_open_questions_publish --example ecbench_leads --example ecbench_leads_findings || exit $?
run_check all-targets cargo +1.98.0 check --release --all-targets || exit $?
run_check clippy cargo +1.98.0 clippy --release --bin ecbench --example ecbench_open_questions --example ecbench_open_questions_run --example ecbench_open_questions_replay_failures --example ecbench_open_questions_publish --example ecbench_leads --example ecbench_leads_findings -- -D warnings || exit $?
run_check native-freshness cargo +1.98.0 run --release --locked --manifest-path research/ecbench_open_questions_20261014/native-publish/Cargo.toml -- --check || exit $?
run_check site python3 scripts/site/test_build.py -v || exit $?
run_check workflow-parser python3 scripts/test_workflows_parse.py || exit $?
