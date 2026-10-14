#!/usr/bin/env bash
# Build/run native verification; retain failures and continue independent checks.
set -u
export PATH=/workspace/toolchains/cargo/bin:$PATH
export RUSTUP_HOME=/workspace/toolchains/rustup
export CARGO_HOME=/workspace/toolchains/cargo
task_records=research/f4_f5_dense_packing_20261012/records/verification
mkdir -p "$task_records"

run_check() {
    task_name=$1
    shift
    printf '%q ' "$@" > "$task_records/$task_name-command.txt"
    printf '\n' >> "$task_records/$task_name-command.txt"
    date -u '+%Y-%m-%dT%H:%M:%SZ' > "$task_records/$task_name-start.txt"
    "$@" > "$task_records/$task_name.log" 2>&1
    task_status=$?
    printf '%s\n' "$task_status" > "$task_records/$task_name-exit-status.txt"
    date -u '+%Y-%m-%dT%H:%M:%SZ' > "$task_records/$task_name-end.txt"
    printf 'Completed %s: exit %s\n' "$task_name" "$task_status"
}

run_check focused-f4 timeout -k 5 1200 cargo test --locked --release --lib \
    cryptanalysis::koblitz_groebner::tests -j2 -- --test-threads=3
run_check focused-f5 timeout -k 5 1200 cargo test --locked --release --lib \
    cryptanalysis::matrix_f5_f2::tests -j2 -- --test-threads=3
run_check focused-pq timeout -k 5 1200 cargo test --locked --release --lib \
    cryptanalysis::pq_groebner_f2::tests -j2 -- --test-threads=3
run_check library-catalog cargo test --locked --release --lib -j2 -- --list
run_check library timeout -k 5 1200 cargo test --locked --release --lib -j2 -- --test-threads=3
run_check examples timeout -k 5 1200 cargo test --locked --release \
    --example f4_f5_hash_compare --example f4_f5_row_build_bench -j2
