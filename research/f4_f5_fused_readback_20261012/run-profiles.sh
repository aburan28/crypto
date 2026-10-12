#!/usr/bin/env bash
# Thin orchestration only: the Rust driver executes and verifies every call.
set -euo pipefail
task_record_root=research/f4_f5_fused_readback_20261012/records
task_profiler=/workspace/scratch/f4-f5-hash/valgrind/usr/bin/valgrind
export VALGRIND_LIB=/workspace/scratch/f4-f5-hash/valgrind/usr/libexec/valgrind

profile() {
    task_name=$1
    task_variant=$2
    task_repeats=$3
    task_destination=$task_record_root/$task_name
    # Never overwrite a completed or failed measurement.
    if [ -e "$task_destination" ]; then
        printf 'Refusing existing run directory: %s\n' "$task_destination" >&2
        return 1
    fi
    mkdir -p "$task_destination"
    date -u '+%Y-%m-%dT%H:%M:%SZ' > "$task_destination/start-utc.txt"
    printf '%s\n' "VALGRIND_LIB=$VALGRIND_LIB timeout -k 5 180 $task_profiler --tool=callgrind --instr-atstart=no --separate-threads=no --callgrind-out-file=$task_destination/callgrind.out /workspace/scratch/f4-f5-fused-readback/$task_variant/bench profile $task_repeats" > "$task_destination/command.txt"
    if timeout -k 5 180 "$task_profiler" --tool=callgrind --instr-atstart=no --separate-threads=no \
        --callgrind-out-file="$task_destination/callgrind.out" \
        "/workspace/scratch/f4-f5-fused-readback/$task_variant/bench" profile "$task_repeats" \
        > "$task_destination/records.jsonl" 2> "$task_destination/profiler.log"; then
        printf '0\n' > "$task_destination/exit-status.txt"
    else
        task_status=$?
        printf '%s\n' "$task_status" > "$task_destination/exit-status.txt"
        return "$task_status"
    fi
    date -u '+%Y-%m-%dT%H:%M:%SZ' > "$task_destination/end-utc.txt"
    printf 'Completed %s\n' "$task_name"
}

if [ "${1:-}" = "baseline" ]; then
    profile primary-a a 3
    profile aa-a a 3
    exit
fi

if [ "${1:-}" = "legacy" ]; then
    export KIC_F5_FUSED_BUILD=0
    profile legacy-b b 1
    exit
fi

# primary-a was captured before any production edits; aa-a repeats that binary.
profile primary-b b 3
for task_round in 0 1 2 3 4; do
    case $((task_round % 3)) in
        0) profile "round-$task_round-a" a 1; profile "round-$task_round-b" b 1; profile "round-$task_round-aprime" a 1 ;;
        1) profile "round-$task_round-b" b 1; profile "round-$task_round-aprime" a 1; profile "round-$task_round-a" a 1 ;;
        2) profile "round-$task_round-aprime" a 1; profile "round-$task_round-a" a 1; profile "round-$task_round-b" b 1 ;;
    esac
done
