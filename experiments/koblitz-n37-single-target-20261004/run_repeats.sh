#!/bin/sh
# Thin orchestration for preregistered pairs r3-r7. The cryptographic work,
# checks, and replay all run in the native Rust binaries.
set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 <release-example-binary-directory>" >&2
    exit 2
fi
bin_dir=$1
study=experiments/koblitz-n37-single-target-20261004
target=$study/target_points.jsonl
timeout_bin=${TIMEOUT_BIN:-gtimeout}

run_ic() {
    out=$1
    KIC_DUMP_BASE="$out/base.jsonl" KIC_DUMP_RANK="$out/rank.jsonl" \
        "$timeout_bin" 900 "$bin_dir/koblitz_orbit_dlp_fast_online" \
        construct:37:0:7 "$target" 3737001 "$out/ic_targets.jsonl" \
        > "$out/ic_summary.jsonl" 2> "$out/ic.stderr"
}

run_rho() {
    out=$1
    "$timeout_bin" 900 "$bin_dir/koblitz_rho_fixture" \
        37 0 signed_frobenius 1 strong 370041 hash:370413 \
        > "$out/rho.jsonl" 2> "$out/rho.stderr"
    jq -c .published_q "$out/rho.jsonl" | cmp - "$target"
}

for repeat in 3 4 5 6 7; do
    out=$study/pairs/r$repeat
    if [ -e "$out" ]; then
        echo "refusing to overwrite $out" >&2
        exit 1
    fi
    mkdir -p "$out"
    if [ $((repeat % 2)) -eq 1 ]; then
        run_ic "$out"
        run_rho "$out"
    else
        run_rho "$out"
        run_ic "$out"
    fi
    "$timeout_bin" 900 "$bin_dir/koblitz_one_target_replay" \
        "$out/base.jsonl" "$out/rank.jsonl" "$out/ic_targets.jsonl" \
        "$out/ic_summary.jsonl" "$out/rho.jsonl" "$target" \
        > "$out/replay.jsonl" 2> "$out/replay.stderr"
done
