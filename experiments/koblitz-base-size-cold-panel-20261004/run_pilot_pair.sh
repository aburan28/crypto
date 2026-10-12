#!/bin/sh
# Frozen thin orchestration; all ECC, rho, and replay work uses native Rust.
set -eu

if [ "$#" -ne 4 ]; then
    echo "usage: $0 <n:41|53> <K> <repeat:1|2|3> <release-example-binary-directory>" >&2
    exit 2
fi
n=$1
columns=$2
repeat=$3
bin_dir=$4
case "$n:$columns" in
    41:20|41:32|41:48|41:64|41:85)
        rank_seed=410041; rho_seed=410041; hash_seed=41261107 ;;
    53:80|53:100|53:128|53:160|53:220)
        rank_seed=530053; rho_seed=530053; hash_seed=53261107 ;;
    *) echo "unsupported frozen pilot cell: $n:$columns" >&2; exit 2 ;;
esac
case "$repeat" in
    1|2|3) ;;
    *) echo "pilot repeat must be 1..3" >&2; exit 2 ;;
esac

study=experiments/koblitz-base-size-cold-panel-20261004
target=$study/targets/n${n}_pilot.jsonl
out=$study/runs/pilot/n${n}/k${columns}/r${repeat}
if [ ! -s "$target" ]; then
    echo "frozen target is missing: $target" >&2
    exit 1
fi
if [ -e "$out" ]; then
    echo "refusing to overwrite $out" >&2
    exit 1
fi
mkdir -p "$out"
timeout_bin=${TIMEOUT_BIN:-timeout}
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 LC_ALL=C

run_rho() {
    /usr/bin/time "$timeout_bin" 30 "$bin_dir/koblitz_rho_fixture" \
        "$n" 0 signed_frobenius 1 strong "$rho_seed" "hash:$hash_seed" \
        > "$out/rho.jsonl" 2> "$out/rho.stderr"
}

run_ic() {
    KIC_DUMP_BASE="$out/base.jsonl" KIC_DUMP_RANK="$out/rank.jsonl" \
        /usr/bin/time "$timeout_bin" 30 "$bin_dir/koblitz_orbit_dlp_fast_online" \
        "construct:$n:0:$columns" "$target" "$rank_seed" "$out/ic_targets.jsonl" \
        > "$out/ic_summary.jsonl" 2> "$out/ic.stderr"
}

if [ $((repeat % 2)) -eq 1 ]; then
    if run_rho; then rho_status=0; else rho_status=$?; fi
    if run_ic; then ic_status=0; else ic_status=$?; fi
else
    if run_ic; then ic_status=0; else ic_status=$?; fi
    if run_rho; then rho_status=0; else rho_status=$?; fi
fi

if [ "$rho_status" -eq 0 ]; then
    if ! jq -ce .published_q "$out/rho.jsonl" | cmp - "$target"; then
        rho_status=98
        echo "rho public point differs from frozen target" > "$out/rho.point_mismatch"
    fi
fi

replay_status=null
if [ "$rho_status" -eq 0 ] && [ "$ic_status" -eq 0 ]; then
    if "$timeout_bin" 30 "$bin_dir/koblitz_one_target_replay" \
        "$out/base.jsonl" "$out/rank.jsonl" "$out/ic_targets.jsonl" \
        "$out/ic_summary.jsonl" "$out/rho.jsonl" "$target" \
        "$hash_seed" "$rho_seed" \
        > "$out/replay.jsonl" 2> "$out/replay.stderr"; then
        replay_status=0
    else
        replay_status=$?
    fi
else
    echo "replay skipped because rho or IC did not complete" > "$out/replay.skipped"
fi

jq -nc --arg phase pilot --argjson n "$n" --argjson K "$columns" \
    --argjson repeat "$repeat" --argjson rho_status "$rho_status" \
    --argjson ic_status "$ic_status" --argjson replay_status "$replay_status" \
    '{phase:$phase,n:$n,K:$K,repeat:$repeat,rho_status:$rho_status,ic_status:$ic_status,replay_status:$replay_status}' \
    > "$out/status.jsonl"
