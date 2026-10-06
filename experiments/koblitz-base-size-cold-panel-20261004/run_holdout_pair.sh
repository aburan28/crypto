#!/bin/sh
# Frozen thin orchestration; all ECC, rho, and replay work uses native Rust.
set -eu

if [ "$#" -ne 4 ]; then
    echo "usage: $0 <n:41|53> <K> <repeat:1..6> <release-example-binary-directory>" >&2
    exit 2
fi
n=$1
columns=$2
repeat=$3
bin_dir=$4
case "$n" in
    41) rank_seed=410041; rho_seed=410041; hash_seed=41261207 ;;
    53) rank_seed=530053; rho_seed=530053; hash_seed=53261207 ;;
    *) echo "unsupported frozen held-out size: $n" >&2; exit 2 ;;
esac
case "$repeat" in
    1|2|3|4|5|6) ;;
    *) echo "held-out repeat must be 1..6" >&2; exit 2 ;;
esac

study=experiments/koblitz-base-size-cold-panel-20261004
selection=$study/PILOT_ANALYSIS.json
baseline=$(jq -er --argjson n "$n" '.cells[] | select(.n==$n) | .baseline_K' "$selection")
selected=$(jq -er --argjson n "$n" '.cells[] | select(.n==$n) | .selected_smaller_K' "$selection")
if [ "$columns" -ne "$baseline" ] && [ "$columns" -ne "$selected" ]; then
    echo "K=$columns differs from published baseline=$baseline and selected=$selected" >&2
    exit 2
fi
target=$study/targets/n${n}_holdout.jsonl
out=$study/runs/holdout/n${n}/k${columns}/r${repeat}
if [ ! -s "$target" ]; then
    echo "frozen held-out target is missing: $target" >&2
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

jq -nc --arg phase holdout --argjson n "$n" --argjson K "$columns" \
    --argjson repeat "$repeat" --argjson rho_status "$rho_status" \
    --argjson ic_status "$ic_status" --argjson replay_status "$replay_status" \
    '{phase:$phase,n:$n,K:$K,repeat:$repeat,rho_status:$rho_status,ic_status:$ic_status,replay_status:$replay_status}' \
    > "$out/status.jsonl"
