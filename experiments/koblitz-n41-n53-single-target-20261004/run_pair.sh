#!/bin/sh
# Thin, no-overwrite orchestration. All cryptographic work and replay use Rust.
set -eu

if [ "$#" -ne 3 ]; then
    echo "usage: $0 <n:41|53> <repeat:1..6> <release-example-binary-directory>" >&2
    exit 2
fi
n=$1
repeat=$2
bin_dir=$3
case "$n" in
    41) columns=85; rank_seed=410041; rho_seed=410041; hash_seed=41261004 ;;
    53) columns=220; rank_seed=530053; rho_seed=530053; hash_seed=53261004 ;;
    *) echo "unsupported frozen cell: $n" >&2; exit 2 ;;
esac
case "$repeat" in
    1|2|3|4|5|6) ;;
    *) echo "repeat must be 1..6" >&2; exit 2 ;;
esac

study=experiments/koblitz-n41-n53-single-target-20261004
target=$study/targets/n$n.jsonl
out=$study/runs/n$n/r$repeat
if [ -e "$out" ]; then
    echo "refusing to overwrite $out" >&2
    exit 1
fi
if [ "$repeat" -eq 1 ]; then
    if [ -e "$target" ]; then
        echo "refusing to replace frozen target $target" >&2
        exit 1
    fi
else
    if [ ! -s "$target" ]; then
        echo "frozen target is missing: $target" >&2
        exit 1
    fi
fi
mkdir -p "$out" "$(dirname "$target")"
timeout_bin=${TIMEOUT_BIN:-gtimeout}
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1 LC_ALL=C

run_rho() {
    /usr/bin/time -l "$timeout_bin" 900 "$bin_dir/koblitz_rho_fixture" \
        "$n" 0 signed_frobenius 1 strong "$rho_seed" "hash:$hash_seed" \
        > "$out/rho.jsonl" 2> "$out/rho.stderr"
    if [ "$repeat" -eq 1 ]; then
        jq -c .published_q "$out/rho.jsonl" > "$target"
    else
        jq -c .published_q "$out/rho.jsonl" | cmp - "$target"
    fi
}

run_ic() {
    KIC_DUMP_BASE="$out/base.jsonl" KIC_DUMP_RANK="$out/rank.jsonl" \
        /usr/bin/time -l "$timeout_bin" 900 "$bin_dir/koblitz_orbit_dlp_fast_online" \
        "construct:$n:0:$columns" "$target" "$rank_seed" "$out/ic_targets.jsonl" \
        > "$out/ic_summary.jsonl" 2> "$out/ic.stderr"
}

if [ $((repeat % 2)) -eq 1 ]; then
    run_rho
    run_ic
else
    run_ic
    run_rho
fi

"$timeout_bin" 900 "$bin_dir/koblitz_one_target_replay" \
    "$out/base.jsonl" "$out/rank.jsonl" "$out/ic_targets.jsonl" \
    "$out/ic_summary.jsonl" "$out/rho.jsonl" "$target" \
    "$hash_seed" "$rho_seed" \
    > "$out/replay.jsonl" 2> "$out/replay.stderr"
