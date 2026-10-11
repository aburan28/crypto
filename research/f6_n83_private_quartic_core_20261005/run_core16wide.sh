#!/bin/sh
# Amendment 2: exact k=16 core with a bounded 6.5m-column reducer.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_private_quartic_core_20261005"
binary="$root/target/bench_bins/f6_n83_private_quartic_core16wide_probe"
export RAYON_NUM_THREADS=1
printf 'mode\tk\texit\trss_kill\n' > "$out/core16wide_status.tsv"
printf 'mode\telapsed_s\trss_kib\n' > "$out/core16wide_rss.tsv"
run_one() {
    mode=$1
    label=$2
    if [ "$mode" = planted ]; then
        gtimeout -k 5s 300s "$binary" planted 16 > "$out/$label.jsonl" 2> "$out/$label.stderr.txt" &
    else
        gtimeout -k 5s 300s "$binary" ordinary 0 16 core16wide > "$out/$label.jsonl" 2> "$out/$label.stderr.txt" &
    fi
    wrapper=$!
    elapsed=0
    rss_kill=0
    while kill -0 "$wrapper" 2>/dev/null; do
        child=$(pgrep -P "$wrapper" | head -n 1)
        if [ -n "$child" ]; then
            rss=$(ps -o rss= -p "$child" | tr -d ' ')
            if [ -n "$rss" ]; then
                printf '%s\t%s\t%s\n' "$mode" "$elapsed" "$rss" >> "$out/core16wide_rss.tsv"
                if [ "$rss" -gt 7340032 ]; then
                    rss_kill=1
                    kill -KILL "$child" "$wrapper" 2>/dev/null
                    break
                fi
            fi
        fi
        sleep 1
        elapsed=$((elapsed + 1))
    done
    wait "$wrapper"
    code=$?
    printf '%s\t16\t%s\t%s\n' "$mode" "$code" "$rss_kill" >> "$out/core16wide_status.tsv"
    printf '%s k=16 exit=%s rss_kill=%s\n' "$mode" "$code" "$rss_kill"
    return "$code"
}
run_one planted core16wide_planted || exit 1
run_one ordinary core16wide_ordinary
