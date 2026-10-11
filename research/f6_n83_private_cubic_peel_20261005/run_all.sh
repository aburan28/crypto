#!/bin/sh
# Amendment 1: every exact cubic monomial, with a live RSS gate.
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_private_cubic_peel_20261005"
binary="$root/target/bench_bins/f6_n83_private_cubic_all_probe"
export RAYON_NUM_THREADS=1
printf 'mode\tk\texit\trss_kill\n' > "$out/all_status.tsv"
printf 'mode\telapsed_s\trss_kib\n' > "$out/all_rss.tsv"
run_one() {
    mode=$1
    if [ "$mode" = planted ]; then
        gtimeout -k 5s 300s "$binary" planted 90 > "$out/all_planted.jsonl" 2> "$out/all_planted.stderr.txt" &
    else
        gtimeout -k 5s 300s "$binary" ordinary 0 90 cubic90all > "$out/all_ordinary.jsonl" 2> "$out/all_ordinary.stderr.txt" &
    fi
    wrapper=$!
    elapsed=0
    rss_kill=0
    while kill -0 "$wrapper" 2>/dev/null; do
        child=$(pgrep -P "$wrapper" | head -n 1)
        if [ -n "$child" ]; then
            rss=$(ps -o rss= -p "$child" | tr -d ' ')
            if [ -n "$rss" ]; then
                printf '%s\t%s\t%s\n' "$mode" "$elapsed" "$rss" >> "$out/all_rss.tsv"
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
    printf '%s\t90\t%s\t%s\n' "$mode" "$code" "$rss_kill" >> "$out/all_status.tsv"
    printf '%s k=90 exit=%s rss_kill=%s\n' "$mode" "$code" "$rss_kill"
    return "$code"
}
run_one planted || exit 1
run_one ordinary
