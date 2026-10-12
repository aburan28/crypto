#!/bin/sh
set -u
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_karatsuba_sat_20261005"
probe="$root/target/bench_bins/f6_n83_karatsuba_sat_probe"
solver=/opt/homebrew/opt/cryptominisat/bin/cryptominisat5
export RAYON_NUM_THREADS=1
printf 'mode\toffset\tpin\tlimit_s\tsolver_exit\trss_kill\tverify_exit\n' > "$out/status.tsv"
printf 'mode\toffset\tpin\telapsed_s\trss_kib\n' > "$out/rss.tsv"
run_one() {
    mode=$1
    offset=$2
    pin=$3
    limit=$4
    label="karatsuba_${mode}_${pin}_${offset}"
    xcnf="$out/$label.xcnf"
    if ! gtimeout -k 5s 300s "$probe" emit_xor_karatsuba "$mode" "$offset" "$xcnf" "$pin" > "$out/$label.emit.jsonl" 2> "$out/$label.emit.stderr.txt"; then
        printf '%s\t%s\t%s\t%s\temit_failed\t0\tnot_run\n' "$mode" "$offset" "$pin" "$limit" >> "$out/status.tsv"
        return 1
    fi
    gtimeout -k 5s "${limit}s" "$solver" --threads=1 --verb=0 "$xcnf" > "$out/$label.solver.stdout" 2> "$out/$label.solver.stderr.txt" &
    wrapper=$!
    elapsed=0
    rss_kill=0
    while kill -0 "$wrapper" 2>/dev/null; do
        child=$(pgrep -P "$wrapper" | head -n 1)
        if [ -n "$child" ]; then
            rss=$(ps -o rss= -p "$child" | tr -d ' ')
            if [ -n "$rss" ]; then
                printf '%s\t%s\t%s\t%s\t%s\n' "$mode" "$offset" "$pin" "$elapsed" "$rss" >> "$out/rss.tsv"
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
    solver_exit=$?
    if "$probe" verify "$mode" "$offset" "$out/$label.solver.stdout" > "$out/$label.verify.jsonl" 2> "$out/$label.verify.stderr.txt"; then
        verify_exit=0
    else
        verify_exit=$?
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$mode" "$offset" "$pin" "$limit" "$solver_exit" "$rss_kill" "$verify_exit" >> "$out/status.tsv"
    printf '%s offset=%s pin=%s solver_exit=%s rss_kill=%s verify_exit=%s\n' "$mode" "$offset" "$pin" "$solver_exit" "$rss_kill" "$verify_exit"
    return 0
}
run_one planted 0 all 30 || exit 1
if ! jq -e '.status == "verified_planted"' "$out/karatsuba_planted_all_0.verify.jsonl" >/dev/null; then
    exit 1
fi
run_one planted 0 source 30 || exit 1
for offset in 0 1 2 3; do
    run_one ordinary "$offset" none 120 || exit 1
done
