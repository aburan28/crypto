#!/usr/bin/env bash
set -u

if [[ "$#" -ne 1 ]]; then
    printf 'usage: %s OUTPUT_DIRECTORY\n' "$0" >&2
    exit 2
fi
output="$1"
mkdir -p "$output"
base=research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/inverse_velu/trend_four_cell_v1
public="$base/P2473_PUBLIC.json"
checker="$base/P2473_CHECKER.json"
solver=target/release/examples/different_j_inverse_velu_oracle
meter=target/release/examples/different_j_order_process_meter
verify=target/release/examples/different_j_inverse_velu_check

timeout -k 5s 300s prlimit --as=536870912 -- \
    "$meter" "$solver" "$public" 0 full-image 1 \
    "$output/run.json" "$output/run.stderr" "$output/meter.json" \
    > "$output/meter.stdout" 2> "$output/meter.stderr"
run_status=$?
printf '%s\n' "$run_status" > "$output/run.exit"
if [[ "$run_status" -ne 0 ]]; then
    printf 'censored or failed full-image child: exit %s\n' "$run_status" > "$output/STATUS.txt"
    (cd "$output" && sha256sum * > artifacts.sha256)
    exit 0
fi

"$verify" check-single-full-image "$public" "$checker" \
    "$output/run.json" "$output/meter.json" \
    > "$output/check.json" 2> "$output/check.stderr"
check_status=$?
printf '%s\n' "$check_status" > "$output/check.exit"
if [[ "$check_status" -eq 0 ]]; then
    printf 'complete and independently replayed\n' > "$output/STATUS.txt"
else
    printf 'checker failure: exit %s\n' "$check_status" > "$output/STATUS.txt"
fi
(cd "$output" && sha256sum * > artifacts.sha256)
exit "$check_status"
