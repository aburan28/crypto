#!/usr/bin/env bash
set -u

base=research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/inverse_velu/trend_four_cell_v1
output="$base/darwin_v1"
mkdir -p "$output"
binary=target/release/examples/different_j_inverse_velu_oracle
meter=target/release/examples/different_j_order_process_meter
arms=(inverse-velu source-bsgs rho-negation rho-negation-aa)

for p in 2473 5521 9619 19213; do
    for round in 0 1 2; do
        for offset in 0 1 2 3; do
            arm_index=$(((offset + round) % 4))
            arm="${arms[$arm_index]}"
            for prefix in 1 64; do
                id="p${p}_r${round}_${arm}_n${prefix}"
                timeout 300 "$meter" "$binary" "$base/P${p}_PUBLIC.json" \
                    "$round" "$arm" "$prefix" "$output/$id.json" \
                    "$output/$id.stderr" "$output/$id.meter.json" \
                    > /dev/null 2> "$output/$id.meter.stderr"
                printf '%s\n' "$?" > "$output/$id.exit"
            done
        done
    done
done
