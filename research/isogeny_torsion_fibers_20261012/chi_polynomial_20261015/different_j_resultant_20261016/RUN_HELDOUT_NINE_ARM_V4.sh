#!/usr/bin/env bash
set -u

base=research/isogeny_torsion_fibers_20261012/chi_polynomial_20261015/different_j_resultant_20261016
output="$base/heldout_nine_arm_v4"
mkdir -p "$output"
public="$base/FOUR_TARGET_PUBLIC_V1.json"
chi_binary=target/release/examples/different_j_chi_resultant
inverse_binary=target/release/examples/different_j_inverse_velu_oracle
meter=target/release/examples/different_j_chi_heldout_meter
arms=(chi-resultant direct-pair source-s3 quotient-s3 inverse-velu full-image source-bsgs rho-negation rho-negation-aa)

for round in 0 1 2; do
    for offset in 0 1 2 3 4 5 6 7 8; do
        arm_index=$(((offset + round) % 9))
        arm="${arms[$arm_index]}"
        binary="$chi_binary"
        if [[ "$arm" == "inverse-velu" || "$arm" == "full-image" ]]; then
            binary="$inverse_binary"
        fi
        for position in 0 1 2 3; do
            id="r${round}_${arm}_t${position}"
            "$meter" "$binary" "$public" "$round" "$arm" "$position" \
                "$output/$id.json" "$output/$id.stderr" "$output/$id.meter.json" 300 \
                > /dev/null 2> "$output/$id.meter.stderr"
            printf '%s\n' "$?" > "$output/$id.exit"
        done
    done
done
