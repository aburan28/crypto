#!/usr/bin/env bash
set -euo pipefail

root_dir=$(git rev-parse --show-toplevel)
out_dir="$root_dir/research/ghs_curve_comparison_20261009/structural-standard"
mkdir -p "$out_dir"
modulus=0x20000000000000000000000000201

"$root_dir/target/debug/ghs_screen" --degree 113 --modulus "$modulus" \
    --a 0x3088250ca6e7c7fe649ce85820f7 \
    --b 0xe8bee4d3e2260744188be0e9c723 --genus-bound 4 \
    --output "$out_dir/sect113r1.screen.json"
"$root_dir/target/debug/ghs_transport" --degree 113 --modulus "$modulus" \
    --a 0x3088250ca6e7c7fe649ce85820f7 \
    --b 0xe8bee4d3e2260744188be0e9c723 --relative-degree 113 \
    --p-x 0x009d73616f35f4ab1407d73562c10f \
    --p-y 0x00a52830277958ee84d1315ed31886 \
    --q-x 0x009d73616f35f4ab1407d73562c10f \
    --q-y 0x00a52830277958ee84d1315ed31886 \
    --order 0x0100000000000000d9ccec8a39e56f \
    > "$out_dir/sect113r1.transport.json"

"$root_dir/target/debug/ghs_screen" --degree 113 --modulus "$modulus" \
    --a 0x689918dbec7e5a0dd6dfc0aa55c7 \
    --b 0x95e9a9ec9b297bd4bf36e059184f --genus-bound 4 \
    --output "$out_dir/sect113r2.screen.json"
"$root_dir/target/debug/ghs_transport" --degree 113 --modulus "$modulus" \
    --a 0x689918dbec7e5a0dd6dfc0aa55c7 \
    --b 0x95e9a9ec9b297bd4bf36e059184f --relative-degree 113 \
    --p-x 0x01a57a6a7b26ca5ef52fcdb8164797 \
    --p-y 0x00b3adc94ed1fe674c06e695baba1d \
    --q-x 0x01a57a6a7b26ca5ef52fcdb8164797 \
    --q-y 0x00b3adc94ed1fe674c06e695baba1d \
    --order 0x010000000000000108789b2496af93 \
    > "$out_dir/sect113r2.transport.json"

printf 'source_commit=%s\n' "$(git rev-parse HEAD)" > "$out_dir/manifest.txt"
printf 'curve_source_sha256=' >> "$out_dir/manifest.txt"
shasum -a 256 "$root_dir/src/binary_ecc/curve.rs" | cut -d' ' -f1 >> "$out_dir/manifest.txt"
