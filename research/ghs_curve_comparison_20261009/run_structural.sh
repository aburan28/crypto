#!/usr/bin/env bash
set -euo pipefail

root_dir=$(git rev-parse --show-toplevel)
out_dir="$root_dir/research/ghs_curve_comparison_20261009/structural"
mkdir -p "$out_dir"

jq -r '.curves[] | select(.family == "koblitz" and ((.params.n == 7 and .params.a == 0) or (.params.n == 9) or (.params.n == 13 and .params.a == 0) or (.params.n == 15 and .params.a == 0))) | [.slug, .params.n, .params.modulus, .params.a, .representations[0].curve.subgroup_order, .representations[0].curve.generator[0], .representations[0].curve.generator[1]] | @tsv' "$root_dir/docs/curves/registry.json" |
while IFS=$'\t' read -r slug degree modulus a subgroup_order gx gy; do
    "$root_dir/target/debug/ghs_screen" \
        --degree "$degree" --modulus "$modulus" --a "$a" --b 1 \
        --genus-bound 4 --output "$out_dir/$slug.screen.json"
    jq -r '.towers[].relative_degree' "$out_dir/$slug.screen.json" |
    while read -r relative_degree; do
        "$root_dir/target/debug/ghs_transport" \
            --degree "$degree" --modulus "$modulus" --a "$a" --b 1 \
            --relative-degree "$relative_degree" \
            --p-x "$gx" --p-y "$gy" --q-x "$gx" --q-y "$gy" \
            --order "$subgroup_order" \
            > "$out_dir/$slug.n$relative_degree.transport.json"
    done
done

printf 'source_commit=%s\n' "$(git rev-parse HEAD)" > "$out_dir/manifest.txt"
printf 'registry_sha256=' >> "$out_dir/manifest.txt"
shasum -a 256 "$root_dir/docs/curves/registry.json" | cut -d' ' -f1 >> "$out_dir/manifest.txt"
printf 'screen_sha256=' >> "$out_dir/manifest.txt"
shasum -a 256 "$root_dir/target/debug/ghs_screen" | cut -d' ' -f1 >> "$out_dir/manifest.txt"
printf 'transport_sha256=' >> "$out_dir/manifest.txt"
shasum -a 256 "$root_dir/target/debug/ghs_transport" | cut -d' ' -f1 >> "$out_dir/manifest.txt"
