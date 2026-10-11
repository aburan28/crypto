#!/usr/bin/env bash
set -euo pipefail

root_dir=$(git rev-parse --show-toplevel)
study_dir="$root_dir/research/ghs_curve_comparison_20261009"
out_dir="$study_dir/structural-non-koblitz"
mkdir -p "$out_dir"

jq -r '.workloads.curves[] | select(.b != 1) | [.n, .modulus, .a, .b, .r, .gx, .gy] | @tsv' "$study_dir/non_koblitz.spec.json" |
while IFS=$'\t' read -r degree modulus a b subgroup_order gx gy; do
    printf -v a_hex '0x%x' "$a"
    printf -v b_hex '0x%x' "$b"
    slug=$(jq -r --argjson n "$degree" --arg a_hex "$a_hex" --arg b_hex "$b_hex" \
        '.curves[] | select(.family == "binary" and .params.m == $n and .params.a == $a_hex and .params.b == $b_hex) | .slug' \
        "$root_dir/docs/curves/registry.json")
    if [[ -z "$slug" ]]; then
        printf 'unregistered curve degree=%s a=%s b=%s\n' "$degree" "$a" "$b" >&2
        exit 1
    fi
    "$root_dir/target/debug/ghs_screen" \
        --degree "$degree" --modulus "$modulus" --a "$a" --b "$b" \
        --genus-bound 4 --output "$out_dir/$slug.screen.json"
    jq -r '.towers[].relative_degree' "$out_dir/$slug.screen.json" |
    while read -r relative_degree; do
        "$root_dir/target/debug/ghs_transport" \
            --degree "$degree" --modulus "$modulus" --a "$a" --b "$b" \
            --relative-degree "$relative_degree" \
            --p-x "$gx" --p-y "$gy" --q-x "$gx" --q-y "$gy" \
            --order "$subgroup_order" \
            > "$out_dir/$slug.n$relative_degree.transport.json"
    done
done

printf 'source_commit=%s\n' "$(git rev-parse HEAD)" > "$out_dir/manifest.txt"
printf 'spec_sha256=' >> "$out_dir/manifest.txt"
shasum -a 256 "$study_dir/non_koblitz.spec.json" | cut -d' ' -f1 >> "$out_dir/manifest.txt"
