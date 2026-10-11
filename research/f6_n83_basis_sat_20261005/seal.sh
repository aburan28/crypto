#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
out="$root/research/f6_n83_basis_sat_20261005"
cd "$root"
printf 'file\traw_bytes\traw_sha256\tgzip_bytes\tgzip_sha256\n' > "$out/XCNF_HASHES.tsv"
for raw in "$out"/*.xcnf; do
    packed="$raw.gz"
    bytes=$(stat -f '%z' "$raw")
    raw_hash=$(shasum -a 256 "$raw" | awk '{print $1}')
    gzip -n -9 -c "$raw" > "$packed"
    gzip -dc "$packed" | cmp - "$raw"
    packed_bytes=$(stat -f '%z' "$packed")
    packed_hash=$(shasum -a 256 "$packed" | awk '{print $1}')
    printf '%s\t%s\t%s\t%s\t%s\n' "$(basename "$raw")" "$bytes" "$raw_hash" "$packed_bytes" "$packed_hash" >> "$out/XCNF_HASHES.tsv"
    rm "$raw"
done
printf 'file\traw_bytes\traw_sha256\tgzip_bytes\tgzip_sha256\n' > "$out/SAT_MODEL_HASHES.tsv"
raw="$out/basis_planted_all_0.solver.stdout"
packed="$raw.gz"
bytes=$(stat -f '%z' "$raw")
raw_hash=$(shasum -a 256 "$raw" | awk '{print $1}')
gzip -n -9 -c "$raw" > "$packed"
gzip -dc "$packed" | cmp - "$raw"
packed_bytes=$(stat -f '%z' "$packed")
packed_hash=$(shasum -a 256 "$packed" | awk '{print $1}')
printf '%s\t%s\t%s\t%s\t%s\n' "$(basename "$raw")" "$bytes" "$raw_hash" "$packed_bytes" "$packed_hash" >> "$out/SAT_MODEL_HASHES.tsv"
rm "$raw"
find research/f6_n83_basis_sat_20261005 -maxdepth 1 -type f ! -name SHA256SUMS ! -name seal.log -print | LC_ALL=C sort | xargs shasum -a 256 > "$out/SHA256SUMS"
shasum -a 256 -c "$out/SHA256SUMS"
