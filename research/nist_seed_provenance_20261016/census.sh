#!/bin/sh
# Rebuild the census of one run: enumerate candidates (Rust), point-count
# them in 4 parallel GP processes, classify (GP). Thin orchestration only.
#   research/nist_seed_provenance_20261016/census.sh <srand48 value> <blocks> <p in decimal>
set -eu
run=$1; blocks=$2; p=$3
here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d)
cargo run --quiet --release --example nist_seed_candidates -- --run "$run" --blocks "$blocks" > "$work/cand.txt"
split -n l/4 -d "$work/cand.txt" "$work/chunk."
for c in "$work"/chunk.*; do
  printf 'P=%s;\nCAND="%s";\n' "$p" "$c" | cat - "$here/sea_census.gp" | gp -q 2>/dev/null > "$c.out" &
done
wait
cat "$work"/chunk.*.out > "$work/orders.txt"
printf 'P=%s;\nORDERS="%s";\nRUN="srand48(%s)";\n' "$p" "$work/orders.txt" "$run" | cat - "$here/classify.gp" | gp -q
echo "orders: $work/orders.txt"
