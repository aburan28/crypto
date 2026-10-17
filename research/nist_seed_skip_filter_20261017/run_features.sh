#!/bin/sh
# One GP process per run (6 in parallel), then the analysis. Thin orchestration only.
#   research/nist_seed_skip_filter_20261017/run_features.sh   (from the repository root)
set -eu
d=research/nist_seed_skip_filter_20261017
for r in 1 2 3 4 5 6; do
  printf 'RSEL=[%s];\n' "$r" | cat - "$d/features.gp" | gp -q 2>/dev/null > "$d/out/features.$r.part" &
done
wait
cat "$d"/out/features.[1-6].part > "$d/out/features.jsonl" && rm "$d"/out/features.[1-6].part
gp -q "$d/analysis.gp" | tee "$d/out/analysis.txt"
