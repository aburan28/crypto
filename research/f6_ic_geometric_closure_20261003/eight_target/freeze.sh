#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
out=research/f6_ic_geometric_closure_20261003/eight_target
base=research/f6_ic_geometric_closure_20261003/three_way
binary=${IC_WORKER_BINARY:-/private/tmp/f6-ic-target/release/examples/ic_tournament_worker}
digest() { shasum -a 256 "$1" | awk '{print $1}'; }
test ! -e "$out/freeze.tsv"
test "$(digest "$binary")" = "$(sed -n 's/^binary_sha256=//p' "$base/freeze.txt")"
for arm in inherited_f4 f5 f6_ic; do
  test "$(jq -r .record.implementation.flags.binary_sha256 "$base/$arm-candidate.json")" = "$(digest "$binary")"
done
printf 'index\tseed\tworkload_id\tfixture_sha256\tf4_input_sha256\tf5_input_sha256\tf6_input_sha256\n' > "$out/freeze.tsv"
for index in 1 2 3 4 5 6 7 8; do
  seed=$((20261004030 + index))
  dir="$out/T$index"
  mkdir "$dir"
  jq -cnS --argjson seed "$seed" \
    '{algorithm_seed:20261004039,config:{},curve_a:1,degree:17,
      factor_base:{dimension:6,kind:"standard_subspace"},mode:"fixture",
      target_seeds:[$seed]}' > "$dir/fixture-input.json"
  "$binary" < "$dir/fixture-input.json" > "$dir/fixture.json"
  test "$(jq -r '.fixture.target_scalar_constructed' "$dir/fixture.json")" = false
  test "$(jq -r '.fixture.targets|length' "$dir/fixture.json")" = 1
  jq -cS --slurpfile fixture "$dir/fixture.json" --argjson seed "$seed" \
    '.algorithm_seed=20261004039 | .target_seeds=[$seed] |
      .targets=[$fixture[0].fixture.targets[0]|map(tonumber)] |
      .resource_envelope.total_wall_limit_seconds=180 |
      .question="prepared-one-target-f4-f5-f6-eight-target-variability-panel"' \
    "$base/workload-record.json" > "$dir/workload-record.json"
  workload_sha=$(jq -jcS . "$dir/workload-record.json" | shasum -a 256 | awk '{print $1}')
  jq -cnS --slurpfile record "$dir/workload-record.json" --arg sha "$workload_sha" \
    '{record:$record[0],record_sha256:$sha,workload_id:$sha[0:12]}' \
    > "$dir/workload.json"
  for arm in inherited_f4 f5 f6_ic; do
    jq -cS --slurpfile fixture "$dir/fixture.json" \
      '.public_targets=$fixture[0].fixture.targets | .algorithm_seed=20261004039' \
      "$base/input-$arm.json" > "$dir/input-$arm.json"
  done
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$index" "$seed" \
    "$(jq -r .workload_id "$dir/workload.json")" "$(digest "$dir/fixture.json")" \
    "$(digest "$dir/input-inherited_f4.json")" "$(digest "$dir/input-f5.json")" \
    "$(digest "$dir/input-f6_ic.json")" >> "$out/freeze.tsv"
done
printf 'source_commit=%s\nbinary_sha256=%s\nfrozen_candidate_f4=%s\nfrozen_candidate_f5=%s\nfrozen_candidate_f6=%s\n' \
  "$(sed -n 's/^source_commit=//p' "$base/freeze.txt")" "$(digest "$binary")" \
  "$(digest "$base/inherited_f4-candidate.json")" \
  "$(digest "$base/f5-candidate.json")" \
  "$(digest "$base/f6_ic-candidate.json")" > "$out/freeze.txt"
