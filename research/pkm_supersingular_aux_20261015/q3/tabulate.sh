#!/bin/bash
# Lay out the sealed Q3 records and ecbench's saved comparisons as TSV (field
# extraction, sums and means only; the ratios and intervals are ecbench's).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
{
  printf 'curve\tarm\truns\tverified\tmean_S\tmean_total_gae\tmean_aux_share\taux_search_gae\taux_build_gae\trelations\ttargets_tried\tyield\tabscissae\tfb_id\n'
  for n in s20 s22 s24; do
    jq -s -r --arg n "$n" '
      map(select(.warmup == false)) | group_by(.arm)[] |
      . as $g |
      ($g | map(([.phases[] | select(.name == "factor_base") | ((.native.aux_search_gae // 0) + (.native.aux_tower_build_gae // 0))][0] // 0) / .cost.total_gae)) as $sh |
      [$n, $g[0].arm, ($g | length), ($g | map(select(.outcome.status == "verified")) | length),
       ($g | map(.cost.s) | add / length),
       ($g | map(.cost.total_gae) | add / length),
       ($sh | add / length),
       ([$g[0].phases[]? | select(.name == "factor_base") | .native.aux_search_gae // empty][0] // ""),
       ([$g[0].phases[]? | select(.name == "factor_base") | .native.aux_tower_build_gae // empty][0] // ""),
       ($g | map(.counters.relations_found // 0) | add),
       ($g | map(.counters.targets_tried // 0) | add),
       (if ($g | map(.counters.targets_tried // 0) | add) > 0 then (($g | map(.counters.relations_found // 0) | add) / ($g | map(.counters.targets_tried // 0) | add)) else "" end),
       ($g[0].factor_base.abscissae // ""), ($g[0].factor_base.fb_id // "")] | @tsv' sessions/$n/records.jsonl
  done
} > tables/summary-arms.tsv
{
  printf 'curve\tB\tA\tratio\tlo\thi\tpairs\n'
  for f in tables/s*.compare.*.txt; do
    n=$(basename "$f" | cut -d. -f1)
    grep -m1 '^ops' "$f" | sed -E "s/^ops \(ok\): ([^ ]+) \/ ([^ ]+) = ([0-9.]+) \[([0-9.]+), ([0-9.]+)\].* over ([0-9]+) pairs.*/$n\t\1\t\2\t\3\t\4\t\5\t\6/"
  done
} > tables/summary-ratios.tsv
