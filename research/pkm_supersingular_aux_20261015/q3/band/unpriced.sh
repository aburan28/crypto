#!/bin/bash
# Every arm's unpriced counters with their counts (AGENTS.md §5a), from the
# sealed records, measured runs only: the same lookup as
# ecbench::band::unpriced_counts (top-level counter, else the phases' native
# `X_uncharged` or `X`). Sums and means only.
set -eu
here=$(cd "$(dirname "$0")/.." && pwd)
cd "$here"
printf 'curve\tarm\tcounter\truns\truns_with_count\ttotal\tmean_per_run\n'
for n in s20 s22 s24; do
  jq -s -r --arg n "$n" '
    map(select(.warmup == false)) | group_by(.arm)[] | . as $g |
    ([$g[] | .cost.unpriced[]] | unique)[] as $c |
    ($c | sub("_uncharged$"; "")) as $b |
    [$g[] | . as $r |
      (if ($r.counters[$c] != null) then $r.counters[$c]
       else ([$r.phases[] | (.native[$c] // .native[$b])] | map(select(. != null))) as $v
            | (if ($v | length) > 0 then ($v | add) else null end) end)] as $vals |
    [$n, $g[0].arm, $c, ($g | length), ($vals | map(select(. != null)) | length),
     ($vals | map(select(. != null)) | add // 0),
     (($vals | map(select(. != null)) | add // 0) / ($g | length))] | @tsv' sessions/$n/records.jsonl
done
