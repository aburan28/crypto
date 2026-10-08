#!/usr/bin/env bash
# Reduce runs/ to summary.json: per size and path, each arm's median, min
# and max of the collect and build phases, the whole-pipeline time and its
# units; the A/A spread; the A/B ratios; the identity check; and how many
# runs isolated_bench marked contended.  Computes nothing a report does not
# contain beyond medians, minima, maxima and their quotients.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
runs="$here/runs"
jq -n \
  --slurpfile aa <(for f in "$runs"/aa/*.price.json; do jq -c --arg f "$(basename "$f")" '{file: $f, r: .}' "$f"; done) \
  --slurpfile ab <(for f in "$runs"/ab/*.price.json; do jq -c --arg f "$(basename "$f")" '{file: $f, r: .}' "$f"; done) \
  --slurpfile iso <(cat "$runs"/aa/isolation.jsonl "$runs"/ab/isolation.jsonl) '
def med: sort | if length % 2 == 1 then .[length/2|floor] else (.[length/2-1] + .[length/2]) / 2 end;
def stats: {median: med, min: min, max: max, n: length};
def parse: .file | capture("^(?<size>icv1-[a-z0-9]+-[a-z0-9]+-[0-9a-f]{8})-(?<path>portable|avx512)-(?<arm>baseline|candidate)-r(?<round>[0-9]+)");
def metrics: {
  collect_ms: (.r.repetitions[0].phases_ns.collect / 1e6),
  build_ms: (.r.repetitions[0].phases_ns.build / 1e6),
  total_ms: (.r.repetitions[0].total_ns / 1e6),
  total_units: .r.repetitions[0].total_units,
  s_per_target: .r.repetitions[0].s_per_target,
  identity: ([.r.counts, .r.recovered, .r.all_verified] | tostring)
};
def table(rows): rows | map(parse + metrics) | group_by([.size, .path, .arm]) | map({
  size: .[0].size, path: .[0].path, arm: .[0].arm,
  collect_ms: (map(.collect_ms) | stats), build_ms: (map(.build_ms) | stats),
  total_ms: (map(.total_ms) | stats), total_units: (map(.total_units) | stats),
  s_per_target: (map(.s_per_target) | stats),
  identities: (map(.identity) | unique | length),
  identity: .[0].identity
});
(table($aa)) as $aat | (table($ab)) as $abt |
{
  aa: ($aat | map({size, path, collect_spread: (.collect_ms.max / .collect_ms.min - 1), total_spread: (.total_ms.max / .total_ms.min - 1), n: .collect_ms.n})),
  ab: ($abt | group_by([.size, .path]) | map(
      (map(select(.arm == "baseline"))[0]) as $b | (map(select(.arm == "candidate"))[0]) as $c | {
        size: $b.size, path: $b.path,
        identical_outputs: ($b.identities == 1 and $c.identities == 1 and $b.identity == $c.identity),
        baseline: ($b | del(.identity, .identities, .size, .path, .arm)),
        candidate: ($c | del(.identity, .identities, .size, .path, .arm)),
        collect_ratio: ($b.collect_ms.median / $c.collect_ms.median),
        build_ratio: ($b.build_ms.median / $c.build_ms.median),
        total_ratio: ($b.total_ms.median / $c.total_ms.median),
        total_units_ratio: ($b.total_units.median / $c.total_units.median)
      })),
  isolation: {runs: ($iso | length), contended: ($iso | map(select(.run.contended == true)) | length), nonzero_exit: ($iso | map(select(.run.exit_status != 0)) | length)}
}' > "$here/summary.json"
echo "wrote $here/summary.json"
