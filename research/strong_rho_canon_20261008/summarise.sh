#!/usr/bin/env bash
# Reduce runs/ to summary.json: for M1, M2 and M3, each arm's median, min and
# max, the A/A spread (max/min - 1 of the baseline's A/A runs), the A/B ratio
# of medians (baseline / candidate), and the identity check the protocol
# declares.  Computes nothing beyond medians, extremes and their quotients.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
runs="$here/runs"
rows() { # dir glob -> one JSON object per file: {file, phase, r}
    for phase in aa ab; do
        for f in "$runs/$1/$phase"/$2; do
            [ -e "$f" ] || continue
            jq -c -s --arg f "$(basename "$f")" --arg p "$phase" '{file: $f, phase: $p, r: .}' "$f"
        done
    done
}
jq -n \
  --slurpfile m1 <(rows m1 'icv1-*.json') \
  --slurpfile m2 <(rows m2 'icv1-*.jsonl') \
  --slurpfile m3 <(rows m3 'icv1-*.price.json') \
  --slurpfile iso <(cat "$runs"/m*/*/isolation.jsonl) '
def med: sort | if length % 2 == 1 then .[length/2|floor] else (.[length/2-1] + .[length/2]) / 2 end;
def stats: {median: med, min: min, max: max, n: length};
def spread: (max / min - 1);
($iso | map({key: .label, value: .run}) | from_entries) as $run |
def wall($l): $run[$l].wall_seconds;
# M1: one JSON object per file.
($m1 | map(.file | capture("^(?<slug>icv1-[^.]+)-(?<arm>baseline|candidate)-r(?<round>[0-9]+)")) ) as $m1k |
([$m1, $m1k] | transpose | map(.[0] + .[1] | {phase, arm, canon_ns: .r[0].canonicalize_ns_per_state,
   walk_ms: (.r[0].walk.wall_ns / 1e6), digest: .r[0].canonical_digest, capped: .r[0].walk.reached_cap})) as $m1r |
# M2: one JSON line per fixture; identity drops the RSS and every *_ms
# timing field, the only fields that vary between runs of one build.
([$m2, ($m2 | map(.file | capture("^(?<slug>icv1-[a-z0-9]+-[a-z0-9]+-[0-9a-f]{8})-(?<arm>baseline|candidate)-r(?<round>[0-9]+)")))] | transpose
   | map(.[0] + .[1] | {phase, slug, arm,
       wall_s: wall("m2/\(.phase)/\(.slug)-\(.arm)-r\(.round)"),
       walk_ms: (.r | map(.walk_ms) | add),
       identity: (.r | map(del(.peak_rss_bytes) | with_entries(select(.key | endswith("_ms") | not))) | tostring)})) as $m2r |
# M3: the pricing reports.
([$m3, ($m3 | map(.file | capture("^(?<slug>icv1-[a-z0-9]+-[a-z0-9]+-[0-9a-f]{8})-(?<path>portable|avx512)-(?<arm>baseline|candidate)-r(?<round>[0-9]+)")))] | transpose
   | map(.[0] + .[1] | {phase, slug, path, arm,
       collect_ms: (.r[0].repetitions[0].phases_ns.collect / 1e6),
       total_ms: (.r[0].repetitions[0].total_ns / 1e6),
       identity: ([.r[0].counts, .r[0].recovered, .r[0].all_verified] | tostring)})) as $m3r |
{
  m1: ({
    aa_spread: {canon: ($m1r | map(select(.phase == "aa") | .canon_ns) | spread),
                walk: ($m1r | map(select(.phase == "aa") | .walk_ms) | spread)},
    baseline: {canon_ns: ($m1r | map(select(.phase == "ab" and .arm == "baseline") | .canon_ns) | stats),
               walk_ms: ($m1r | map(select(.phase == "ab" and .arm == "baseline") | .walk_ms) | stats)},
    candidate: {canon_ns: ($m1r | map(select(.phase == "ab" and .arm == "candidate") | .canon_ns) | stats),
                walk_ms: ($m1r | map(select(.phase == "ab" and .arm == "candidate") | .walk_ms) | stats)},
    identical: (($m1r | map(.digest) | unique | length) == 1 and ($m1r | all(.capped)))
  } | .canon_ratio = (.baseline.canon_ns.median / .candidate.canon_ns.median)
    | .walk_ratio = (.baseline.walk_ms.median / .candidate.walk_ms.median)),
  m2: ($m2r | group_by(.slug) | map({
    slug: .[0].slug,
    aa_spread: (map(select(.phase == "aa") | .wall_s) | spread),
    baseline_s: (map(select(.phase == "ab" and .arm == "baseline") | .wall_s) | stats),
    candidate_s: (map(select(.phase == "ab" and .arm == "candidate") | .wall_s) | stats),
    walk_aa_spread: (map(select(.phase == "aa") | .walk_ms) | spread),
    baseline_walk_ms: (map(select(.phase == "ab" and .arm == "baseline") | .walk_ms) | stats),
    candidate_walk_ms: (map(select(.phase == "ab" and .arm == "candidate") | .walk_ms) | stats),
    identical: (map(.identity) | unique | length == 1)
  } | .ratio = (.baseline_s.median / .candidate_s.median)
    | .walk_ratio = (.baseline_walk_ms.median / .candidate_walk_ms.median))),
  m3: ($m3r | group_by([.slug, .path]) | map({
    slug: .[0].slug, path: .[0].path,
    aa_spread: {collect: (map(select(.phase == "aa") | .collect_ms) | spread),
                total: (map(select(.phase == "aa") | .total_ms) | spread)},
    baseline: {collect_ms: (map(select(.phase == "ab" and .arm == "baseline") | .collect_ms) | stats),
               total_ms: (map(select(.phase == "ab" and .arm == "baseline") | .total_ms) | stats)},
    candidate: {collect_ms: (map(select(.phase == "ab" and .arm == "candidate") | .collect_ms) | stats),
                total_ms: (map(select(.phase == "ab" and .arm == "candidate") | .total_ms) | stats)},
    identical: (map(.identity) | unique | length == 1)
  } | .collect_ratio = (.baseline.collect_ms.median / .candidate.collect_ms.median)
    | .total_ratio = (.baseline.total_ms.median / .candidate.total_ms.median))),
  isolation: {runs: ($iso | length), contended: ($iso | map(select(.run.contended == true) | .label)),
              nonzero_exit: ($iso | map(select(.run.exit_status != 0)) | length)}
}' > "$here/summary.json"
echo "wrote $here/summary.json"
