#!/usr/bin/env bash
# Reduce runs/ to summary.json: Callgrind instructions per walk step for each
# arm (with an output digest per arm), and native ns per walk step per curve
# and arm from the uncontended A/A and A/B runs.  Medians, extremes and their
# quotients only.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
runs="$here/runs"
collected() { grep -o 'Collected : [0-9]*' "$1" | grep -o '[0-9]*$'; }
digest() { jq -c 'del(.peak_rss_bytes) | with_entries(select(.key | endswith("_ms") | not))' "$1" | sha256sum | cut -c1-16; }

ir=$(for arm in baseline candidate1-4 candidate final; do
    f="$runs/callgrind/fixture-n53-4-$arm"
    [ -e "$f.jsonl" ] || continue
    jq -nc --arg arm "$arm" --argjson ir "$(collected "$f.valgrind.log")" \
        --argjson steps "$(jq -s 'map(.walk_steps) | add' "$f.jsonl")" --arg id "$(digest "$f.jsonl")" \
        '{arm: $arm, ir: $ir, walk_steps: $steps, ir_per_step: ($ir / $steps), output_digest: $id}'
done | jq -sc '. as $a | {arms: $a, identical: ($a | map(.output_digest) | unique | length == 1),
    baseline_over_final: (($a | map(select(.arm == "baseline"))[0].ir) / ($a | map(select(.arm == "final"))[0].ir))}')

time=null
if ls "$runs"/time/*/*.jsonl > /dev/null 2>&1; then
    contended=$(cat "$runs"/time/*/isolation.jsonl | jq -sc 'map(select(.run.contended == true) | .label)')
    time=$(for f in "$runs"/time/*/*.jsonl; do
        jq -sc --arg file "$(basename "$f")" --arg phase "$(basename "$(dirname "$f")")" --argjson contended "$contended" '
            {phase: $phase, file: $file,
             contended: (("time/" + $phase + "/" + ($file | sub("\\.jsonl$"; ""))) as $l | $contended | index($l) != null),
             arm: ($file | capture("-(?<a>baseline|candidate)-r").a),
             slug: ($file | capture("^(?<s>icv1-[^-]+-[^-]+-[0-9a-f]{8})").s),
             ns_per_step: ((map(.walk_ms) | add) * 1e6 / (map(.walk_steps) | add)),
             identity: (map(del(.peak_rss_bytes) | with_entries(select(.key | endswith("_ms") | not))) | tostring)}' "$f"
    done | jq -sc '
        def med: sort | if length % 2 == 1 then .[length / 2 | floor] else (.[length / 2 - 1] + .[length / 2]) / 2 end;
        (map(.identity) | group_by(.) | length) as $ids |
        map(select(.contended | not)) | group_by(.slug) | map({slug: .[0].slug,
            aa_spread: (map(select(.phase == "aa") | .ns_per_step) | (max / min - 1)),
            baseline_ns: (map(select(.phase == "ab" and .arm == "baseline") | .ns_per_step) | {median: med, min: min, max: max, n: length}),
            candidate_ns: (map(select(.phase == "ab" and .arm == "candidate") | .ns_per_step) | {median: med, min: min, max: max, n: length}),
            identical: (map(.identity) | unique | length == 1)}
          | .ratio = (.baseline_ns.median / .candidate_ns.median))')
    iso=$(cat "$runs"/time/*/isolation.jsonl | jq -sc '{runs: length,
        contended: map(select(.run.contended == true) | .label), nonzero_exit: (map(select(.run.exit_status != 0)) | length)}')
    time=$(jq -nc --argjson t "$time" --argjson iso "$iso" '{cells: $t, isolation: $iso}')
fi

jq -n --argjson ir "$ir" --argjson time "$time" '{instructions: $ir, native_time: $time}' > "$here/summary.json"
echo "wrote $here/summary.json"
