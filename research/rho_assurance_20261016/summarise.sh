#!/usr/bin/env bash
# Reduce runs/ to summary.json: Q2's step-count ratios per size, Q3's
# instruction counts per walk step and per field operation, and Q4's
# instructions and native time per step for each arm.  Means, standard
# deviations, medians, extremes and their quotients only.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
runs="$here/runs"
collected() { grep -o 'Collected : [0-9]*' "$1" | grep -o '[0-9]*$'; }

q2=$(for f in "$runs"/q2/icv1-*.jsonl; do
    jq -sc --arg slug "$(basename "$f" .jsonl)" '
        def mean: add / length;
        def sd: mean as $m | (map((. - $m) * (. - $m)) | add / (length - 1)) | sqrt;
        (map(.walk_steps / .ideal_steps)) as $r |
        (map((.walk_steps + .walks) / .ideal_steps)) as $v |
        {slug: $slug, n: .[0].n, fixtures: length, all_verified: (map(.verified) | all),
         ideal_steps: .[0].ideal_steps, subgroup_order: .[0].subgroup_order, orbit_size: .[0].automorphism_size,
         ratio: {mean: ($r | mean), sd: ($r | sd),
                 ci95: [($r | mean) - 1.96 * ($r | sd) / ($r | length | sqrt), ($r | mean) + 1.96 * ($r | sd) / ($r | length | sqrt)],
                 median: ($r | sort | .[length / 2 | floor]), min: ($r | min), max: ($r | max)},
         visited_ratio: {mean: ($v | mean),
                 ci95: [($v | mean) - 1.96 * ($v | sd) / ($v | length | sqrt), ($v | mean) + 1.96 * ($v | sd) / ($v | length | sqrt)]},
         per_solve_means: {walks: (map(.walks) | mean), fruitless: (map(.fruitless_cycles) | mean),
                 capped: (map(.capped_walks) | mean), wasted_merges: (map(.wasted_merges) | mean)},
         steps_per_walk: ((map(.walk_steps) | add) / (map(.walks) | add)),
         ratios: $r}' "$f"
done | jq -sc .)

q3_cost() { # muls sqrs
    collected "$runs/q3/gf2-mul-cost-n53-$1-$2.valgrind.log"
}
base=$(q3_cost 1000000 1000000)
q3=$(jq -nc --argjson base "$base" --argjson m2 "$(q3_cost 2000000 1000000)" --argjson s2 "$(q3_cost 1000000 2000000)" \
    --argjson fix "$(collected "$runs/q3/fixture-n53-4.valgrind.log")" \
    --argjson steps "$(jq -s 'map(.walk_steps) | add' "$runs/q3/fixture-n53-4.jsonl")" \
    '{field: {ir_per_mul: (($m2 - $base) / 1e6), ir_per_sqr: (($s2 - $base) / 1e6)},
      fixture_n53_4: {ir: $fix, walk_steps: $steps, ir_per_step: ($fix / $steps)}}
     | .fixture_n53_4.step_in_mul_units = (.fixture_n53_4.ir_per_step / .field.ir_per_mul)')

q4_ir=$(for arm in baseline candidate candidate2; do
    f="$runs/q4/callgrind/fixture-n53-4-$arm"
    jq -nc --arg arm "$arm" --argjson ir "$(collected "$f.valgrind.log")" \
        --argjson steps "$(jq -s 'map(.walk_steps) | add' "$f.jsonl")" \
        --arg id "$(jq -c 'del(.peak_rss_bytes) | with_entries(select(.key | endswith("_ms") | not))' "$f.jsonl" | sha256sum | cut -c1-16)" \
        '{arm: $arm, ir: $ir, walk_steps: $steps, ir_per_step: ($ir / $steps), output_digest: $id}'
done | jq -sc '. as $a | {arms: $a, identical: ($a | map(.output_digest) | unique | length == 1),
    final_ratio: ($a[0].ir / $a[2].ir), first_part_ratio: ($a[0].ir / $a[1].ir)}')

q4_time=null
if ls "$runs"/q4/time/*/*.jsonl > /dev/null 2>&1; then
    contended=$(cat "$runs"/q4/time/*/isolation.jsonl | jq -sc 'map(select(.run.contended == true) | .label)')
    q4_time=$(for f in "$runs"/q4/time/*/*.jsonl; do
        jq -sc --arg file "$(basename "$f")" --arg phase "$(basename "$(dirname "$f")")" --argjson contended "$contended" '
            {phase: $phase, file: $file,
             contended: (("q4/" + $phase + "/" + ($file | sub("\\.jsonl$"; ""))) as $l | $contended | index($l) != null),
             arm: ($file | capture("-(?<a>baseline|candidate)-r").a),
             slug: ($file | capture("^(?<s>icv1-[^-]+-[^-]+-[0-9a-f]{8})").s),
             ns_per_step: ((map(.walk_ms) | add) * 1e6 / (map(.walk_steps) | add)),
             identity: (map(del(.peak_rss_bytes) | with_entries(select(.key | endswith("_ms") | not))) | tostring)}' "$f"
    done | jq -sc '
        def med: sort | if length % 2 == 1 then .[length / 2 | floor] else (.[length / 2 - 1] + .[length / 2]) / 2 end;
        map(select(.contended | not)) |
        group_by(.slug) | map({slug: .[0].slug,
            aa_spread: (map(select(.phase == "aa") | .ns_per_step) | (max / min - 1)),
            baseline_ns: (map(select(.phase == "ab" and .arm == "baseline") | .ns_per_step) | {median: med, min: min, max: max, n: length}),
            candidate_ns: (map(select(.phase == "ab" and .arm == "candidate") | .ns_per_step) | {median: med, min: min, max: max, n: length}),
            identical: (map(.identity) | unique | length == 1)}
          | .ratio = (.baseline_ns.median / .candidate_ns.median))')
    iso=$(cat "$runs"/q4/time/*/isolation.jsonl | jq -sc '{runs: length,
        contended: map(select(.run.contended == true) | .label), nonzero_exit: (map(select(.run.exit_status != 0)) | length)}')
    q4_time=$(jq -nc --argjson t "$q4_time" --argjson iso "$iso" '{cells: $t, isolation: $iso}')
fi

jq -n --argjson q2 "$q2" --argjson q3 "$q3" --argjson q4_ir "$q4_ir" --argjson q4_time "$q4_time" \
    '{q2_step_counts: $q2, q3_step_cost: $q3, q4_instructions: $q4_ir, q4_native_time: $q4_time}' > "$here/summary.json"
echo "wrote $here/summary.json"
