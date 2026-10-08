#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_generic_direct_pack_20261006
rows="$out/measurements.jsonl"
: > "$rows"
while IFS="$(printf '\t')" read -r target repetition arm run_id input_sha candidate_sha workload_id; do
    [ "$target" = target ] && continue
    stem="$out/runs/T${target}-R${repetition}-${arm}"
    rc=$(cat "$stem.exit-status")
    jq -n -c \
        --arg run_id "$run_id" --arg workload_id "$workload_id" \
        --arg candidate_id "${run_id%%W*}" --arg arm "$arm" \
        --argjson target "$target" --argjson repetition "$repetition" \
        --argjson exit_status "$rc" --arg input_sha256 "$input_sha" \
        --arg candidate_sha256 "$candidate_sha" \
        --slurpfile report "$stem.stdout.json" '
      $report[0] as $r |
      {
        candidate_id:$candidate_id, workload_id:$workload_id, run_id:$run_id,
        target:$target, repetition:$repetition, arm:$arm,
        input_sha256:$input_sha256, candidate_sha256:$candidate_sha256,
        exit_status:$exit_status, status:$r.status,
        direct_fused_pack:$r.effective_config.direct_fused_pack,
        scalar_verified:$r.scalar_verified, recovered_scalar:$r.solutions[0].recovered,
        target_point:$r.fixture.targets[0],
        online_wall_ns:$r.online_wall_ns,
        online_phases_ns:$r.generic_phase_timing.online_phases_ns,
        f4_stage_profile:$r.f4_stage_profile,
        f4_layout_hits_online:$r.f4_layout_hits_online,
        f4_layout_misses_online:$r.f4_layout_misses_online,
        f4_stage_ns:($r.f4_stage_profile.build_ns +
          $r.f4_stage_profile.reduce_ns + $r.f4_stage_profile.readback_ns),
        attempts:($r.solutions[0].attempts|length),
        oracle_calls:([$r.solutions[0].attempts[].pdp.stats.stats.geometric_oracle_calls]|add),
        oracle_ns:([$r.solutions[0].attempts[].pdp.stats.stats.geometric_oracle_ns]|add),
        geometric_additions:([$r.solutions[0].attempts[].pdp.stats.stats.geometric_group_additions]|add),
        pair_index_builds:([$r.solutions[0].attempts[].pdp.stats.stats.geometric_pair_index_builds]|add),
        pair_index_lookups:([$r.solutions[0].attempts[].pdp.stats.stats.geometric_pair_index_lookups]|add),
        reductions:([$r.solutions[0].attempts[].pdp.stats.stats.reductions]|add),
        memory_peak_bytes:null,
        isolation:"unisolated-macos-exploratory"
      }' >> "$rows"
done < "$out/runs/INDEX.tsv"
jq -e -s '
    length == 8 and all(.[];
      .exit_status == 0 and .status == "complete" and .scalar_verified == true and
      .f4_stage_profile.calls > 0 and
      .f4_stage_ns <= .online_phases_ns.target_pdp and
      .oracle_ns <= .online_phases_ns.target_pdp and
      .oracle_calls > 0 and
      (.direct_fused_pack == (.arm|endswith("_direct"))) and
      (.online_phases_ns.target_query + .online_phases_ns.target_pdp +
       .online_phases_ns.target_relation_check + .online_phases_ns.target_descent +
       .online_phases_ns.recovery_check) == .online_wall_ns
    ) and
    (group_by(.target) | all(.[];
      ([.[].target_point]|unique|length) == 1 and
      ([.[].recovered_scalar]|unique|length) == 1 and
      ([.[].attempts]|unique|length) == 1
    )) and
    (group_by(.target,.repetition,(.arm|sub("_direct$";""))) | all(.[];
      length == 2 and
      ([.[].direct_fused_pack]|sort) == [false,true] and
      ([.[].reductions]|unique|length) == 1 and
      ([.[].geometric_additions]|unique|length) == 1 and
      ([.[].f4_layout_hits_online]|unique|length) == 1 and
      ([.[].f4_layout_misses_online]|unique|length) == 1 and
      ([.[].f4_stage_profile.calls]|unique|length) == 1 and
      ([.[].f4_stage_profile.rows]|unique|length) == 1 and
      ([.[].f4_stage_profile.cols]|unique|length) == 1 and
      ([.[].f4_stage_profile.word_ops]|unique|length) == 1
    ))
' "$rows" > "$out/DERIVATION_CHECK.json"
shasum -a 256 "$rows" "$out/DERIVATION_CHECK.json" > "$out/DERIVED_SHA256SUMS"
