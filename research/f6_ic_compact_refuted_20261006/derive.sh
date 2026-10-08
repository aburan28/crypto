#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_compact_refuted_20261006
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
        direct_fused_pack:($r.effective_config.direct_fused_pack // false),
        active_multipliers:($r.effective_config.active_multipliers // false),
        support_local_stream:($r.effective_config.support_local_stream // false),
        support_local_bitmap_columns:($r.effective_config.support_local_bitmap_columns // false),
        support_local_profile:($r.effective_config.support_local_profile // false),
        inherited_basis_profile_enabled:($r.effective_config.inherited_basis_profile // false),
        compact_refuted:($r.effective_config.compact_refuted // false),
        scalar_verified:$r.scalar_verified, recovered_scalar:$r.solutions[0].recovered,
        target_point:$r.fixture.targets[0],
        online_wall_ns:$r.online_wall_ns,
        online_phases_ns:$r.generic_phase_timing.online_phases_ns,
        f4_stage_profile:$r.f4_stage_profile,
        support_local_build_profile:$r.support_local_build_profile,
        inherited_basis_profile:$r.inherited_basis_profile,
        inherited_unassigned_ns:($r.f4_stage_profile.build_ns -
          $r.inherited_basis_profile.root_ns -
          $r.inherited_basis_profile.specialise_total_ns),
        support_local_component_ns:($r.support_local_build_profile.rows_ns +
          $r.support_local_build_profile.columns_ns +
          $r.support_local_build_profile.pack_ns),
        f4_layout_hits_online:$r.f4_layout_hits_online,
        f4_layout_misses_online:$r.f4_layout_misses_online,
        f4_stage_ns:($r.f4_stage_profile.build_ns +
          $r.f4_stage_profile.reduce_ns + $r.f4_stage_profile.readback_ns),
        attempts:($r.solutions[0].attempts|length),
        attempt_outcomes:[$r.solutions[0].attempts[].pdp.outcome],
        oracle_calls:([$r.solutions[0].attempts[].pdp.stats.stats.geometric_oracle_calls]|add),
        oracle_ns:([$r.solutions[0].attempts[].pdp.stats.stats.geometric_oracle_ns]|add),
        geometric_refutations:([$r.solutions[0].attempts[].pdp.stats.stats.geometric_refutations]|add),
        geometric_witnesses:([$r.solutions[0].attempts[].pdp.stats.stats.geometric_witnesses]|add),
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
      ((.target == 1 and .recovered_scalar == "4785") or
       (.target == 7 and .recovered_scalar == "2391")) and
      .f4_stage_profile.calls > 0 and
      .f4_stage_ns <= .online_phases_ns.target_pdp and
      .oracle_ns <= .online_phases_ns.target_pdp and
      .oracle_calls > 0 and
      .direct_fused_pack == false and .active_multipliers == false and
      .support_local_stream == false and .support_local_bitmap_columns == false and
      .support_local_profile == false and .inherited_basis_profile_enabled == false and
      ((.arm == "compact_baseline" and .compact_refuted == false) or
       (.arm == "compact_candidate" and .compact_refuted == true)) and
      .inherited_unassigned_ns >= 0 and
      (.online_phases_ns.target_query + .online_phases_ns.target_pdp +
       .online_phases_ns.target_relation_check + .online_phases_ns.target_descent +
       .online_phases_ns.recovery_check) == .online_wall_ns
    ) and
    ([.[].target]|unique|sort) == [1,7] and
    (group_by([.target,.repetition]) | all(.[];
      length == 2 and
      ([.[].arm]|sort) == ["compact_baseline","compact_candidate"] and
      ([.[].target_point]|unique|length) == 1 and
      ([.[].recovered_scalar]|unique|length) == 1 and
      ([.[].attempts]|unique|length) == 1 and
      ([.[].attempt_outcomes]|unique|length) == 1 and
      ([.[].oracle_calls]|unique|length) == 1 and
      ([.[].geometric_refutations]|unique|length) == 1 and
      ([.[].geometric_witnesses]|unique|length) == 1 and
      ([.[].reductions]|unique|length) == 1 and
      ([.[].geometric_additions]|unique|length) == 1
    ))
' "$rows" > "$out/DERIVATION_CHECK.json"
jq -s -c '
  group_by([.target,.repetition])[] |
  (map(select(.arm == "compact_baseline"))[0]) as $b |
  (map(select(.arm == "compact_candidate"))[0]) as $c |
  {
    target:$b.target, repetition:$b.repetition,
    baseline_candidate_id:$b.candidate_id, candidate_candidate_id:$c.candidate_id,
    baseline_online_ms:($b.online_wall_ns/1000000),
    candidate_online_ms:($c.online_wall_ns/1000000),
    online_ratio:($b.online_wall_ns/$c.online_wall_ns),
    baseline_build_ms:($b.f4_stage_profile.build_ns/1000000),
    candidate_build_ms:($c.f4_stage_profile.build_ns/1000000),
    build_ratio:($b.f4_stage_profile.build_ns/$c.f4_stage_profile.build_ns),
    baseline_f4_word_ops:$b.f4_stage_profile.word_ops,
    candidate_f4_word_ops:$c.f4_stage_profile.word_ops,
    baseline_oracle_calls:$b.oracle_calls,
    candidate_oracle_calls:$c.oracle_calls,
    baseline_refutations:$b.geometric_refutations,
    candidate_refutations:$c.geometric_refutations
  }
' "$rows" > "$out/pairs.jsonl"
shasum -a 256 "$rows" "$out/DERIVATION_CHECK.json" "$out/pairs.jsonl" > "$out/DERIVED_SHA256SUMS"
