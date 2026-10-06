#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_specialise_breakdown_20261006
rows="$out/measurements.jsonl"
: > "$rows"
while IFS="$(printf '\t')" read -r target repetition run_id input_sha candidate_sha workload_id; do
    [ "$target" = target ] && continue
    stem="$out/runs/T${target}-R${repetition}"
    rc=$(cat "$stem.exit-status")
    if jq -e . "$stem.stdout.json" > /dev/null 2>&1; then
        jq -n -c --arg run_id "$run_id" --arg workload_id "$workload_id" \
            --arg input_sha "$input_sha" --arg candidate_sha "$candidate_sha" \
            --argjson target "$target" --argjson repetition "$repetition" \
            --argjson rc "$rc" --slurpfile report "$stem.stdout.json" '
            $report[0] as $r | $r.specialise_profile as $p |
            ($p.layout_ns+$p.bookkeeping_ns+$p.rewrite_ns+$p.reduction_ns+
             $p.completion_ns+$p.closure_ns) as $component_sum |
            {target:$target,repetition:$repetition,run_id:$run_id,
             candidate_id:($run_id|split("W")[0]),workload_id:$workload_id,
             input_sha256:$input_sha,candidate_sha256:$candidate_sha,
             exit_status:$rc,status:$r.status,verified:$r.scalar_verified,
             recovered_scalar:$r.solutions[0].recovered,
             attempts:($r.solutions[0].attempts|length),
             online_wall_ns:$r.online_wall_ns,
             online_phases_ns:$r.generic_phase_timing.online_phases_ns,
             phase_sum_ns:($r.generic_phase_timing.online_phases_ns.target_query+
               $r.generic_phase_timing.online_phases_ns.target_pdp+
               $r.generic_phase_timing.online_phases_ns.target_relation_check+
               $r.generic_phase_timing.online_phases_ns.target_descent+
               $r.generic_phase_timing.online_phases_ns.recovery_check),
             f4_build_ns:$r.f4_stage_profile.build_ns,
             f4_word_ops:$r.f4_stage_profile.word_ops,
             inherited_basis_profile:$r.inherited_basis_profile,
             specialise_profile:$p,
             component_sum_ns:$component_sum,
             instrumentation_residual_ns:($p.total_ns-$component_sum),
             memory_peak_bytes:null,isolation:"unisolated-macos-exploratory"}' >> "$rows"
    else
        jq -n -c --arg run_id "$run_id" --arg workload_id "$workload_id" \
            --argjson target "$target" --argjson repetition "$repetition" \
            --argjson rc "$rc" \
            '{target:$target,repetition:$repetition,run_id:$run_id,
              workload_id:$workload_id,exit_status:$rc,status:"invalid_worker_json",
              verified:false,online_wall_ns:null}' >> "$rows"
    fi
done < "$out/runs/INDEX.tsv"
jq -s -c '{rows:length,all_verified:(all(.[]; .exit_status==0 and
  .status=="complete" and .verified==true)),
  all_phase_sums:(all(.[]; .phase_sum_ns==.online_wall_ns)),
  all_nested_sums:(all(.[]; .component_sum_ns<=.specialise_profile.total_ns)),
  expected_targets:([.[].target]|sort==[1,1,7,7]),
  exact_counted_controls:(all(.[];
    if .target==1 then .recovered_scalar=="4785" and .attempts==1 and
      .f4_word_ops==422488 and .specialise_profile.calls==35
    else .recovered_scalar=="2391" and .attempts==11 and
      .f4_word_ops==17167040 and .specialise_profile.calls==2470 end))}' \
  "$rows" > "$out/DERIVATION_CHECK.json"
shasum -a 256 "$rows" "$out/DERIVATION_CHECK.json" > "$out/DERIVED_SHA256SUMS"
