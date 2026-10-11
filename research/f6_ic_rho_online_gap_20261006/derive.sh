#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_rho_online_gap_20261006
rows="$out/measurements.jsonl"
: > "$rows"
while IFS="$(printf '\t')" read -r target repetition arm run_id workload input_sha point_sha; do
    [ "$target" = target ] && continue
    stem="$out/runs/T${target}-R${repetition}-${arm}"
    rc=$(cat "$stem.exit-status")
    if jq -e . "$stem.stdout.json" > /dev/null 2>&1; then
        jq -n -c \
          --arg run_id "$run_id" --arg workload "$workload" --arg arm "$arm" \
          --arg point_sha "$point_sha" --arg input_sha "$input_sha" \
          --argjson target "$target" --argjson repetition "$repetition" \
          --argjson rc "$rc" --slurpfile report "$stem.stdout.json" '
          $report[0] as $r |
          {
            target:$target,repetition:$repetition,arm:$arm,run_id:$run_id,
            workload_id:$workload,input_sha256:$input_sha,point_sha256:$point_sha,
            exit_status:$rc,status:$r.status,candidate_id:
              (if $arm=="ic" then "IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0h9628dfd41b76" else null end),
            target_point:$r.fixture.targets[0],recovered_scalar:$r.solutions[0].recovered,
            verified:(if $arm=="ic" then $r.scalar_verified
              else ($r.status=="complete" and $r.solutions[0].verified==true
                and $r.scalar_replay_included==true) end),
            online_wall_ns:$r.online_wall_ns,
            online_phases_ns:$r.generic_phase_timing.online_phases_ns,
            phase_sum_ns:(if $arm=="ic" then
                ($r.generic_phase_timing.online_phases_ns.target_query +
                 $r.generic_phase_timing.online_phases_ns.target_pdp +
                 $r.generic_phase_timing.online_phases_ns.target_relation_check +
                 $r.generic_phase_timing.online_phases_ns.target_descent +
                 $r.generic_phase_timing.online_phases_ns.recovery_check)
              else ($r.generic_phase_timing.online_phases_ns.rho_solve +
                    $r.generic_phase_timing.online_phases_ns.recovery_check) end),
            ic_attempts:(if $arm=="ic" then ($r.solutions[0].attempts|length) else null end),
            ic_f4_profile:(if $arm=="ic" then $r.f4_stage_profile else null end),
            rho_dispatch:(if $arm=="rho" then $r.rho_dispatch else null end),
            rho_iterations:(if $arm=="rho" then $r.solutions[0].iterations else null end),
            rho_restarts:(if $arm=="rho" then $r.solutions[0].restarts else null end),
            rho_walk_additions:(if $arm=="rho" then $r.solutions[0].walk_group_additions else null end),
            rho_effective_walks:(if $arm=="rho" then $r.solutions[0].effective_walks else null end),
            memory_peak_bytes:null,isolation:"unisolated-macos-exploratory"
          }' >> "$rows"
    else
        jq -n -c --arg run_id "$run_id" --arg workload "$workload" --arg arm "$arm" \
          --arg point_sha "$point_sha" --arg input_sha "$input_sha" \
          --argjson target "$target" --argjson repetition "$repetition" --argjson rc "$rc" \
          '{target:$target,repetition:$repetition,arm:$arm,run_id:$run_id,
            workload_id:$workload,input_sha256:$input_sha,point_sha256:$point_sha,
            exit_status:$rc,status:"invalid_worker_json",verified:false,
            online_wall_ns:null}' >> "$rows"
    fi
done < "$out/runs/INDEX.tsv"
jq -s -c '
  group_by([.target,.repetition])[] |
  (map(select(.arm=="ic"))[0]) as $ic |
  (map(select(.arm=="rho"))[0]) as $rho |
  (($ic.exit_status==0 and $rho.exit_status==0 and
    $ic.status=="complete" and $rho.status=="complete" and
    $ic.verified==true and $rho.verified==true and
    $ic.recovered_scalar==$rho.recovered_scalar and
    $ic.point_sha256==$rho.point_sha256 and
    $ic.phase_sum_ns==$ic.online_wall_ns and
    $rho.phase_sum_ns==$rho.online_wall_ns)) as $paired |
  {target:$ic.target,repetition:$ic.repetition,workload_id:$ic.workload_id,
    point_sha256:$ic.point_sha256,ic_run_id:$ic.run_id,rho_run_id:$rho.run_id,
    ic_status:$ic.status,rho_status:$rho.status,paired_verified:$paired,
    ic_online_ms:(if $paired then $ic.online_wall_ns/1000000 else null end),
    rho_online_ms:(if $paired then $rho.online_wall_ns/1000000 else null end),
    rho_over_ic:(if $paired then $rho.online_wall_ns/$ic.online_wall_ns else null end),
    recovered_scalar:(if $paired then $ic.recovered_scalar else null end),
    ic_attempts:$ic.ic_attempts,rho_iterations:$rho.rho_iterations,
    rho_restarts:$rho.rho_restarts,rho_walk_additions:$rho.rho_walk_additions}
' "$rows" > "$out/pairs.jsonl"
jq -s -c '{rows:length,all_phase_sums:(all(.[]; .phase_sum_ns==.online_wall_ns)),
  all_verified:(all(.[]; .exit_status==0 and .status=="complete" and .verified==true)),
  target_set:([.[].target]|unique|sort),arm_set:([.[].arm]|unique|sort),
  same_point_per_target:(group_by(.target)|all(.[];([.[].point_sha256]|unique|length)==1))}' \
  "$rows" > "$out/DERIVATION_CHECK.json"
shasum -a 256 "$rows" "$out/pairs.jsonl" "$out/DERIVATION_CHECK.json" > "$out/DERIVED_SHA256SUMS"
