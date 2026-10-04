#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
root=research/ic_solver_online_20261003
out="$root/confirmation"
workload=$(jq -r .workload_id "$out/workload.json")
: > "$out/measurements.jsonl"
: > "$out/rho-reference.jsonl"
for repeat in 1 2 3; do
  for arm in f5 inherited_f4; do
    id=$(jq -r .candidate_id "$root/$arm-candidate.json")
    stem="$out/runs/R$repeat-$arm"
    status=$(cat "$stem.exit-status")
    jq -c --arg id "$id" --arg workload "$workload" \
      --argjson repeat "$repeat" --argjson exit "$status" \
      '{candidate_id:$id,workload_id:$workload,
        run_id:($id+"W"+$workload+"R"+($repeat|tostring)),
        repetition:$repeat,process_exit:$exit,status:.status,
        online_wall_ns:.online_wall_ns,
        online_phases_ns:.generic_phase_timing.online_phases_ns,
        recovered:.solutions[0].recovered,scalar_verified:.scalar_verified,
        attempts:(.solutions[0].attempts|length),
        outcome_mix:(.solutions[0].attempts|map(.pdp.outcome)|group_by(.)|map({outcome:.[0],count:length})),
        usable_factor_base_points:62,folded_columns:.columns,
        peak_rss_bytes:null}' "$stem.stdout.json" >> "$out/measurements.jsonl"
  done
  stem="$out/runs/R$repeat-rho"
  status=$(cat "$stem.exit-status")
  jq -c --arg workload "$workload" --argjson repeat "$repeat" \
    --argjson exit "$status" \
    '{reference_id:"rho-signed-frobenius-n17a1-one-walk",
      workload_id:$workload,repetition:$repeat,process_exit:$exit,
      status:.status,online_wall_ns:.online_wall_ns,
      online_phases_ns:.generic_phase_timing.online_phases_ns,
      recovered:.solutions[0].recovered,
      scalar_verified:.solutions[0].verified,
      iterations:.solutions[0].iterations,
      walk_group_additions:.solutions[0].walk_group_additions,
      peak_rss_bytes:null}' "$stem.stdout.json" >> "$out/rho-reference.jsonl"
done

# Admission is strict: no absent ledger field is interpreted as zero.
jq -es 'length == 6 and all(.[];
  .process_exit == 0 and .status == "complete" and .scalar_verified == true and
  .recovered == "852" and .attempts == 17 and
  (.online_phases_ns | [.target_query,.target_pdp,.target_relation_check,.target_descent,.recovery_check] | all(. != null)) and
  .online_wall_ns == ([.online_phases_ns.target_query,.online_phases_ns.target_pdp,
    .online_phases_ns.target_relation_check,.online_phases_ns.target_descent,
    .online_phases_ns.recovery_check] | add))' \
  "$out/measurements.jsonl" > "$out/ic-admission.json"
jq -es 'length == 3 and all(.[];
  .process_exit == 0 and .status == "complete" and .scalar_verified == true and
  .recovered == "852" and .online_wall_ns > 0)' \
  "$out/rho-reference.jsonl" > "$out/rho-admission.json"
