#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
out=research/f6_ic_geometric_closure_20261003/small_cold
test ! -e "$out/measurements.jsonl"
workload=$(jq -r .workload_id "$out/workload.json")
for arm in inherited_f4 f6_ic f5; do
  stem="$out/runs/R1-$arm"
  candidate=$(jq -r .candidate_id "$out/$arm-candidate.json")
  raw_sha=$(shasum -a 256 "$stem.stdout.json" | awk '{print $1}')
  exit_status=$(cat "$stem.exit-status")
  jq -cnS --slurpfile raw "$stem.stdout.json" \
    --arg arm "$arm" --arg candidate "$candidate" \
    --arg workload "$workload" --arg raw_sha "$raw_sha" \
    --argjson exit_status "$exit_status" \
    '{arm:$arm,candidate_id:$candidate,workload_id:$workload,
      run_id:($candidate+"W"+$workload+"R1"),
      exit_status:$exit_status,raw_sha256:$raw_sha,
      cpu_isolation:"L0-unverified-mac",memory_peak_bytes:null,
      status:$raw[0].status,recovered:$raw[0].solutions[0].recovered,
      target:($raw[0].fixture.targets[0]|map(tonumber)),
      actual_usable_base_points:14,folded_columns:$raw[0].columns,
      ordinary_queries:$raw[0].trials,
      verified_relations:$raw[0].accepted_relations,
      duplicate_relations:$raw[0].duplicate_relations,
      rejected_relations:$raw[0].rejected_relations,
      rank_certified:$raw[0].log_table_report.verified,
      certified_column_logs:($raw[0].column_logs|length),
      cold_inside_worker_ns:$raw[0].generic_phase_timing.observed_wall_ns,
      cold_phase_sum_ns:($raw[0].generic_phase_timing.phases_ns|[.[]|select(.!=null)]|add),
      cold_phases_ns:$raw[0].generic_phase_timing.phases_ns,
      online_wall_ns:$raw[0].online_wall_ns,
      online_phase_sum_ns:($raw[0].generic_phase_timing.online_phases_ns|[.[]|select(.!=null)]|add),
      online_phases_ns:$raw[0].generic_phase_timing.online_phases_ns,
      target_attempts:$raw[0].solutions[0].attempts,
      scalar_replay_included:$raw[0].scalar_replay_included}' \
    >> "$out/measurements.jsonl"
done
jq -es 'length==3 and all(.[];
  .exit_status==0 and .status=="complete" and .recovered=="4"
  and .target==[305,466] and .rank_certified==true
  and .scalar_replay_included==true
  and .cold_inside_worker_ns==.cold_phase_sum_ns
  and .online_wall_ns==.online_phase_sum_ns)' "$out/measurements.jsonl" >/dev/null
