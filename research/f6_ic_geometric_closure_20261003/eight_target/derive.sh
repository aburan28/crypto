#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
out=research/f6_ic_geometric_closure_20261003/eight_target
base=research/f6_ic_geometric_closure_20261003/three_way
test ! -e "$out/measurements.jsonl"
for index in 1 2 3 4 5 6 7 8; do
  workload=$(jq -r .workload_id "$out/T$index/workload.json")
  for arm in inherited_f4 f5 f6_ic; do
    candidate=$(jq -r .candidate_id "$base/$arm-candidate.json")
    if [ "$arm" = f5 ]; then reps='1'; else reps='1 2'; fi
    for repeat in $reps; do
      stem="$out/runs/T${index}-R${repeat}-${arm}"
      exit_status=$(cat "$stem.exit-status")
      raw_sha=$(shasum -a 256 "$stem.stdout.json" | awk '{print $1}')
      jq -cnS --slurpfile raw "$stem.stdout.json" \
        --arg arm "$arm" --arg candidate "$candidate" \
        --arg workload "$workload" --arg raw_sha "$raw_sha" \
        --argjson index "$index" --argjson repeat "$repeat" \
        --argjson exit_status "$exit_status" \
        '{target_index:$index,arm:$arm,candidate_id:$candidate,
          workload_id:$workload,
          run_id:($candidate+"W"+$workload+"R"+($repeat|tostring)),
          repeat:$repeat,exit_status:$exit_status,raw_sha256:$raw_sha,
          cpu_isolation:"L0-unverified-mac",memory_peak_bytes:null,
          status:$raw[0].status,scalar_verified:$raw[0].scalar_verified,
          recovered:$raw[0].solutions[0].recovered,
          target:($raw[0].fixture.targets[0]|map(tonumber)),
          online_wall_ns:$raw[0].online_wall_ns,
          online_phases_ns:$raw[0].generic_phase_timing.online_phases_ns,
          phase_sum_ns:($raw[0].generic_phase_timing.online_phases_ns |
            [.target_query,.target_pdp,.target_relation_check,.target_descent,.recovery_check]|add),
          attempts:$raw[0].solutions[0].attempts}' \
        >> "$out/measurements.jsonl"
    done
  done
done
jq -es 'length==40 and all(.[];
  .exit_status==0 and .status=="complete" and .scalar_verified==true
  and .online_wall_ns==.phase_sum_ns)' "$out/measurements.jsonl" >/dev/null
jq -sS '
  group_by(.target_index) | map(
    . as $rows |
    ($rows|map(select(.arm=="inherited_f4"))) as $f4 |
    ($rows|map(select(.arm=="f6_ic"))) as $f6 |
    ($rows|map(select(.arm=="f5"))) as $f5 |
    ($f4|map(.online_wall_ns)|add/length) as $f4_ns |
    ($f6|map(.online_wall_ns)|add/length) as $f6_ns |
    {target_index:$rows[0].target_index,
      workload_id:$rows[0].workload_id,target:$rows[0].target,
      recovered:$rows[0].recovered,
      attempts:($rows[0].attempts|length),
      f4_online_ns:$f4_ns,f5_online_ns:$f5[0].online_wall_ns,
      f6_online_ns:$f6_ns,
      f4_over_f6:($f4_ns/$f6_ns),
      f5_over_f6:($f5[0].online_wall_ns/$f6_ns),
      f4_aa_spread_ratio:((($f4|map(.online_wall_ns)|max)/($f4|map(.online_wall_ns)|min))-1),
      f4_reductions:([$f4[0].attempts[].pdp.stats.stats.reductions]|add),
      f5_reductions:([$f5[0].attempts[].pdp.stats.stats.reductions]|add),
      f6_reductions:([$f6[0].attempts[].pdp.stats.stats.reductions]|add),
      f6_group_additions:([$f6[0].attempts[].pdp.stats.stats.geometric_group_additions]|add),
      f6_residual_lookups:([$f6[0].attempts[].pdp.stats.stats.geometric_residual_lookups]|add),
      f6_batch_groups:([$f6[0].attempts[].pdp.stats.stats.geometric_batch_groups]|add)}
  )' "$out/measurements.jsonl" > "$out/target-summary.json"
jq -e 'length==8 and all(.[];.f4_over_f6>1 and .f5_over_f6>1)' "$out/target-summary.json" >/dev/null
