#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
out=research/f6_ic_geometric_closure_20261003/three_way
test ! -e "$out/measurements.jsonl"
workload=$(jq -r .workload_id "$out/workload.json")
for repeat in 1 2 3; do
  for arm in inherited_f4 f5 f6_ic; do
    stem="$out/runs/R${repeat}-${arm}"
    candidate=$(jq -r .candidate_id "$out/$arm-candidate.json")
    exit_status=$(cat "$stem.exit-status")
    raw_sha=$(shasum -a 256 "$stem.stdout.json" | awk '{print $1}')
    jq -cnS --slurpfile raw "$stem.stdout.json" \
      --arg arm "$arm" --arg candidate "$candidate" \
      --arg workload "$workload" --arg raw_sha "$raw_sha" \
      --argjson repeat "$repeat" --argjson exit_status "$exit_status" \
      '{arm:$arm,candidate_id:$candidate,workload_id:$workload,
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
        preparation_phases_ns:$raw[0].generic_phase_timing.phases_ns,
        attempts:$raw[0].solutions[0].attempts}' \
      >> "$out/measurements.jsonl"
  done
done
jq -es 'length==9 and
  all(.[]; .exit_status==0 and .status=="complete" and .scalar_verified==true
    and .recovered=="32917" and .target==[79391,5777]
    and .online_wall_ns==.phase_sum_ns)' "$out/measurements.jsonl" >/dev/null
