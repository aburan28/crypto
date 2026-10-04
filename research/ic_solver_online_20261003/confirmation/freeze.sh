#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
root=research/ic_solver_online_20261003
out="$root/confirmation"
state=research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/mathematics.json
binary=${IC_WORKER_BINARY:-/private/tmp/ic-solver-e2e-target/release/examples/ic_tournament_worker}
actual=$(shasum -a 256 "$binary" | awk '{print $1}')
expected=$(sed -n 's/^binary_sha256=//p' "$root/freeze.txt")
test "$actual" = "$expected"
test "$(jq -r '.fixture.target_seeds[0]' "$out/fixture.json")" = 20261003018
jq -cnS --slurpfile prior "$root/workload.json" --slurpfile fixture "$out/fixture.json" \
  '$prior[0].record | .algorithm_seed=20261003033 |
    .target_seeds=[20261003018] |
    .targets=[$fixture[0].fixture.targets[0] | map(tonumber)] |
    .question="fresh-target-production-f5-versus-inherited-f4-confirmation"' \
  > "$out/workload-record.json"
sha=$(jq -jcS . "$out/workload-record.json" | shasum -a 256 | awk '{print $1}')
jq -cnS --slurpfile record "$out/workload-record.json" --arg sha "$sha" \
  '{record:$record[0],record_sha256:$sha,workload_id:$sha[0:12]}' > "$out/workload.json"
for solver in f5 inherited_f4; do
  jq -cnS --arg solver "$solver" --slurpfile fixture "$out/fixture.json" \
    --slurpfile prepared "$state" \
    '{mode:"ic",degree:17,curve_a:1,target_seeds:[],
      public_targets:$fixture[0].fixture.targets,algorithm_seed:20261003033,
      factor_base:{kind:"standard_subspace",dimension:6},
      config:{solver:$solver,linear_algebra:"dense",summands:3,
        groebner_degree:3,node_budget:8192,batch_trials:1,max_trials:32,
        rho_parallel_walks:1},exclusive_phases:true,prepared:$prepared[0]}' \
    > "$out/input-$solver.json"
done
jq -cnS --slurpfile fixture "$out/fixture.json" \
  '{mode:"rho",degree:17,curve_a:1,target_seeds:[],
    public_targets:$fixture[0].fixture.targets,algorithm_seed:20261003033,
    factor_base:{kind:"standard_subspace",dimension:6},
    config:{rho_parallel_walks:1,max_trials:65536,batch_trials:1},
    exclusive_phases:true}' > "$out/input-rho.json"
printf '%s\n' "binary_sha256=$actual" "workload_id=$(jq -r .workload_id "$out/workload.json")" \
  "f5=$(jq -r .candidate_id "$root/f5-candidate.json")" \
  "inherited_f4=$(jq -r .candidate_id "$root/inherited_f4-candidate.json")" \
  > "$out/freeze.txt"
