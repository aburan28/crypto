#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
out=research/f6_ic_geometric_closure_20261003/pilot
test ! -e "$out/freeze.txt"
state=research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/mathematics.json
prior=research/ic_solver_online_20261003/inherited_f4-record.json
binary=${IC_WORKER_BINARY:-/private/tmp/f6-ic-target/release/examples/ic_tournament_worker}
digest() { shasum -a 256 "$1" | awk '{print $1}'; }
commit=$(git rev-parse HEAD)
worker_sha=$(digest examples/ic_tournament_worker.rs)
groebner_sha=$(digest src/cryptanalysis/koblitz_groebner.rs)
ic_sha=$(digest src/cryptanalysis/koblitz_index_calculus.rs)
measurement_sha=$(digest src/cryptanalysis/ic_measurement.rs)
state_sha=$(digest "$state")
binary_sha=$(digest "$binary")
test "$(jq -r '.fixture.target_seeds[0]' "$out/fixture.json")" = 20261003019
test "$(jq -r '.fixture.target_scalar_constructed' "$out/fixture.json")" = false
test "$(jq -r '.fixture.targets | length' "$out/fixture.json")" = 1

jq -cS --slurpfile fixture "$out/fixture.json" \
  '.algorithm_seed=20261003034 |
   .target_seeds=[20261003019] |
   .targets=[$fixture[0].fixture.targets[0] | map(tonumber)] |
   .question="prepared-one-target-inherited-f4-versus-f6-ic-pilot"' \
  research/ic_solver_online_20261003/confirmation/workload-record.json \
  > "$out/workload-record.json"
workload_sha=$(jq -jcS . "$out/workload-record.json" | shasum -a 256 | awk '{print $1}')
jq -cnS --slurpfile record "$out/workload-record.json" --arg sha "$workload_sha" \
  '{record:$record[0],record_sha256:$sha,workload_id:$sha[0:12]}' \
  > "$out/workload.json"

for solver in inherited_f4 f6_ic; do
  jq -cS --arg solver "$solver" --arg commit "$commit" \
    --arg worker_sha "$worker_sha" --arg groebner_sha "$groebner_sha" \
    --arg ic_sha "$ic_sha" --arg measurement_sha "$measurement_sha" \
    --arg binary_sha "$binary_sha" --arg state_sha "$state_sha" \
    '.point_decomposition.solver=$solver |
     .point_decomposition.source_sha256=$groebner_sha |
     .point_decomposition.limits.split_rule="highest-free" |
     .implementation={code_snapshot_commit:$commit,
       components:[{role:"ic-worker",sha256:$worker_sha},
                   {role:"boolean-groebner",sha256:$groebner_sha},
                   {role:"ic-engine",sha256:$ic_sha},
                   {role:"exclusive-measurement",sha256:$measurement_sha},
                   {role:"prepared-log-state",sha256:$state_sha}],
       flags:{binary_sha256:$binary_sha,rayon_threads:1,
         exclusive_phases:true,
         preparation_mode:"imported-certified-log-table-v1",
         reusable_template:"fixed-zero-symbolic-specialization-before-online",
         f6_ic:($solver=="f6_ic")}} |
     if $solver=="f6_ic" then
       .point_decomposition.geometric_oracle={
         version:0,
         support:"exact-usable-factor-base-x-index",
         closure:"m-minus-one-fixed-x-exact-residual-group-witness",
         fallback:"inherited-f4-if-encoding-unproved",
         oracle_cost:"charged-to-target-pdp"}
     else . end' "$prior" > "$out/$solver-record.json"
  sha=$(jq -jcS . "$out/$solver-record.json" | shasum -a 256 | awk '{print $1}')
  if [ "$solver" = f6_ic ]; then stage=f6; else stage=f4; fi
  id="IC1N17Ckb1fb62PDP3${stage}RCsampleLAgaussTDpdpISO0h$(printf '%.12s' "$sha")"
  jq -cnS --slurpfile record "$out/$solver-record.json" --arg sha "$sha" --arg id "$id" \
    '{record:$record[0],record_sha256:$sha,candidate_id:$id}' \
    > "$out/$solver-candidate.json"
  jq -cnS --arg solver "$solver" --slurpfile fixture "$out/fixture.json" \
    --slurpfile prepared "$state" \
    '{mode:"ic",degree:17,curve_a:1,target_seeds:[],
      public_targets:$fixture[0].fixture.targets,algorithm_seed:20261003034,
      factor_base:{kind:"standard_subspace",dimension:6},
      config:{solver:$solver,linear_algebra:"dense",summands:3,
        groebner_degree:3,node_budget:8192,batch_trials:1,max_trials:32,
        rho_parallel_walks:1},exclusive_phases:true,prepared:$prepared[0]}' \
    > "$out/input-$solver.json"
done

printf '%s\n' "source_commit=$commit" "binary_sha256=$binary_sha" \
  "fixture_sha256=$(digest "$out/fixture.json")" \
  "input_inherited_f4_sha256=$(digest "$out/input-inherited_f4.json")" \
  "input_f6_ic_sha256=$(digest "$out/input-f6_ic.json")" \
  "workload_id=$(jq -r .workload_id "$out/workload.json")" \
  "inherited_f4=$(jq -r .candidate_id "$out/inherited_f4-candidate.json")" \
  "f6_ic=$(jq -r .candidate_id "$out/f6_ic-candidate.json")" \
  > "$out/freeze.txt"
