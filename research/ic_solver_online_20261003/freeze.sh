#!/bin/sh
set -eu

repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/ic_solver_online_20261003
template=research/ic_candidate_tournament_20260915/goal_20260924/prepared-one-target-controls-v1/f5-candidate.json
state=research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/mathematics.json
binary=${IC_WORKER_BINARY:-/private/tmp/ic-solver-e2e-target/release/examples/ic_tournament_worker}
commit=$(git rev-parse HEAD)
digest() { shasum -a 256 "$1" | awk '{print $1}'; }
worker_sha=$(digest examples/ic_tournament_worker.rs)
groebner_sha=$(digest src/cryptanalysis/koblitz_groebner.rs)
ic_sha=$(digest src/cryptanalysis/koblitz_index_calculus.rs)
measurement_sha=$(digest src/cryptanalysis/ic_measurement.rs)
state_sha=$(digest "$state")
binary_sha=$(digest "$binary")

jq -cnS '{schema_version:1,status:"fixture",fixture:{cofactor:"2",curve_a:1,degree:17,generator:["43693","23339"],group_order:"131174",irreducible:{degree:17,low_terms:[0,3]},lambda:"17184",subgroup_order:"65587",target_scalar_constructed:false,target_seeds:[20261003017],targets:[["36389","93850"]]}}' > "$out/fixture.json"

jq -cS --argjson seed 20261003032 \
  '.record | .algorithm_seed=$seed | .cache_policy="warm-prepared" | .input_law="one-new-public-hash-to-curve-point; fixture-construction-excluded" | .point_was_previously_supplied=false | .question="production-prepared-f5-versus-inherited-f4-online" | .resource_envelope.cpu_workers=1 | .resource_envelope.host_class="physical-macos-arm64" | .resource_envelope.total_wall_limit_seconds=600 | .target_seeds=[20261003017] | .targets=[[36389,93850]]' \
  research/ic_candidate_tournament_20260915/goal_20260924/prepared-one-target-controls-v1/f5-workload.json > "$out/workload-record.json"
workload_sha=$(jq -jcS . "$out/workload-record.json" | shasum -a 256 | awk '{print $1}')
jq -cnS --slurpfile record "$out/workload-record.json" --arg sha "$workload_sha" \
  '{record:$record[0],record_sha256:$sha,workload_id:$sha[0:12]}' > "$out/workload.json"

for solver in f5 inherited_f4; do
  if [ "$solver" = f5 ]; then stage=f5; else stage=f4; fi
  jq -cS --arg solver "$solver" --arg commit "$commit" \
    --arg worker_sha "$worker_sha" --arg groebner_sha "$groebner_sha" \
    --arg ic_sha "$ic_sha" --arg measurement_sha "$measurement_sha" \
    --arg binary_sha "$binary_sha" --arg state_sha "$state_sha" \
    '.record |
      .point_decomposition.solver=$solver |
      .point_decomposition.source_sha256=$groebner_sha |
      .point_decomposition.limits.node_budget=8192 |
      .relation_collection.source_sha256=$ic_sha |
      .relation_collection.query_rule="import-certified-ordinary-rows; no new ordinary queries" |
      .relation_linear_algebra.source_sha256=$ic_sha |
      .target_descent.source_sha256=$ic_sha |
      .target_descent.stop_rule.max_trials=32 |
      .implementation={code_snapshot_commit:$commit,
        components:[{role:"ic-worker",sha256:$worker_sha},{role:"boolean-groebner",sha256:$groebner_sha},{role:"ic-engine",sha256:$ic_sha},{role:"exclusive-measurement",sha256:$measurement_sha},{role:"prepared-log-state",sha256:$state_sha}],
        flags:{binary_sha256:$binary_sha,rayon_threads:1,exclusive_phases:true,preparation_mode:"imported-certified-log-table-v1",reusable_template:"fixed-zero-symbolic-specialization-before-online"}}' \
    "$template" > "$out/$solver-record.json"
  sha=$(jq -jcS . "$out/$solver-record.json" | shasum -a 256 | awk '{print $1}')
  id="IC1N17Ckb1fb62PDP3${stage}RCsampleLAgaussTDpdpISO0h$(printf '%.12s' "$sha")"
  jq -cnS --slurpfile record "$out/$solver-record.json" --arg sha "$sha" --arg id "$id" \
    '{record:$record[0],record_sha256:$sha,candidate_id:$id}' > "$out/$solver-candidate.json"
done

printf '%s\n' "source_commit=$commit" "binary_sha256=$binary_sha" "workload_id=$(jq -r .workload_id "$out/workload.json")" \
  "f5=$(jq -r .candidate_id "$out/f5-candidate.json")" \
  "inherited_f4=$(jq -r .candidate_id "$out/inherited_f4-candidate.json")" > "$out/freeze.txt"
