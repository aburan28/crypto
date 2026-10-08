#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/ic_solver_online_20261003
state=research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/mathematics.json
for solver in f5 inherited_f4; do
  jq -cnS --arg solver "$solver" --slurpfile prepared "$state" \
    '{mode:"ic",degree:17,curve_a:1,target_seeds:[],
      public_targets:[["36389","93850"]],algorithm_seed:20261003032,
      factor_base:{kind:"standard_subspace",dimension:6},
      config:{solver:$solver,linear_algebra:"dense",summands:3,
        groebner_degree:3,node_budget:8192,batch_trials:1,max_trials:32,
        rho_parallel_walks:1},exclusive_phases:true,prepared:$prepared[0]}' \
    > "$out/input-$solver.json"
done
jq -cnS '{mode:"rho",degree:17,curve_a:1,target_seeds:[],
  public_targets:[["36389","93850"]],algorithm_seed:20261003032,
  factor_base:{kind:"standard_subspace",dimension:6},
  config:{rho_parallel_walks:1,max_trials:65536,batch_trials:1},
  exclusive_phases:true}' > "$out/input-rho.json"
