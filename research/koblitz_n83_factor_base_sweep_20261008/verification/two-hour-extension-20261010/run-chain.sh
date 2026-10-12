#!/bin/sh
# Thin shell resource guard; field arithmetic and search remain native Rust.
set -eu
[ "$#" -eq 7 ] || exit 2
out=$1; binary=$2; checkout=$3; panel=$4; fixtures=$5; image=$6; stage=$7
case "$stage" in capacity|preflight) wall=60; worker_cap=50; trials=0;; search) wall=600; worker_cap=590; trials=1;; *) exit 2;; esac
cd "$checkout"
[ -z "$(git status --porcelain)" ] || exit 2
commit=$(git rev-parse HEAD)
common=$(git rev-parse --path-format=absolute --git-common-dir)
mkdir "$out"
name="n83-extension-chain-$stage-$$"
shasum -a 256 "$binary" "$panel/manifest.json" "$panel/replay.json" "$panel/upload-receipt.json" "$fixtures/probe-corpus.json" "$fixtures/probe-validation.json" > "$out/frozen-sha256.txt"
jq -n --arg stage "$stage" --arg source_commit "$commit" --arg image "$image" --argjson trials "$trials" --argjson worker_cap "$worker_cap" --argjson wall "$wall" '{schema:"n83.native-chain-config/v1",stage:$stage,source_commit:$source_commit,image_id:$image,columns:64,policy:"public_x_hash",seed:2026100801,summands:5,identity_mask:0,max_variables:150000,max_domain_clauses:3000000,max_models:1,conflict_budget:10000,max_trials:$trials,worker_wall_seconds:$worker_cap,outer_wall_seconds:$wall,memory_limit_bytes:4294967296,swap_limit_bytes:0,cpu_limit:1,cpuset:"0",total_index_calculus_runtime_ms:null,selected_best_total_runtime:null}' > "$out/config.json"
set -- docker run --detach --pull never --name "$name" --network none --cpus 1 --cpuset-cpus 0 --memory 4096m --memory-swap 4096m --pids-limit 64 --cap-drop ALL --security-opt no-new-privileges --read-only --tmpfs /tmp:rw,noexec,nosuid,size=64m --user "$(id -u):$(id -g)" --env GIT_OPTIONAL_LOCKS=0 --env GIT_CONFIG_COUNT=1 --env GIT_CONFIG_KEY_0=safe.directory --env "GIT_CONFIG_VALUE_0=$checkout" --mount "type=bind,src=$checkout,dst=$checkout,readonly" --mount "type=bind,src=$common,dst=$common,readonly" --mount "type=bind,src=$binary,dst=/worker,readonly" --mount "type=bind,src=$panel,dst=/panel,readonly" --mount "type=bind,src=$fixtures,dst=/fixtures,readonly" --mount "type=bind,src=$out,dst=/out" --workdir "$checkout" "$image" /worker
if [ "$stage" = capacity ]; then
  set -- "$@" primary-chain-build /panel 64 public_x_hash 2026100801 5 0 150000 3000000 4096 /out/worker.json
else
  set -- "$@" primary-chain-cold /panel /fixtures 64 public_x_hash 2026100801 5 "$trials" 150000 3000000 1 10000 "$worker_cap" 4096 /out/run
fi
jq -n --args '$ARGS.positional' -- "$@" > "$out/command.json"
cleanup() { docker kill "$name" >/dev/null 2>&1 || true; docker rm -f "$name" >/dev/null 2>&1 || true; }
trap cleanup EXIT HUP INT TERM
started=$(date +%s)
if "$@" > "$out/container-id.txt" 2> "$out/launch.stderr"; then launch=0; else launch=$?; fi
wait_rc=0; worker_rc=-1; timed_out=false; cleanup_ok=false
if [ "$launch" -eq 0 ]; then
  if gtimeout "$wall" docker wait "$name" > "$out/exit.txt" 2> "$out/wait.stderr"; then wait_rc=0; else wait_rc=$?; fi
  if [ "$wait_rc" -eq 124 ]; then
    timed_out=true
    docker kill "$name" > "$out/kill.stdout" 2> "$out/kill.stderr" || true
    gtimeout 20 docker wait "$name" > "$out/exit.txt" 2>> "$out/wait.stderr" || true
  fi
  docker logs "$name" > "$out/stdout.log" 2> "$out/stderr.log" || true
  docker inspect --format '{{json .State}}' "$name" > "$out/container-state.json" 2> "$out/inspect.stderr" || true
  if [ -s "$out/exit.txt" ]; then worker_rc=$(tail -1 "$out/exit.txt"); fi
  docker rm -f "$name" > "$out/remove.stdout" 2> "$out/remove.stderr" || true
fi
[ -z "$(docker ps -aq --filter "name=^${name}$")" ] && cleanup_ok=true
ended=$(date +%s)
result=PRODUCER_FAILURE
if [ "$cleanup_ok" != true ]; then result=PRODUCER_FAILURE_cleanup;
elif [ "$timed_out" = true ]; then result=UNKNOWN_wall_cap;
elif [ "$launch" -ne 0 ]; then result=PRODUCER_FAILURE_launch;
elif [ "$worker_rc" = 124 ]; then result=UNKNOWN_worker_cap;
elif [ "$worker_rc" = 137 ] || [ "$worker_rc" = 134 ]; then result=UNKNOWN_resource_or_worker_exit;
elif [ "$worker_rc" = 0 ]; then
  if [ "$stage" = capacity ]; then
    if jq -e --arg commit "$commit" '.schema=="n83.chain-s3-capacity-worker/v1" and .source_commit==$commit and .orbit_columns==64 and .summands==5 and .identity_mask==0 and .max_variables==150000 and .max_domain_clauses==3000000 and .memory_cgroup_limit_bytes==4294967296 and .memory_cgroup_swap_limit_bytes==0 and .solver_search_executed==false and .relation_stage_executed==false and .rank_stage_executed==false and .total_index_calculus_runtime_ms==null and .selected_best_total_runtime==null and (.status=="PASS_model_construction_only" or .status=="UNKNOWN_variable_cap" or .status=="UNKNOWN_domain_clause_cap") and (.status!="PASS_model_construction_only" or (.sat_variables>0 and .sat_clauses>0))' "$out/worker.json" > "$out/receipt-check.stdout" 2> "$out/receipt-check.stderr"; then result=$(jq -r .status "$out/worker.json"); fi
  fi
  if [ "$stage" != capacity ] && jq -e --arg commit "$commit" --argjson trials "$trials" --argjson cap "$worker_cap" '.schema=="n83.primary-cold-result/v1" and .source_commit==$commit and .orbit_columns==64 and .summands==5 and .strategy=="chain-s3" and .policy=="public_x_hash" and .seed==2026100801 and .chain_s3_limits=={max_variables:150000,max_domain_clauses:3000000,max_models:1,conflict_budget:10000} and .source_attestation_mode=="clean_checkout" and .point_set_blake3=="6dee0500cd0c94191903ecb4558fe7e91f32841b95fadf66aa9f25f4000879f1" and .max_trials==$trials and .budget_seconds==$cap and .memory_cgroup_limit_bytes==4294967296 and .memory_cgroup_swap_limit_bytes==0 and .column_log_verification==false and .total_index_calculus_runtime_ms==null and .selected_best_total_runtime==null and .report.trials<=$trials and .report.orbit_count==64 and .report.factor_base_size==10624 and .report.relations==(.report.independent_relations+.report.dependent_relations) and .report.inconsistent_relations==0 and .report.verification_failures==0 and .report.sat_invalid_models==0 and .solver_stage_executed==(.report.sat_calls>0) and (.status=="PREFLIGHT_ONLY" or .status=="UNKNOWN_solver_cap" or .status=="UNKNOWN_trial_cap" or .status=="INADMISSIBLE_cofactor_class" or .status=="PASS_verified_target_only") and (.status!="PREFLIGHT_ONLY" or ($trials==0 and .report.trials==0 and .report.sat_calls==0)) and (.status!="PASS_verified_target_only" or .verified_log!=null)' "$out/run/summary.json" > "$out/receipt-check.stdout" 2> "$out/receipt-check.stderr"; then result=$(jq -r .status "$out/run/summary.json"); fi
fi
if [ "$(git rev-parse HEAD)" != "$commit" ] || [ -n "$(git status --porcelain)" ] || ! shasum -a 256 -c "$out/frozen-sha256.txt" > "$out/frozen-check.log" 2>&1; then result=PRODUCER_FAILURE_source_or_input_changed; fi
jq -n --arg stage "$stage" --arg result "$result" --argjson launch "$launch" --argjson wait_rc "$wait_rc" --argjson worker_rc "$worker_rc" --argjson timed_out "$timed_out" --argjson cleanup "$cleanup_ok" --argjson charge "$((ended-started+1))" '{schema:"n83.native-chain-outer/v1",stage:$stage,status:$result,launch_status:$launch,wait_status:$wait_rc,worker_exit_code:$worker_rc,timed_out:$timed_out,cleanup_verified:$cleanup,conservative_charged_seconds:$charge,total_index_calculus_runtime_ms:null,selected_best_total_runtime:null}' > "$out/outer.json"
cat "$out/outer.json"
case "$result" in PASS_*|PREFLIGHT_ONLY|UNKNOWN_*|INADMISSIBLE_*) exit 0;; *) exit 1;; esac
