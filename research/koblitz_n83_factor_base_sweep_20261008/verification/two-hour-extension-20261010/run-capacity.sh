#!/bin/sh
# Thin, bounded Docker orchestration of the native Rust construction worker.
set -eu
[ "$#" -eq 5 ] || { echo 'usage: run-capacity.sh NEW_OUTPUT BINARY PANEL WALL_SECONDS IMAGE_ID' >&2; exit 2; }
out=$1; binary=$2; panel=$3; wall=$4; image=$5
mkdir "$out"
name="n83-extension-capacity-$$"
source_commit=$(git rev-parse HEAD)
binary_hash=$(shasum -a 256 "$binary" | cut -d ' ' -f1)
panel_hash=$(shasum -a 256 "$panel/manifest.json" | cut -d ' ' -f1)
jq -n --arg source_commit "$source_commit" --arg binary "$binary" --arg binary_sha256 "$binary_hash" --arg panel "$panel" --arg panel_manifest_sha256 "$panel_hash" --arg image "$image" --arg name "$name" --argjson wall "$wall" '{schema:"n83.native-capacity-config/v1",source_commit:$source_commit,binary:$binary,binary_sha256:$binary_sha256,panel:$panel,panel_manifest_sha256:$panel_manifest_sha256,container_image_id:$image,container_name:$name,wall_cap_seconds:$wall,cpu_limit:1,cpuset:"0",memory_limit_bytes:4294967296,swap_limit_bytes:0,orbit_columns:64,solver_search_executed:false,selected_best_total_runtime:null}' > "$out/config.json"
started=$(date +%s)
cleanup() { docker kill "$name" >/dev/null 2>&1 || true; docker rm -f "$name" >/dev/null 2>&1 || true; }
trap cleanup EXIT HUP INT TERM
if docker run --detach --pull never --name "$name" --network none --cpus 1 --cpuset-cpus 0 --memory 4096m --memory-swap 4096m --pids-limit 64 --cap-drop ALL --security-opt no-new-privileges --read-only --tmpfs /tmp:rw,noexec,nosuid,size=64m --user "$(id -u):$(id -g)" --mount "type=bind,src=$panel,dst=/panel,readonly" --mount "type=bind,src=$binary,dst=/worker,readonly" --mount "type=bind,src=$out,dst=/out" "$image" /worker primary-sat-build /panel 64 4096 /out/worker.json > "$out/container-id.txt" 2> "$out/launch.stderr"; then launch=0; else launch=$?; fi
wait_status=0; worker_exit=-1; timed_out=false; cleanup_ok=false
if [ "$launch" -eq 0 ]; then
  if gtimeout "$wall" docker wait "$name" > "$out/exit.txt" 2> "$out/wait.stderr"; then wait_status=0; else wait_status=$?; fi
  if [ "$wait_status" -eq 124 ]; then
    timed_out=true
    docker kill "$name" > "$out/kill.stdout" 2> "$out/kill.stderr" || true
    gtimeout 20 docker wait "$name" > "$out/exit.txt" 2>> "$out/wait.stderr" || true
  fi
  docker logs "$name" > "$out/stdout.log" 2> "$out/stderr.log" || true
  docker inspect --format '{{json .State}}' "$name" > "$out/container-state.json" 2> "$out/inspect.stderr" || true
  if [ -s "$out/exit.txt" ]; then worker_exit=$(tail -1 "$out/exit.txt"); fi
  docker rm -f "$name" > "$out/remove.stdout" 2> "$out/remove.stderr" || true
fi
remaining=$(docker ps -aq --filter "name=^${name}$")
[ -z "$remaining" ] && cleanup_ok=true
ended=$(date +%s)
status=PRODUCER_FAILURE
if [ "$cleanup_ok" != true ]; then status=PRODUCER_FAILURE_cleanup;
elif [ "$timed_out" = true ]; then status=UNKNOWN_wall_cap;
elif [ "$launch" -ne 0 ]; then status=PRODUCER_FAILURE_launch;
elif [ "$worker_exit" = 137 ] || [ "$worker_exit" = 134 ]; then status=UNKNOWN_resource_or_worker_exit;
elif [ "$worker_exit" = 0 ] && jq -e '.schema=="n83.factored-s4-capacity-worker/v1" and .status=="PASS_model_construction_only" and .orbit_columns==64 and .memory_cgroup_limit_bytes==4294967296 and .memory_cgroup_swap_limit_bytes==0 and .sat_variables>0 and .sat_clauses>0 and .solver_search_executed==false and .model_lifting_executed==false and .total_index_calculus_runtime_ms==null and .selected_best_total_runtime==null' "$out/worker.json" >/dev/null 2> "$out/receipt-check.stderr"; then status=PASS_model_construction_only; fi
jq -n --arg status "$status" --argjson launched "$launch" --argjson wait_status "$wait_status" --argjson worker_exit "$worker_exit" --argjson timed_out "$timed_out" --argjson cleanup_ok "$cleanup_ok" --argjson charge "$((ended-started+1))" '{schema:"n83.native-capacity-outer/v1",status:$status,launch_status:$launched,wait_status:$wait_status,worker_exit_code:$worker_exit,timed_out:$timed_out,cleanup_verified:$cleanup_ok,conservative_charged_seconds:$charge,measurement_scope:"construction only; shared host; no comparative runtime",solver_search_executed:false,total_index_calculus_runtime_ms:null,selected_best_total_runtime:null}' > "$out/outer.json"
cat "$out/outer.json"
[ "$status" = PASS_model_construction_only ]
