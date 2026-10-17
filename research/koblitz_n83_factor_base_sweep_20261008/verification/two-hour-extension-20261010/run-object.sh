#!/bin/sh
# Thin native-worker launch guard. All field arithmetic remains in Rust.
set -eu
[ "$#" -eq 8 ] || { echo 'usage: run-object.sh NEW_STAGE BINARY CHECKOUT IMAGE construct|replay OBJECT_DIR WORKER_CAP OUTER_CAP' >&2; exit 2; }
out=$1; binary=$2; checkout=$3; image=$4; stage=$5; object=$6; worker_cap=$7; wall=$8
[ "$stage" = construct ] || [ "$stage" = replay ] || exit 2
cd "$checkout"
[ -z "$(git status --porcelain)" ] || { echo 'source checkout is dirty' >&2; exit 2; }
commit=$(git rev-parse HEAD)
common=$(git rev-parse --path-format=absolute --git-common-dir)
mkdir "$out"
name="n83-extension-$stage-$$"
binary_hash=$(shasum -a 256 "$binary" | cut -d ' ' -f1)
jq -n --arg stage "$stage" --arg source_commit "$commit" --arg binary_sha256 "$binary_hash" --arg image "$image" --argjson worker_cap "$worker_cap" --argjson wall "$wall" '{schema:"n83.native-v2-stage-config/v1",stage:$stage,source_commit:$source_commit,binary_sha256:$binary_sha256,container_image_id:$image,worker_wall_seconds:$worker_cap,outer_wall_seconds:$wall,cpu_limit:1,cpuset:"0",memory_limit_bytes:4294967296,swap_limit_bytes:0,curve_a:0,policy:"public_x_hash",columns:1182,seed:2026100801,total_index_calculus_runtime_ms:null,selected_best_total_runtime:null}' > "$out/config.json"
set -- docker run --detach --pull never --name "$name" --network none --cpus 1 --cpuset-cpus 0 --memory 4096m --memory-swap 4096m --pids-limit 64 --cap-drop ALL --security-opt no-new-privileges --read-only --tmpfs /tmp:rw,noexec,nosuid,size=64m --user "$(id -u):$(id -g)" --env GIT_OPTIONAL_LOCKS=0 --env GIT_CONFIG_COUNT=1 --env GIT_CONFIG_KEY_0=safe.directory --env "GIT_CONFIG_VALUE_0=$checkout" --mount "type=bind,src=$checkout,dst=$checkout,readonly" --mount "type=bind,src=$common,dst=$common,readonly" --mount "type=bind,src=$binary,dst=/worker,readonly" --mount "type=bind,src=$out,dst=/out" --workdir "$checkout"
if [ "$stage" = construct ]; then
  [ "$object" = "$out/object" ] || exit 2
  set -- "$@" "$image" /worker v2-construct-one /out/object 0 public_x_hash 1182 2026100801 "$worker_cap"
else
  [ -d "$object" ] || exit 2
  set -- "$@" --mount "type=bind,src=$object,dst=/object" "$image" /worker v2-replay-one /object "$worker_cap"
fi
jq -n --args '$ARGS.positional' -- "$@" > "$out/command.json"
cleanup() { docker kill "$name" >/dev/null 2>&1 || true; docker rm -f "$name" >/dev/null 2>&1 || true; }
trap cleanup EXIT HUP INT TERM
started=$(date +%s)
if "$@" > "$out/container-id.txt" 2> "$out/launch.stderr"; then launch=0; else launch=$?; fi
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
[ -z "$(docker ps -aq --filter "name=^${name}$")" ] && cleanup_ok=true
ended=$(date +%s)
status=PRODUCER_FAILURE
if [ "$cleanup_ok" != true ]; then status=PRODUCER_FAILURE_cleanup;
elif [ "$timed_out" = true ]; then status=UNKNOWN_wall_cap;
elif [ "$launch" -ne 0 ]; then status=PRODUCER_FAILURE_launch;
elif [ "$worker_exit" = 124 ]; then status=UNKNOWN_worker_cap;
elif [ "$worker_exit" = 137 ] || [ "$worker_exit" = 134 ]; then status=UNKNOWN_resource_or_worker_exit;
elif [ "$worker_exit" = 0 ]; then
  if [ "$stage" = construct ]; then
    if jq -e --arg commit "$commit" --argjson cap "$worker_cap" '.schema=="n83.factor-base-panel/v2-size-frontier" and .status=="completed_factor_base_object" and .source_commit==$commit and .budget_seconds==$cap and .completed_base_count==1 and .selected_best_total_runtime==null and (.bases|length)==1 and .bases[0].a==0 and .bases[0].columns==1182 and .bases[0].policy=="public_x_hash" and .bases[0].seed==2026100801 and .bases[0].points==196212 and .bases[0].total_index_calculus_runtime_ms==null' "$object/manifest.json" > "$out/receipt-check.stdout" 2> "$out/receipt-check.stderr"; then status=PASS_construction_only; fi
  elif jq -e '.schema=="n83.factor-base-replay/v2-size-frontier" and .status=="PASS"' "$object/replay.json" > "$out/receipt-check.stdout" 2> "$out/receipt-check.stderr"; then status=PASS_generic_replay_only; fi
fi
if [ "$(git rev-parse HEAD)" != "$commit" ] || [ -n "$(git status --porcelain)" ] || [ "$(shasum -a 256 "$binary" | cut -d ' ' -f1)" != "$binary_hash" ]; then status=PRODUCER_FAILURE_source_changed; fi
jq -n --arg stage "$stage" --arg status "$status" --argjson launch "$launch" --argjson wait_status "$wait_status" --argjson worker_exit "$worker_exit" --argjson timed_out "$timed_out" --argjson cleanup "$cleanup_ok" --argjson charge "$((ended-started+1))" '{schema:"n83.native-v2-stage-outer/v1",stage:$stage,status:$status,launch_status:$launch,wait_status:$wait_status,worker_exit_code:$worker_exit,timed_out:$timed_out,cleanup_verified:$cleanup,conservative_charged_seconds:$charge,total_index_calculus_runtime_ms:null,selected_best_total_runtime:null}' > "$out/outer.json"
cat "$out/outer.json"
case "$status" in PASS_*) exit 0;; *) exit 1;; esac
