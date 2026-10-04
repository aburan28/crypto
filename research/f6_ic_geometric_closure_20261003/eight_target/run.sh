#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
out=research/f6_ic_geometric_closure_20261003/eight_target
base=research/f6_ic_geometric_closure_20261003/three_way
binary=${IC_WORKER_BINARY:-/private/tmp/f6-ic-target/release/examples/ic_tournament_worker}
digest() { shasum -a 256 "$1" | awk '{print $1}'; }
frozen() { sed -n "s/^$1=//p" "$out/freeze.txt"; }
test "$(digest "$binary")" = "$(frozen binary_sha256)"
for arm in inherited_f4 f5 f6_ic; do
  case "$arm" in
    inherited_f4) field=f4 ;;
    f5) field=f5 ;;
    f6_ic) field=f6 ;;
  esac
  test "$(digest "$base/$arm-candidate.json")" = "$(frozen frozen_candidate_$field)"
done
for role in ic-worker boolean-groebner ic-engine exclusive-measurement prepared-log-state; do
  case "$role" in
    ic-worker) source=examples/ic_tournament_worker.rs ;;
    boolean-groebner) source=src/cryptanalysis/koblitz_groebner.rs ;;
    ic-engine) source=src/cryptanalysis/koblitz_index_calculus.rs ;;
    exclusive-measurement) source=src/cryptanalysis/ic_measurement.rs ;;
    prepared-log-state) source=research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/mathematics.json ;;
  esac
  expected=$(jq -r --arg role "$role" '.record.implementation.components[]|select(.role==$role)|.sha256' "$base/f6_ic-candidate.json")
  test "$(digest "$source")" = "$expected"
done
while IFS="$(printf '\t')" read -r index seed workload fixture f4 f5 f6; do
  [ "$index" = index ] && continue
  dir="$out/T$index"
  test "$(jq -r .workload_id "$dir/workload.json")" = "$workload"
  test "$(jq -r .fixture.target_seeds[0] "$dir/fixture.json")" = "$seed"
  test "$(digest "$dir/fixture.json")" = "$fixture"
  test "$(digest "$dir/input-inherited_f4.json")" = "$f4"
  test "$(digest "$dir/input-f5.json")" = "$f5"
  test "$(digest "$dir/input-f6_ic.json")" = "$f6"
done < "$out/freeze.tsv"
test ! -e "$out/runs"
mkdir "$out/runs"
uname -srm > "$out/runs/host.txt"
rustc --version > "$out/runs/rustc.txt"
for index in 1 2 3 4 5 6 7 8; do
  if [ $((index % 2)) -eq 1 ]; then
    schedule='inherited_f4:1 f6_ic:1 f5:1 f6_ic:2 inherited_f4:2'
  else
    schedule='f6_ic:1 inherited_f4:1 f5:1 inherited_f4:2 f6_ic:2'
  fi
  for item in $schedule; do
    arm=${item%:*}
    repeat=${item#*:}
    stem="$out/runs/T${index}-R${repeat}-${arm}"
    date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.start-utc"
    if RAYON_NUM_THREADS=1 IC_ARTIFACT_CACHE=off IC_F2_BACKEND=cpu \
      gtimeout 180 "$binary" < "$out/T$index/input-$arm.json" \
        > "$stem.stdout.json" 2> "$stem.stderr.txt"; then
      rc=0
    else
      rc=$?
    fi
    printf '%s\n' "$rc" > "$stem.exit-status"
    date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.end-utc"
  done
done
