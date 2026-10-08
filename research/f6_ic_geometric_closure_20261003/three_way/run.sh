#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
out=research/f6_ic_geometric_closure_20261003/three_way
binary=${IC_WORKER_BINARY:-/private/tmp/f6-ic-target/release/examples/ic_tournament_worker}
digest() { shasum -a 256 "$1" | awk '{print $1}'; }
frozen() { sed -n "s/^$1=//p" "$out/freeze.txt"; }
test "$(digest "$binary")" = "$(frozen binary_sha256)"
test "$(digest "$out/fixture.json")" = "$(frozen fixture_sha256)"
for arm in inherited_f4 f5 f6_ic; do
  test "$(digest "$out/input-$arm.json")" = "$(frozen input_${arm}_sha256)"
done
for role in ic-worker boolean-groebner ic-engine exclusive-measurement prepared-log-state; do
  case "$role" in
    ic-worker) source=examples/ic_tournament_worker.rs ;;
    boolean-groebner) source=src/cryptanalysis/koblitz_groebner.rs ;;
    ic-engine) source=src/cryptanalysis/koblitz_index_calculus.rs ;;
    exclusive-measurement) source=src/cryptanalysis/ic_measurement.rs ;;
    prepared-log-state) source=research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/mathematics.json ;;
  esac
  expected=$(jq -r --arg role "$role" '.record.implementation.components[] | select(.role==$role) | .sha256' "$out/f6_ic-candidate.json")
  test "$(digest "$source")" = "$expected"
done
test ! -e "$out/runs"
mkdir "$out/runs"
uname -srm > "$out/runs/host.txt"
rustc --version > "$out/runs/rustc.txt"
printf '%s\n' '1 inherited_f4' '1 f5' '1 f6_ic' \
  '2 f6_ic' '2 inherited_f4' '2 f5' \
  '3 f5' '3 f6_ic' '3 inherited_f4' |
while read -r repeat arm; do
  stem="$out/runs/R${repeat}-${arm}"
  date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.start-utc"
  if RAYON_NUM_THREADS=1 IC_ARTIFACT_CACHE=off IC_F2_BACKEND=cpu \
    gtimeout 600 "$binary" < "$out/input-$arm.json" \
      > "$stem.stdout.json" 2> "$stem.stderr.txt"; then
    rc=0
  else
    rc=$?
  fi
  printf '%s\n' "$rc" > "$stem.exit-status"
  date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.end-utc"
done
