#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
out=research/f6_ic_geometric_closure_20261003/small_cold
binary=${IC_WORKER_BINARY:-/private/tmp/f6-ic-target/release/examples/ic_tournament_worker}
digest() { shasum -a 256 "$1" | awk '{print $1}'; }
frozen() { sed -n "s/^$1=//p" "$out/freeze.txt"; }
test "$(digest "$binary")" = "$(frozen binary_sha256)"
test "$(digest examples/ic_f6_small_inventory.rs)" = "$(frozen inventory_helper_sha256)"
test "$(digest "$out/fixture.json")" = "$(frozen fixture_sha256)"
test "$(digest "$out/inventory.json")" = "$(frozen inventory_sha256)"
test "$(digest "$out/usable-inventory.json")" = "$(frozen usable_inventory_sha256)"
test "$(jq -r .workload_id "$out/workload.json")" = "$(frozen workload_id)"
for arm in inherited_f4 f6_ic f5; do
  test "$(digest "$out/input-$arm.json")" = "$(frozen input_${arm}_sha256)"
  test "$(jq -r .candidate_id "$out/$arm-candidate.json")" = "$(frozen "$arm")"
  got=$(jq -jcS .record "$out/$arm-candidate.json" | shasum -a 256 | awk '{print $1}')
  test "$got" = "$(jq -r .record_sha256 "$out/$arm-candidate.json")"
done
for role in ic-worker boolean-groebner ic-engine exclusive-measurement; do
  case "$role" in
    ic-worker) source=examples/ic_tournament_worker.rs ;;
    boolean-groebner) source=src/cryptanalysis/koblitz_groebner.rs ;;
    ic-engine) source=src/cryptanalysis/koblitz_index_calculus.rs ;;
    exclusive-measurement) source=src/cryptanalysis/ic_measurement.rs ;;
  esac
  expected=$(jq -r --arg role "$role" '.record.implementation.components[]|select(.role==$role)|.sha256' "$out/f6_ic-candidate.json")
  test "$(digest "$source")" = "$expected"
done
test ! -e "$out/runs"
mkdir "$out/runs"
uname -srm > "$out/runs/host.txt"
rustc --version > "$out/runs/rustc.txt"
for arm in inherited_f4 f6_ic f5; do
  stem="$out/runs/R1-$arm"
  date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.start-utc"
  if RAYON_NUM_THREADS=1 IC_ARTIFACT_CACHE=off IC_F2_BACKEND=cpu \
    gtimeout 180 "$binary" < "$out/input-$arm.json" \
      > "$stem.stdout.json" 2> "$stem.stderr.txt"; then
    rc=0
  else
    rc=$?
  fi
  printf '%s\n' "$rc" > "$stem.exit-status"
  date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.end-utc"
done
