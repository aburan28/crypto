#!/usr/bin/env bash
# Thin native-study build/run orchestration. The producer and verifier are Rust.
set -euo pipefail

phase=${1:?discovery or full}
destination=${2:?new output directory}
discovery=${3:-}
case "$phase" in discovery|full) ;; *) echo 'expected discovery or full' >&2; exit 2 ;; esac
if [[ $(uname -s) != Linux ]]; then
  echo 'qualified runs require Linux CPU affinity and pressure records' >&2
  exit 2
fi
if [[ -e "$destination" ]]; then
  echo 'output directory already exists' >&2
  exit 2
fi
study_dir=$(cd "$(dirname "$0")" && pwd)
repo_root=$(cd "$study_dir/../.." && pwd)
mkdir -p -- "$destination"
out=$(cd "$destination" && pwd)
isolation="$repo_root/tools/isolated_bench.py"
manifest="$study_dir/Cargo.toml"

cp "$study_dir"/{worker.rs,verify.rs,Cargo.toml,Cargo.lock,protocol.json,PROTOCOL.md,README.md,run.sh} "$out/"
cp "$study_dir/baseline_worker.rs" "$out/baseline_worker.rs"
cp /proc/cpuinfo "$out/cpuinfo.txt"
cp /proc/meminfo "$out/meminfo.txt"
rustc --version --verbose > "$out/rustc.txt"
git -C "$repo_root" rev-parse HEAD > "$out/git-head.txt"

python3 "$isolation" busy -- cargo build --release --locked --manifest-path "$manifest" > "$out/build.stdout" 2> "$out/build.stderr"
python3 "$isolation" busy -- cargo test --release --locked --manifest-path "$manifest" -- --test-threads=1 > "$out/test.stdout" 2> "$out/test.stderr"
python3 "$isolation" busy -- rustc --edition=2021 -O "$out/baseline_worker.rs" -o "$out/baseline" > "$out/baseline-build.stdout" 2> "$out/baseline-build.stderr"
python3 "$isolation" busy -- rustc --edition=2021 -O --test "$out/baseline_worker.rs" -o "$out/baseline-tests" > "$out/baseline-test-build.stdout" 2> "$out/baseline-test-build.stderr"
python3 "$isolation" busy -- "$out/baseline-tests" --test-threads=1 > "$out/baseline-test.stdout" 2> "$out/baseline-test.stderr"
cp "$study_dir/target/release/boolean-envelope-indexed" "$out/boolean-envelope-indexed"
worker="$out/boolean-envelope-indexed"
: > "$out/raw.jsonl"
: > "$out/worker.stderr"

if [[ "$phase" == full ]]; then
  if [[ -z "$discovery" ]]; then
    "$worker" --failure "$out" discovery-binding-failed
    "$worker" --seal "$out"
    echo 'full phase requires the sealed qualified discovery' >&2
    exit 2
  fi
  if ! "$worker" --check-discovery "$discovery" "$study_dir" "$out/binding.json" > "$out/binding.stdout" 2> "$out/binding.stderr"; then
    "$worker" --failure "$out" discovery-binding-failed
    "$worker" --seal "$out"
    exit 1
  fi
fi

set +e
"$worker" --wait-quiet "$out/readiness.json" > "$out/readiness.stdout" 2> "$out/readiness.stderr"
readiness_status=$?
set -e
if [[ $readiness_status -ne 0 ]]; then
  "$worker" --failure "$out" resource-run-failed
  "$worker" --seal "$out"
  exit "$readiness_status"
fi

set +e
python3 "$isolation" run --cpus 3 --out "$out/conditions.jsonl" --label "envelope-indexed/$phase" \
  --settle 2 --max-other-cpu 0.10 --max-psi 5.0 -- \
  "$worker" --campaign "$phase" "$out/protocol.json" > "$out/raw.jsonl" 2> "$out/worker.stderr"
resource_status=$?
set -e
if [[ $resource_status -ne 0 ]]; then
  "$worker" --failure "$out" resource-run-failed
  "$worker" --seal "$out"
  exit "$resource_status"
fi

set +e
"$worker" --verify "$phase" "$out/protocol.json" "$out/raw.jsonl" "$out/conditions.jsonl" "$out/results.json" > "$out/verify.stdout" 2> "$out/verify.stderr"
verify_status=$?
set -e
if [[ $verify_status -ne 0 ]]; then
  "$worker" --failure "$out" verifier-failed
  "$worker" --seal "$out"
  exit "$verify_status"
fi
"$worker" --seal "$out"
