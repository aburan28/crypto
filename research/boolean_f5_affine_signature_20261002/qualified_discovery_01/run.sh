#!/usr/bin/env bash
# Thin native structural campaign orchestration.
set -euo pipefail

phase=${1:?discovery or holdout}
destination=${2:?new output directory}
discovery=${3:-}
case "$phase" in discovery|holdout) ;; *) echo 'expected discovery or holdout' >&2; exit 2 ;; esac
if [[ -e "$destination" ]]; then
  echo 'refuse to overwrite an attempt directory' >&2
  exit 2
fi
study_dir=$(cd "$(dirname "$0")" && pwd)
repo_root=$(cd "$study_dir/../.." && pwd)
mkdir -p -- "$destination"
out=$(cd "$destination" && pwd)
manifest="$study_dir/Cargo.toml"
isolation="$repo_root/tools/isolated_bench.py"
cp "$study_dir"/{worker.rs,verify.rs,Cargo.toml,Cargo.lock,protocol.json,PROTOCOL.md,README.md,run.sh} "$out/"
git -C "$repo_root" rev-parse HEAD > "$out/git-head.txt"
rustc --version --verbose > "$out/rustc.txt"
uname -a > "$out/host.txt"
: > "$out/raw.jsonl"

python3 "$isolation" busy -- cargo test --offline --release --locked --manifest-path "$manifest" -- --test-threads=1 > "$out/test.stdout" 2> "$out/test.stderr"
python3 "$isolation" busy -- cargo build --offline --release --locked --manifest-path "$manifest" > "$out/build.stdout" 2> "$out/build.stderr"
cp "$study_dir/target/release/boolean-f5-affine-signature" "$out/boolean-f5-affine-signature"
worker="$out/boolean-f5-affine-signature"

if [[ "$phase" == holdout ]]; then
  if [[ -z "$discovery" ]]; then
    "$worker" --failure "$out" missing-discovery
    "$worker" --seal "$out"
    exit 2
  fi
  if ! "$worker" --check-discovery "$discovery" "$out/binding.json" > "$out/binding.stdout" 2> "$out/binding.stderr"; then
    "$worker" --failure "$out" discovery-binding-failed
    "$worker" --seal "$out"
    exit 1
  fi
fi

set +e
python3 "$isolation" busy -- "$worker" --campaign "$phase" "$out/protocol.json" > "$out/raw.jsonl" 2> "$out/worker.stderr"
worker_status=$?
set -e
if [[ $worker_status -ne 0 ]]; then
  "$worker" --failure "$out" worker-failed-or-censored
  "$worker" --seal "$out"
  exit "$worker_status"
fi

set +e
"$worker" --verify "$phase" "$out/raw.jsonl" "$out/results.json" > "$out/verify.stdout" 2> "$out/verify.stderr"
verify_status=$?
set -e
if [[ $verify_status -ne 0 ]]; then
  "$worker" --failure "$out" verifier-failed
  "$worker" --seal "$out"
  exit "$verify_status"
fi
"$worker" --seal "$out"
set +e
"$worker" --verify-bundle "$out"
postseal_status=$?
set -e
if [[ $postseal_status -ne 0 ]]; then
  mv "$out/manifest.json" "$out/manifest-before-postcheck.json"
  "$worker" --failure "$out" post-seal-replay-failed
  "$worker" --seal "$out"
  "$worker" --verify-bundle "$out"
  exit "$postseal_status"
fi
