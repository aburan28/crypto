#!/usr/bin/env bash
# Thin native orchestration; the existing repository controller owns isolation.
set -euo pipefail

phase=${1:?discovery or holdout}
destination=${2:?new output directory}
discovery=${3:-}
case "$phase" in discovery|holdout) ;; *) echo 'expected discovery or holdout' >&2; exit 2 ;; esac
if [[ $(uname -s) != Linux || $(uname -m) != x86_64 ]]; then
  echo 'qualified phase requires Linux x86-64 AVX2 and PSI receipts' >&2
  exit 2
fi
if [[ -e "$destination" ]]; then
  echo 'refuse to overwrite an attempt directory' >&2
  exit 2
fi
study_dir=$(cd "$(dirname "$0")" && pwd)
repo_root=$(cd "$study_dir/../.." && pwd)
if [[ -n $(git -C "$repo_root" status --porcelain --untracked-files=normal) ]]; then
  echo 'qualified execution requires a committed clean source tree' >&2
  exit 2
fi
mkdir -p -- "$destination"
out=$(cd "$destination" && pwd)
manifest="$study_dir/Cargo.toml"
isolation="$repo_root/tools/isolated_bench.py"

cp "$study_dir"/{worker.rs,campaign.rs,Cargo.toml,Cargo.lock,protocol.json,PROTOCOL.md,README.md,run.sh} "$out/"
cp "$repo_root/src/cryptanalysis/gf2_elim.rs" "$out/gf2_elim.rs"
cp /proc/cpuinfo "$out/cpuinfo.txt"
cp /proc/meminfo "$out/meminfo.txt"
cat /sys/devices/system/cpu/cpu2/topology/thread_siblings_list > "$out/cpu-siblings.txt"
rustc --version --verbose > "$out/rustc.txt"
git -C "$repo_root" rev-parse HEAD > "$out/git-head.txt"
: > "$out/raw.jsonl"

unset KIC_F5_FUSED_BUILD KIC_GF2_SIMD KIC_GF2_PARALLEL_WORDS KIC_GF2_BRANCHLESS_STRIP
unset KIC_GF2_FORCE_AVX2 KIC_GF2_DEFER_ABOVE KIC_DEV_PROFILE KIC_DEV_ASSERT_PREFIX
export KIC_F5_DIRECT_PACK=1
export KIC_F5_UNPACK_DIRECT=1
export KIC_GF2_TABLES=4
export KIC_F5_AVX512_UNPACK=0
export KIC_GF2_REUSE_TABLE=0
export RAYON_NUM_THREADS=1

python3 "$isolation" busy -- cargo build --release --locked --manifest-path "$manifest" > "$out/build.stdout" 2> "$out/build.stderr"
python3 "$isolation" busy -- cargo test --release --locked --manifest-path "$manifest" -- --test-threads=1 > "$out/test.stdout" 2> "$out/test.stderr"
cp "$study_dir/target/release/boolean-f5-support-aware" "$out/boolean-f5-support-aware"
worker="$out/boolean-f5-support-aware"

if [[ "$phase" == holdout ]]; then
  if [[ -z "$discovery" ]] || ! "$worker" --check-discovery "$discovery" "$out/binding.json" > "$out/binding.stdout" 2> "$out/binding.stderr"; then
    "$worker" --failure "$out" discovery-binding-failed
    "$worker" --seal "$out"
    exit 2
  fi
fi

set +e
"$worker" --wait-quiet "$out/readiness.json" > "$out/readiness.stdout" 2> "$out/readiness.stderr"
readiness_status=$?
set -e
if [[ $readiness_status -ne 0 ]]; then
  "$worker" --failure "$out" resource-readiness-failed
  "$worker" --seal "$out"
  exit "$readiness_status"
fi

set +e
python3 "$isolation" run --cpus "$(cat "$out/cpu-siblings.txt")" --out "$out/conditions.jsonl" --label "f5-support-aware/$phase" \
  --settle 2 --max-other-cpu 0.10 --max-psi 5.0 -- \
  timeout --signal=TERM --kill-after=5s 2100 \
  "$worker" --campaign "$phase" "$out/protocol.json" > "$out/raw.jsonl" 2> "$out/worker.stderr"
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
