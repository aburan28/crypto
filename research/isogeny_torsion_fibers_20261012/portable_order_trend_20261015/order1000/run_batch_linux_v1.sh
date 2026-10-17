#!/usr/bin/env bash
# Frozen selected order-1201 second-host runner. Invoke from the repository root after building.
set -u
study=research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/order1000
spec="$study/RUN_SPEC_V2.json"
host_spec="$study/HOST_SPEC_V1.json"
out="${1:?usage: run_batch_linux_v1.sh OUTPUT_DIRECTORY}"
binary=target/release/examples/different_j_trend_batch
checker=target/release/examples/different_j_trend_batch_check
meter=target/release/examples/different_j_order_process_meter
mkdir -p "$out"

verify_hash() {
  local actual expected
  actual=$(sha256sum "$1" | cut -d ' ' -f 1)
  expected=$(jq -r "$2" "${3:-$spec}")
  if [[ "$actual" != "$expected" ]]; then
    printf 'hash mismatch: %s: %s versus %s\n' "$1" "$actual" "$expected" >&2
    exit 99
  fi
}

verify_hash examples/different_j_trend_batch.rs '.runner_source_sha256'
verify_hash examples/different_j_trend_batch_check.rs '.checker_source_sha256'
verify_hash examples/different_j_order_process_meter.rs '.meter_source_sha256' "$host_spec"
verify_hash "$spec" '.run_spec_sha256' "$host_spec"
verify_hash .github/workflows/isogeny-fiber-order1201-replay.yml '.workflow_sha256' "$host_spec"
verify_hash "$study/run_batch_linux_v1.sh" '.linux_script_sha256' "$host_spec"
stem="$study/BATCH_1237_1201_V2"
verify_hash "$stem.public.json" '.public_sha256'
verify_hash "$stem.checker.json" '.checker_input_sha256'
verify_hash "$study/FIBER_1237_1201_V1.json" '.fiber_sha256'

(
  set -e
  git rev-parse HEAD
  rustc --version
  cargo --version
  uname -a
  lscpu
  nproc
  cat /proc/self/status
) > "$out/host.txt" || exit 98
sha256sum "$binary" "$checker" "$meter" examples/different_j_trend_batch.rs \
  examples/different_j_trend_batch_check.rs \
  examples/different_j_order_process_meter.rs "$spec" "$host_spec" "$study/run_batch_linux_v1.sh" \
  > "$out/build-sha256.txt" || exit 98
cat /proc/stat > "$out/proc-stat-before.txt"

arms=(source isogeny-image destination-direct full-image source-bsgs rho-negation rho-negation-aa)
any_failed=0
for cell in 1237:1201; do
  p=${cell%%:*}
  r=${cell#*:}
  stem="$study/BATCH_${p}_${r}_V2"
  cell_out="$out/$p-$r"
  mkdir -p "$cell_out/raw"
  : > "$cell_out/run.jsonl"
  : > "$cell_out/meter.jsonl"
  cell_failed=0
  for prefix in 1 4 16 64; do
    for round in 0 1 2; do
      for offset in 0 1 2 3 4 5 6; do
        arm=${arms[$(((round + offset) % 7))]}
        child="$cell_out/raw/prefix-$prefix-round-$round-$arm"
        timeout 300 "$meter" "$binary" "$stem.public.json" "$round" "$arm" "$prefix" \
          "$child.stdout" "$child.stderr" "$child.meter.json" \
          > "$child.meter.stdout" 2> "$child.meter.stderr"
        child_exit=$?
        printf '%s\n' "$child_exit" > "$child.exit"
        if (( child_exit != 0 )); then cell_failed=1; any_failed=1; fi
        if [[ -e "$child.stdout" ]]; then cat "$child.stdout" >> "$cell_out/run.jsonl"; fi
        if [[ -e "$child.meter.json" ]]; then
          cat "$child.meter.json" >> "$cell_out/meter.jsonl"
          printf '\n' >> "$cell_out/meter.jsonl"
        fi
      done
    done
  done
  if (( cell_failed == 0 )); then
    timeout 300 "$checker" "$stem.public.json" "$stem.checker.json" \
      "$cell_out/run.jsonl" "$cell_out/meter.jsonl" \
      > "$cell_out/replay.json" 2> "$cell_out/replay.stderr"
    replay_exit=$?
    printf '%s\n' "$replay_exit" > "$cell_out/replay.exit"
    if (( replay_exit != 0 )); then any_failed=1; fi
  else
    printf 'skipped_after_failed_children\n' > "$cell_out/replay.status"
  fi
done
cat /proc/stat > "$out/proc-stat-after.txt"
printf '%s\n' "$any_failed" > "$out/exit"
(cd "$out" && set -o pipefail && find . -type f ! -name sha256sums-before-manifest.txt -print0 | sort -z | xargs -0 sha256sum) \
  > "$out/sha256sums-before-manifest.txt" || exit 98
exit "$any_failed"
