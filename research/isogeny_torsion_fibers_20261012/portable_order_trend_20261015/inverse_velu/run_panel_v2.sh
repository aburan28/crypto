#!/usr/bin/env bash
# Frozen native V2 inverse-Velu comparison. Run from the repository root.
set -u
study=research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/inverse_velu
old=research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/order1000
spec="$study/RUN_SPEC_V2.json"
out="${1:?usage: run_panel_v2.sh NEW_OUTPUT_DIRECTORY}"
test ! -e "$out" || { echo 'output directory already exists' >&2; exit 97; }

verify_hash() {
  local actual expected
  actual=$(shasum -a 256 "$1" | awk '{print $1}')
  expected=$(jq -r "$2" "$spec")
  if [[ "$actual" != "$expected" ]]; then
    printf 'hash mismatch: %s: %s versus %s\n' "$1" "$actual" "$expected" >&2
    exit 99
  fi
}

binary=target/release/examples/different_j_inverse_velu_oracle
checker=target/release/examples/different_j_inverse_velu_check
meter=target/release/examples/different_j_order_process_meter
public="$old/BATCH_1237_1201_V2.public.json"
answers="$old/BATCH_1237_1201_V2.checker.json"
verify_hash examples/different_j_inverse_velu_oracle.rs '.runner_source_sha256'
verify_hash examples/different_j_inverse_velu_check.rs '.checker_source_sha256'
verify_hash examples/different_j_order_process_meter.rs '.meter_source_sha256'
verify_hash "$binary" '.runner_binary_sha256'
verify_hash "$checker" '.checker_binary_sha256'
verify_hash "$meter" '.meter_binary_sha256'
verify_hash "$public" '.public_sha256'
verify_hash "$answers" '.checker_input_sha256'
verify_hash "$study/run_panel_v2.sh" '.script_sha256'

mkdir -p "$out/raw"
(
  set -e
  git rev-parse HEAD
  date -u '+%Y-%m-%dT%H:%M:%SZ'
  uname -a
  rustc --version
  cargo --version
  sysctl -n machdep.cpu.brand_string
  sysctl -n hw.logicalcpu
  sysctl -n hw.physicalcpu
) > "$out/host.txt"
shasum -a 256 "$binary" "$checker" "$meter" "$spec" "$study/run_panel_v2.sh" > "$out/build-sha256.txt"
: > "$out/run.jsonl"
: > "$out/meter.jsonl"

arms=(inverse-velu full-image source-bsgs rho-negation rho-negation-aa)
any_failed=0
for prefix in 1 4 16 64; do
  for round in 0 1 2; do
    for offset in 0 1 2 3 4; do
      arm=${arms[$(((round + offset) % 5))]}
      child="$out/raw/prefix-$prefix-round-$round-$arm"
      timeout 300 "$meter" "$binary" "$public" "$round" "$arm" "$prefix" \
        "$child.stdout" "$child.stderr" "$child.meter.json" \
        > "$child.meter.stdout" 2> "$child.meter.stderr"
      child_exit=$?
      printf '%s\n' "$child_exit" > "$child.exit"
      if (( child_exit != 0 )); then any_failed=1; fi
      if [[ -e "$child.stdout" ]]; then cat "$child.stdout" >> "$out/run.jsonl"; fi
      if [[ -e "$child.meter.json" ]]; then
        cat "$child.meter.json" >> "$out/meter.jsonl"
        printf '\n' >> "$out/meter.jsonl"
      fi
    done
  done
done
if (( any_failed == 0 )); then
  timeout 300 "$checker" "$public" "$answers" "$out/run.jsonl" "$out/meter.jsonl" \
    > "$out/replay.json" 2> "$out/replay.stderr"
  replay_exit=$?
  printf '%s\n' "$replay_exit" > "$out/replay.exit"
  if (( replay_exit != 0 )); then any_failed=1; fi
else
  printf 'skipped_after_failed_children\n' > "$out/replay.status"
fi
printf '%s\n' "$any_failed" > "$out/exit"
(cd "$out" && find . -type f ! -name ARTIFACT_FILES.sha256 -print0 | sort -z | xargs -0 shasum -a 256) \
  > "$out/ARTIFACT_FILES.sha256"
exit "$any_failed"
