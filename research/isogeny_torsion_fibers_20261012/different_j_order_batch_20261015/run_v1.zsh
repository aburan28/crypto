#!/bin/zsh
set -u

repo=/Volumes/SSD990-2/.codex/worktrees/d92d/crypto
study=$repo/research/isogeny_torsion_fibers_20261012/different_j_order_batch_20261015
binary=$repo/target/order-batch-check/release/examples/different_j_order_batch
checker_binary=$repo/target/order-batch-check/release/examples/different_j_order_batch_check
meter_binary=$repo/target/order-batch-check/release/examples/different_j_order_process_meter
public=$study/PUBLIC_V1.json
checker_input=$study/CHECKER_V1.json
spec=$study/RUN_SPEC_V1.json
raw=$study/raw-v1

verify_hash() {
  local actual expected
  actual=$(shasum -a 256 "$1" | cut -d ' ' -f 1)
  expected=$(jq -r "$2" "$spec")
  if [[ "$actual" != "$expected" ]]; then
    print -u2 -- "hash mismatch for $1: $actual versus $expected"
    exit 99
  fi
}

verify_hash "$repo/examples/different_j_order_batch.rs" '.runner_source_sha256'
verify_hash "$binary" '.runner_binary_sha256'
verify_hash "$repo/examples/different_j_order_batch_check.rs" '.checker_source_sha256'
verify_hash "$checker_binary" '.checker_binary_sha256'
verify_hash "$repo/examples/different_j_order_process_meter.rs" '.meter_source_sha256'
verify_hash "$meter_binary" '.meter_binary_sha256'
verify_hash "$public" '.public_sha256'
verify_hash "$checker_input" '.checker_input_sha256'
verify_hash "$study/run_v1.zsh" '.script_sha256'

if [[ -e "$study/CALIBRATION_V1.jsonl" || -e "$study/RUN_V1.jsonl" || -e "$study/METERS_V1.jsonl" || -e "$raw" ]]; then
  print -u2 -- "V1 output exists; preserve it and use a new version"
  exit 98
fi
gtimeout 120 "$binary" calibrate > "$study/CALIBRATION_V1.jsonl" 2> "$study/CALIBRATION_V1.stderr"
calibration_exit=$?
print -r -- "$calibration_exit" > "$study/CALIBRATION_V1.exit"
if (( calibration_exit )); then
  print -u2 -- "calibration failed; raw receipt preserved"
  exit "$calibration_exit"
fi
mkdir -p "$raw"
: > "$study/RUN_V1.jsonl"
: > "$study/METERS_V1.jsonl"
typeset -a arms
arms=(source isogeny-image destination-direct full-image source-bsgs rho-negation rho-negation-aa)
failed=0
for prefix in 1 4 16 64; do
  for round in 0 1 2; do
    for offset in 0 1 2 3 4 5 6; do
      index=$(( (round + offset) % 7 + 1 ))
      arm=${arms[$index]}
      stem=$raw/prefix-$prefix-round-$round-$arm
      if gtimeout 120 "$meter_binary" "$binary" "$public" "$round" "$arm" "$prefix" \
        "$stem.stdout" "$stem.stderr" "$stem.meter.json" > "$stem.meter.stdout" 2> "$stem.meter.stderr"; then
        run_exit_code=0
      else
        run_exit_code=$?
        failed=1
      fi
      print -r -- "$run_exit_code" > "$stem.exit"
      if [[ -e "$stem.stdout" ]]; then
        cat "$stem.stdout" >> "$study/RUN_V1.jsonl"
      fi
      if [[ -e "$stem.meter.json" ]]; then
        cat "$stem.meter.json" >> "$study/METERS_V1.jsonl"
        print -- >> "$study/METERS_V1.jsonl"
      fi
    done
  done
done

if (( failed )); then
  print -u2 -- "one or more arm processes failed; raw receipts preserved"
  exit 1
fi

gtimeout 120 "$checker_binary" "$public" "$checker_input" "$study/RUN_V1.jsonl" "$study/METERS_V1.jsonl" \
  > "$study/REPLAY_V1.json" 2> "$study/REPLAY_V1.stderr"
