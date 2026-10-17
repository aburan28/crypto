#!/bin/zsh
# Frozen native seven-arm/three-cell execution. Invoke from checkout root.
set -u

study=research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015
spec=$study/BATCH_RUN_SPEC_V1.json
binary=target/order-batch-check/release/examples/different_j_trend_batch
checker_binary=target/order-batch-check/release/examples/different_j_trend_batch_check
meter_binary=target/order-batch-check/release/examples/different_j_order_process_meter
raw=$study/raw-batch-v1

verify_hash() {
  local actual expected
  actual=$(shasum -a 256 "$1" | cut -d ' ' -f 1)
  expected=$(jq -r "$2" "$spec")
  if [[ "$actual" != "$expected" ]]; then
    print -u2 -- "hash mismatch for $1: $actual versus $expected"
    exit 99
  fi
}

verify_hash examples/different_j_trend_batch.rs '.runner_source_sha256'
verify_hash "$binary" '.runner_binary_sha256'
verify_hash examples/different_j_trend_batch_check.rs '.checker_source_sha256'
verify_hash "$checker_binary" '.checker_binary_sha256'
verify_hash "$meter_binary" '.meter_binary_sha256'
verify_hash "$study/run_batch_v1.zsh" '.script_sha256'
for cell in 409:433 421:457 457:421; do
  cell_p=${cell%%:*}
  cell_r=${cell#*:}
  stem="$study/BATCH_${cell_p}_${cell_r}_V1"
  verify_hash "$stem.public.json" ".public_sha256.\"$cell_p:$cell_r\""
  verify_hash "$stem.checker.json" ".checker_sha256.\"$cell_p:$cell_r\""
  verify_hash "$study/FIBER_${cell_p}_${cell_r}_V1.json" ".fiber_sha256.\"$cell_p:$cell_r\""
  if [[ -e "$stem.run.jsonl" || -e "$stem.meter.jsonl" || -e "$raw/$cell_p-$cell_r" ]]; then
    print -u2 -- "V1 output exists for $cell; preserve it and use a new run version"
    exit 98
  fi
done

mkdir -p "$raw"
typeset -a arms
arms=(source isogeny-image destination-direct full-image source-bsgs rho-negation rho-negation-aa)
any_failed=0
for cell in 409:433 421:457 457:421; do
  cell_p=${cell%%:*}
  cell_r=${cell#*:}
  stem="$study/BATCH_${cell_p}_${cell_r}_V1"
  cell_raw="$raw/$cell_p-$cell_r"
  mkdir -p "$cell_raw"
  : > "$stem.run.jsonl"
  : > "$stem.meter.jsonl"
  cell_failed=0
  for prefix in 1 4 16 64; do
    for round in 0 1 2; do
      for offset in 0 1 2 3 4 5 6; do
        arm_index=$(( (round + offset) % 7 + 1 ))
        arm=${arms[$arm_index]}
        child="$cell_raw/prefix-$prefix-round-$round-$arm"
        gtimeout 300 "$meter_binary" "$binary" "$stem.public.json" "$round" "$arm" "$prefix" \
          "$child.stdout" "$child.stderr" "$child.meter.json" > "$child.meter.stdout" 2> "$child.meter.stderr"
        child_exit=$?
        print -r -- "$child_exit" > "$child.exit"
        if (( child_exit != 0 )); then
          cell_failed=1
          any_failed=1
        fi
        if [[ -e "$child.stdout" ]]; then
          cat "$child.stdout" >> "$stem.run.jsonl"
        fi
        if [[ -e "$child.meter.json" ]]; then
          cat "$child.meter.json" >> "$stem.meter.jsonl"
          print -- >> "$stem.meter.jsonl"
        fi
      done
    done
  done
  if (( cell_failed == 0 )); then
    gtimeout 300 "$checker_binary" "$stem.public.json" "$stem.checker.json" \
      "$stem.run.jsonl" "$stem.meter.jsonl" > "$stem.replay.json" 2> "$stem.replay.stderr"
    replay_exit=$?
    print -r -- "$replay_exit" > "$stem.replay.exit"
    if (( replay_exit != 0 )); then any_failed=1; fi
  else
    print -r -- "skipped_after_failed_children" > "$stem.replay.status"
  fi
done
print -r -- "$any_failed" > "$study/BATCH_RUN_V1.exit"
exit "$any_failed"
