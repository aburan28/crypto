#!/bin/zsh
# Frozen order-1201, seven-arm, three-round, four-prefix native execution.
# Invoke from the checkout root after RUN_SPEC_V2.json has been committed.
set -u

study=research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/order1000
spec=$study/RUN_SPEC_V2.json
stem=$study/BATCH_1237_1201_V2
binary=target/order-batch-check/release/examples/different_j_trend_batch
checker=target/order-batch-check/release/examples/different_j_trend_batch_check
meter=target/order-batch-check/release/examples/different_j_order_process_meter
raw=$study/raw-v2

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
verify_hash "$checker" '.checker_binary_sha256'
verify_hash "$meter" '.meter_binary_sha256'
verify_hash "$study/run_order1201_v2.zsh" '.script_sha256'
verify_hash "$stem.public.json" '.public_sha256'
verify_hash "$stem.checker.json" '.checker_input_sha256'
verify_hash "$study/FIBER_1237_1201_V1.json" '.fiber_sha256'
if [[ -e "$stem.run.jsonl" || -e "$stem.meter.jsonl" || -e "$raw" ]]; then
  print -u2 -- "V2 output already exists; preserve it and use a new run version"
  exit 98
fi

mkdir -p "$raw"
: > "$stem.run.jsonl"
: > "$stem.meter.jsonl"
typeset -a arms
arms=(source isogeny-image destination-direct full-image source-bsgs rho-negation rho-negation-aa)
any_failed=0
for prefix in 1 4 16 64; do
  for round in 0 1 2; do
    for offset in 0 1 2 3 4 5 6; do
      arm_index=$(( (round + offset) % 7 + 1 ))
      arm=${arms[$arm_index]}
      child="$raw/prefix-$prefix-round-$round-$arm"
      gtimeout 300 "$meter" "$binary" "$stem.public.json" "$round" "$arm" "$prefix" \
        "$child.stdout" "$child.stderr" "$child.meter.json" > "$child.meter.stdout" 2> "$child.meter.stderr"
      child_exit=$?
      print -r -- "$child_exit" > "$child.exit"
      if (( child_exit != 0 )); then any_failed=1; fi
      if [[ -e "$child.stdout" ]]; then cat "$child.stdout" >> "$stem.run.jsonl"; fi
      if [[ -e "$child.meter.json" ]]; then
        cat "$child.meter.json" >> "$stem.meter.jsonl"
        print -- >> "$stem.meter.jsonl"
      fi
    done
  done
done
if (( any_failed == 0 )); then
  gtimeout 300 "$checker" "$stem.public.json" "$stem.checker.json" \
    "$stem.run.jsonl" "$stem.meter.jsonl" > "$stem.replay.json" 2> "$stem.replay.stderr"
  replay_exit=$?
  print -r -- "$replay_exit" > "$stem.replay.exit"
  if (( replay_exit != 0 )); then any_failed=1; fi
else
  print -r -- "skipped_after_failed_children" > "$stem.replay.status"
fi
print -r -- "$any_failed" > "$study/BATCH_RUN_V2.exit"
exit "$any_failed"
