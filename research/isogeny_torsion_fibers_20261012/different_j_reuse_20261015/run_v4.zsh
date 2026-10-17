#!/bin/zsh
set -u

repo=/Volumes/SSD990-2/.codex/worktrees/d92d/crypto
study=$repo/research/isogeny_torsion_fibers_20261012/different_j_reuse_20261015
binary=$repo/target/release/examples/nonisomorphic_extension_isogeny
source=$repo/examples/nonisomorphic_extension_isogeny.rs
meter_source=$repo/examples/different_j_process_meter.rs
meter_binary=$repo/target/release/examples/different_j_process_meter
public=$study/PUBLIC_BATCH_V2.json
checker=$study/CHECKER_BATCH_V2.json
spec=$study/RUN_SPEC_V4.json
raw=$study/raw-v4

verify_hash() {
  local actual expected
  actual=$(shasum -a 256 "$1" | cut -d ' ' -f 1)
  expected=$(jq -r "$2" "$spec")
  if [[ "$actual" != "$expected" ]]; then
    print -u2 -- "hash mismatch for $1: $actual versus $expected"
    exit 99
  fi
}

verify_hash "$source" '.source_sha256'
verify_hash "$binary" '.binary_sha256'
verify_hash "$meter_source" '.meter_source_sha256'
verify_hash "$meter_binary" '.meter_binary_sha256'
verify_hash "$public" '.public_sha256'
verify_hash "$checker" '.checker_sha256'
verify_hash "$study/run_v4.zsh" '.script_sha256'

if [[ -e "$study/RUN_V4.jsonl" || -e "$study/METERS_V4.jsonl" || -e "$raw" ]]; then
  print -u2 -- "V4 raw output already exists; preserve it and use a new version"
  exit 98
fi
mkdir -p "$raw"
: > "$study/RUN_V4.jsonl"
: > "$study/METERS_V4.jsonl"
typeset -a arms
arms=(source isogeny-image destination-direct full-image source-bsgs rho-negation rho-negation-aa)
failed=0
for round in 0 1 2; do
  for offset in 0 1 2 3 4 5 6; do
    index=$(( (round + offset) % 7 + 1 ))
    arm=${arms[$index]}
    stem=$raw/round-$round-$arm
    if gtimeout 120 "$meter_binary" "$binary" "$public" "$round" "$arm" \
      "$stem.stdout" "$stem.stderr" "$stem.meter.json" > "$stem.meter.stdout" 2> "$stem.meter.stderr"; then
      run_exit_code=0
    else
      run_exit_code=$?
      failed=1
    fi
    print -r -- "$run_exit_code" > "$stem.exit"
    cat "$stem.stdout" >> "$study/RUN_V4.jsonl"
    if [[ -e "$stem.meter.json" ]]; then
      cat "$stem.meter.json" >> "$study/METERS_V4.jsonl"
      print -- >> "$study/METERS_V4.jsonl"
    fi
  done
done

if (( failed )); then
  print -u2 -- "one or more batch arm processes failed; raw receipts preserved"
  exit 1
fi

"$binary" batch-replay "$public" \
  "$repo/research/isogeny_torsion_fibers_20261012/different_j_full_image_20261015/PUBLIC_V2.json" \
  "$checker" "$study/RUN_V4.jsonl" > "$study/REPLAY_V4.jsonl" 2> "$study/REPLAY_V4.stderr"
