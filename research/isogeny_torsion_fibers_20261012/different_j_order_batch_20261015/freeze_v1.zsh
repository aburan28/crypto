#!/bin/zsh
set -u

repo=/Volumes/SSD990-2/.codex/worktrees/d92d/crypto
study=$repo/research/isogeny_torsion_fibers_20261012/different_j_order_batch_20261015
fiber=$repo/research/isogeny_torsion_fibers_20261012/different_j_order_screen_20261015/FIBER_V1.json
binary=$repo/target/order-batch-check/release/examples/different_j_order_batch
spec=$study/FREEZE_SPEC_V1.json

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
verify_hash "$fiber" '.fiber_receipt_sha256'
verify_hash "$study/freeze_v1.zsh" '.script_sha256'

if [[ -e "$study/PUBLIC_V1.json" || -e "$study/CHECKER_V1.json" || -e "$study/FREEZE_V1.stdout" ]]; then
  print -u2 -- "V1 freeze output exists; preserve it and use a new version"
  exit 98
fi

if gtimeout 120 "$binary" freeze "$fiber" "$study/PUBLIC_V1.json" "$study/CHECKER_V1.json" \
  > "$study/FREEZE_V1.stdout" 2> "$study/FREEZE_V1.stderr"; then
  freeze_exit=0
else
  freeze_exit=$?
fi
print -r -- "$freeze_exit" > "$study/FREEZE_V1.exit"
exit "$freeze_exit"
