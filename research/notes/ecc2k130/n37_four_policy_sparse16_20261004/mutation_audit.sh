#!/usr/bin/env bash
# Run from the repository root after building both release replay examples.
set -euo pipefail

note=research/notes/ecc2k130/n37_four_policy_sparse16_20261004
raw=$note/RESULT.json.gz
replay=target/release/examples/n37_four_policy_sparse16_pdp_replay
input_replay=target/release/examples/n37_four_policy_sparse16_inputs_replay
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
receipts=${1:-$note/mutations}
if test -e "$receipts"; then
  echo "refuse to overwrite mutation receipts: $receipts" >&2
  exit 1
fi
mkdir -p "$receipts"
gzip -dc "$raw" > "$scratch/result.json"

run_mutation() {
  local name=$1 filter=$2
  jq "$filter" "$scratch/result.json" > "$scratch/$name.json"
  if "$replay" "$scratch/$name.json" "$receipts/$name.receipt.json" \
      > "$scratch/$name.stdout" 2> "$scratch/$name.stderr"; then
    echo "mutation $name unexpectedly passed" >&2
    exit 1
  fi
  jq -e '.status == "FAIL" and (.error | length > 0)' \
    "$receipts/$name.receipt.json" > /dev/null
}

# Mutate both paired arms to test arithmetic replay beyond pairwise equality.
run_mutation witness '(.policies[0].targets[0].m3.indices[0]) += 1 | (.policies[1].targets[0].m3.indices[0]) += 1'
run_mutation rank_row '(.policies[0].rank_trace[0].row[0]) += 1 | (.policies[1].rank_trace[0].row[0]) += 1'
run_mutation hit_flag '.policies[0].targets[0].m3.status = "proved_miss" | .policies[1].targets[0].m3.status = "proved_miss"'
run_mutation decision '.selection_comparison.decision = "FIXED_BLOCK_M3_SELECTION_LEAD"'

cp -R "$note/inputs" "$scratch/inputs"
sed -n '1p' "$note/inputs/points-b0.jsonl" \
  | jq -c '[.[0] + 1, .[1]]' > "$scratch/first-point.json"
tail -n +2 "$scratch/inputs/points-b0.jsonl" > "$scratch/remaining-points.jsonl"
cat "$scratch/first-point.json" "$scratch/remaining-points.jsonl" \
  > "$scratch/inputs/points-b0.jsonl"
if "$input_replay" "$scratch/inputs" "$scratch/input-replay.json" \
    > "$scratch/target.stdout" 2> "$receipts/target.stderr.txt"; then
  echo "mutated target input unexpectedly passed" >&2
  exit 1
fi
if ! rg -q 'SHA-256 mismatch' "$receipts/target.stderr.txt"; then
  echo "mutated target did not fail the input digest gate" >&2
  exit 1
fi
echo "PASS witness rank_row hit_flag decision target"
