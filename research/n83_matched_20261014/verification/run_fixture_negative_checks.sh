#!/usr/bin/env bash
set -eu

# Native programs make every arithmetic/replay decision. This shell script only
# creates altered file copies and retains the program exit statuses and logs.
repo_dir=$(pwd)
checker=/workspace/n83_matched/verification/check_fixtures
generator=/workspace/n83_matched/verification/generate_targets
fixture_dir="$repo_dir/research/n83_matched_20261014/verification/fixtures"
archive="$repo_dir/experiments/koblitz-single-target-n83-20261006"
prior="$repo_dir/research/n83_evaluation_20261012/optimization/verification/fixtures"
speed="$repo_dir/research/n83_speed_20261012/verification/fixtures"
scratch=/workspace/n83_matched/verification/negative_checks
retained="$repo_dir/research/n83_matched_20261014/verification/negative_checks"
test ! -e "$scratch"
test ! -e "$retained"
mkdir -p "$scratch" "$retained"

set +e
"$generator" "$archive" "$prior" "$speed" "$fixture_dir" > "$retained/existing_directory.stdout" 2> "$retained/existing_directory.stderr"
existing_dir_status=$?
"$checker" "$archive" "$prior" "$speed" "$fixture_dir" "$repo_dir/research/n83_matched_20261014/verification/fixture_replay_receipt.json" > "$retained/existing_receipt.stdout" 2> "$retained/existing_receipt.stderr"
existing_receipt_status=$?
set -e

cp -R "$fixture_dir" "$scratch/wrong_point_fixtures"
cp "$scratch/wrong_point_fixtures/target_02.jsonl" "$scratch/wrong_point_fixtures/target_01.jsonl"
set +e
"$checker" "$archive" "$prior" "$speed" "$scratch/wrong_point_fixtures" "$scratch/wrong_point_receipt.json" > "$retained/wrong_point.stdout" 2> "$retained/wrong_point.stderr"
wrong_point_status=$?
set -e

cp -R "$fixture_dir" "$scratch/wrong_scalar_fixtures"
jq '.entries[0].scalar_validation_only="1"' "$scratch/wrong_scalar_fixtures/validation_sidecar.json" > "$scratch/altered_sidecar.json"
mv "$scratch/altered_sidecar.json" "$scratch/wrong_scalar_fixtures/validation_sidecar.json"
read -r altered_sidecar_hash ignored_path < <(sha256sum "$scratch/wrong_scalar_fixtures/validation_sidecar.json")
jq --arg hash "$altered_sidecar_hash" '.validation_sidecar_sha256=$hash' "$scratch/wrong_scalar_fixtures/generator_receipt.json" > "$scratch/altered_generator.json"
mv "$scratch/altered_generator.json" "$scratch/wrong_scalar_fixtures/generator_receipt.json"
set +e
"$checker" "$archive" "$prior" "$speed" "$scratch/wrong_scalar_fixtures" "$scratch/wrong_scalar_receipt.json" > "$retained/wrong_scalar.stdout" 2> "$retained/wrong_scalar.stderr"
wrong_scalar_status=$?
set -e

test "$existing_dir_status" -ne 0
test "$existing_receipt_status" -ne 0
test "$wrong_point_status" -ne 0
test "$wrong_scalar_status" -ne 0
rg -q 'Refusing existing output directory' "$retained/existing_directory.stderr"
rg -q 'Refusing existing verification receipt' "$retained/existing_receipt.stderr"
rg -q 'Matched point-only file bytes/schema mismatch' "$retained/wrong_point.stderr"
rg -q 'Matched validation-only scalar replay failure' "$retained/wrong_scalar.stderr"
jq -n --argjson existing_dir "$existing_dir_status" --argjson existing_receipt "$existing_receipt_status" --argjson wrong_point "$wrong_point_status" --argjson wrong_scalar "$wrong_scalar_status" \
  '{kind:"native_fixture_negative_checks",status:"PASS",existing_output_directory_rejected:$existing_dir,existing_verification_receipt_rejected:$existing_receipt,wrong_point_bytes_rejected:$wrong_point,wrong_scalar_rejected_after_updating_sidecar_hash:$wrong_scalar,arithmetic_and_replay_decisions:"native C++ checker; shell/jq prepare copies only"}' > "$retained/receipt.json"
