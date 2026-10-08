#!/usr/bin/env bash
# Thin native build/replay orchestration. Historical Python is never executed.
set -euo pipefail
project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
search_binary=${1:-"$project_root/tools/endomorphism-search/target/release/endomorphism-search"}
archive="$project_root/research/endomorphism_search_20261005/historical/python-source-and-results.zip"
replay_dir=$(mktemp -d)
trap 'rm -rf "$replay_dir"' EXIT

expected_archive=89d4f459ef440dd238864d4a39a52c9b59d2c6e086d6f2566a81aad7f82c83d6
actual_archive=$(sha256sum "$archive")
if [[ ${actual_archive%% *} != "$expected_archive" ]]; then
  echo 'historical archive digest mismatch' >&2
  exit 1
fi
unzip -q "$archive" -d "$replay_dir"
historical_dir="$replay_dir/endomorphism-search"
"$search_binary" replay --historical-dir "$historical_dir" --out "$replay_dir/replay.json"
cat "$replay_dir/replay.json"

# Corrupt evidence without updating its manifest: a nonzero exit is required.
printf '\n' >> "$historical_dir/results/demo.json"
if "$search_binary" replay --historical-dir "$historical_dir" > "$replay_dir/rejected.json" 2> "$replay_dir/rejected.stderr"; then
  echo 'altered evidence was accepted' >&2
  exit 1
fi
if [[ -s "$replay_dir/rejected.json" ]]; then
  echo 'failed replay emitted a success result' >&2
  exit 1
fi
if ! grep -q 'historical hash mismatch: results/demo.json' "$replay_dir/rejected.stderr"; then
  echo 'altered evidence failed for an unexpected reason' >&2
  cat "$replay_dir/rejected.stderr" >&2
  exit 1
fi
printf '\n' >> "$historical_dir/results/artifact-manifest.json"
if "$search_binary" replay --historical-dir "$historical_dir" > "$replay_dir/rejected.json" 2> "$replay_dir/rejected.stderr"; then
  echo 'altered manifest was accepted' >&2
  exit 1
fi
if ! grep -q 'historical manifest digest mismatch' "$replay_dir/rejected.stderr"; then
  echo 'altered manifest failed for an unexpected reason' >&2
  exit 1
fi
echo 'PASS: frozen correctness parity and altered-evidence rejection controls'
