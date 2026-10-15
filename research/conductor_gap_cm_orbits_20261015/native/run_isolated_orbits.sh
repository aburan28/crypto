#!/usr/bin/env bash
set -euo pipefail
study=research/conductor_gap_cm_orbits_20261015
exec /workspace/native-deps/cm-folded-recurrence \
  "$study/prepared/recurrence-orbits.csv" 8192 1 -1 1 /tmp \
  "$study/prepared/cm-params.txt" > "$study/results/isolated-orbits.jsonl"
