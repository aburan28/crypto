#!/usr/bin/env bash
# Thin orchestration only. Specs contain the frozen scientific choices.
set -u
study=research/conductor_gap_cm_orbits_20261015
mkdir -p "$study/results/replication-claims" "$study/campaign-logs"
for spec in "$study"/prepared/replication/*.json; do
  task=$(basename "$spec" .json)
  case "$task" in control-*) continue;; esac
  if ! mkdir "$study/results/replication-claims/$task" 2>/dev/null; then continue; fi
  target/release/examples/conductor_gap_ecbench run --spec "$spec" \
    --out "$study/sessions/rep-$task" --cpus none \
    --lock "/tmp/cm-rep-$task.lock" --allow-busy --settle 0.1 --quiet \
    > "$study/campaign-logs/rep-$task.log" 2>&1
  status=$?
  printf '%s\n' "$status" > "$study/results/replication-claims/$task/exit"
done
