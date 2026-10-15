#!/usr/bin/env bash
set -euo pipefail
study=research/conductor_gap_cm_orbits_20261015
build_dir=/workspace/native-deps/cm-cold-orbits
mkdir -p "$build_dir" "$study/results/cold-passes"
sed '/^int main(/,$d' "$study/native/folded_recurrence.cpp" > "$build_dir/cm-orbit-kernel.hpp"
g++ -std=c++20 -O3 -pthread -I"$build_dir" -I"$study/native" "$study/native/orbit_process.cpp" -o "$build_dir/orbit-process"
sha256sum "$build_dir/cm-orbit-kernel.hpp" "$build_dir/orbit-process" \
  "$study/native/orbit_process.cpp" "$study/native/folded_recurrence.cpp" \
  "$study/native/run_cold_orbits.sh" "$study/prepared/recurrence-orbits.csv" \
  "$study/prepared/cm-params.txt" > "$study/results/cold-pre-run.sha256"
for pass in $(seq 0 19); do
  /workspace/native-deps/root/usr/bin/time -v -o "$study/results/cold-passes/$pass-time.txt" \
    "$build_dir/orbit-process" "$study/prepared/recurrence-orbits.csv" \
    "$study/prepared/cm-params.txt" "$pass" /tmp \
    > "$study/results/cold-passes/$pass.jsonl"
done
cat "$study"/results/cold-passes/*.jsonl > "$study/results/cold-orbits.jsonl"
