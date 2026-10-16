#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
study=research/conductor_gap_followups_20261016
mkdir -p /workspace/native-deps
RUSTUP_HOME=/workspace/.rustup CARGO_HOME=/workspace/.cargo /workspace/.cargo/bin/cargo +1.98 build --release --example conductor_gap_followups
cp target/release/examples/conductor_gap_followups /workspace/native-deps/cm-followups-v1
for degree in 5 11; do
  for kind in events noevents; do
    c++ -O3 -std=c++17 -pthread -DCM_Q="$degree" -I research/conductor_gap_cm_orbits_20261015/native "$study/native/folded_$kind.cpp" -o "/workspace/native-deps/cm-folded-$kind-q$degree"
  done
done
