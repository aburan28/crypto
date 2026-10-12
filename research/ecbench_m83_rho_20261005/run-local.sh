#!/usr/bin/env bash
# PROTOCOL.md amendment 1: part B on the development host.  Eight one-target
# m = 83 sessions side by side, kept awake by caffeinate; each writes its own
# sealed session.  Thin orchestration only.
set -euo pipefail
cd "$(dirname "$0")/../.."
R=research/ecbench_m83_rho_20261005
mkdir -p "$R/host"
{
  echo "commit $(git rev-parse HEAD)"
  rustc -Vv
  uname -a
  sysctl -n machdep.cpu.brand_string hw.perflevel0.physicalcpu hw.perflevel1.physicalcpu hw.memsize
  shasum -a 256 target/release/ecbench
  shasum -a 256 "$R"/specs/*.json
} > "$R/host/PROVENANCE.txt"
./target/release/ecbench host > "$R/host/HOST.json"
pids=()
for i in 1 2 3 4 5 6 7 8; do
  ./target/release/ecbench run --spec "$R/specs/target-$i.json" \
      --out "$R/sessions/m83-target-$i" --cpus none \
      --lock "/tmp/ecbench-m83-$i.lock" > "$R/host/run-target-$i.log" 2>&1 &
  pids+=($!)
done
echo "started ${#pids[@]} sessions: ${pids[*]}"
wait
date -u > "$R/host/DONE"
echo "all sessions ended"
