#!/bin/bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
sha256sum -c FROZEN.sha256 --quiet
if [[ -e STOP ]]; then echo 'STOP exists; remove only when intentionally resuming'; exit 1; fi
mkdir -p logs
rho_family() {
  "$root/worker" "$root" rho 0 100 > "$root/logs/rho_coverage_0.log" 2>&1 &
  first=$!
  "$root/worker" "$root" rho 1 100 > "$root/logs/rho_coverage_1.log" 2>&1 &
  second=$!
  wait "$first"; wait "$second"
  "$root/worker" "$root" rho 0 2 > "$root/logs/rho_full_0.log" 2>&1 &
  first=$!
  "$root/worker" "$root" rho 1 2 > "$root/logs/rho_full_1.log" 2>&1 &
  second=$!
  wait "$first"; wait "$second"
}
bsgs_family() {
  "$root/worker" "$root" bsgs 0 20 > "$root/logs/bsgs_coverage_0.log" 2>&1 &
  first=$!
  "$root/worker" "$root" bsgs 10 20 > "$root/logs/bsgs_coverage_1.log" 2>&1 &
  second=$!
  wait "$first"; wait "$second"
  "$root/worker" "$root" bsgs 0 2 > "$root/logs/bsgs_full_0.log" 2>&1 &
  first=$!
  "$root/worker" "$root" bsgs 1 2 > "$root/logs/bsgs_full_1.log" 2>&1 &
  second=$!
  wait "$first"; wait "$second"
}
rho_family &
rho_pid=$!
bsgs_family &
bsgs_pid=$!
printf '{"rho_family_pid":%s,"bsgs_family_pid":%s,"expected_native_workers":4}\n' "$rho_pid" "$bsgs_pid" > "$root/LAUNCHED.json"
wait "$rho_pid"; wait "$bsgs_pid"
echo 'Rho and BSGS workers finished; inspect every status before declaring campaign completion.'
