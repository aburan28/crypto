#!/bin/bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
exec 9>"$root/launcher.lock"
flock -n 9 || { echo "Another resumable launcher holds the lock"; exit 3; }
if pgrep -f "^$root/worker $root " >/dev/null; then echo "Existing campaign workers are active; no duplicate launch"; exit 3; fi
sha256sum -c FROZEN.sha256 --quiet
if [[ -e STOP ]]; then echo 'STOP exists; remove only when intentionally resuming'; exit 1; fi
if [[ "${1:-}" == --check ]]; then echo "Frozen inputs verified; no workers active"; exit 0; fi
logroot="$root/logs/launch_$(date -u +%Y%m%dT%H%M%SZ)_$"
mkdir -p "$logroot"
rho_family() {
  "$root/worker" "$root" rho 0 100 > "$logroot/rho_coverage_0.log" 2>&1 &
  first=$!
  "$root/worker" "$root" rho 1 100 > "$logroot/rho_coverage_1.log" 2>&1 &
  second=$!
  wait "$first"; wait "$second"
  "$root/worker" "$root" rho 0 2 > "$logroot/rho_full_0.log" 2>&1 &
  first=$!
  "$root/worker" "$root" rho 1 2 > "$logroot/rho_full_1.log" 2>&1 &
  second=$!
  wait "$first"; wait "$second"
}
bsgs_family() {
  "$root/worker" "$root" bsgs 0 20 > "$logroot/bsgs_coverage_0.log" 2>&1 &
  first=$!
  "$root/worker" "$root" bsgs 10 20 > "$logroot/bsgs_coverage_1.log" 2>&1 &
  second=$!
  wait "$first"; wait "$second"
  "$root/worker" "$root" bsgs 0 2 > "$logroot/bsgs_full_0.log" 2>&1 &
  first=$!
  "$root/worker" "$root" bsgs 1 2 > "$logroot/bsgs_full_1.log" 2>&1 &
  second=$!
  wait "$first"; wait "$second"
}
rho_family &
rho_pid=$!
bsgs_family &
bsgs_pid=$!
printf '{"rho_family_pid":%s,"bsgs_family_pid":%s,"expected_native_workers":4}\n' "$rho_pid" "$bsgs_pid" > "$logroot/LAUNCHED.json"
wait "$rho_pid"; wait "$bsgs_pid"
echo 'Rho and BSGS workers finished; inspect every status before declaring campaign completion.'
