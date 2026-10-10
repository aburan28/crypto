#!/usr/bin/env bash
# Lead-ladder baseline runner (docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md section 2 and 3).
# Operation counts only: this host cannot earn an ICMS isolation level, so wall time in the
# reports is a practicality note, never a result.
#
#   research/ic_leads_20261010/run_baseline.sh [sweep.json] [n ...]
#
# Runs `ic bench` on K_0 (a = 0) at each rung with the sweep's configurations, writing one JSON
# report per rung under runs/<date>/ and the exact command beside it.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
ic="${IC_BIN:-$root/target/release/ic}"
sweep="${1:-$here/sweeps/baseline_l1_l4.json}"
shift || true
rungs=("$@")
[ ${#rungs[@]} -eq 0 ] && rungs=(17 23 31 41)
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
out="$here/runs/$stamp"
mkdir -p "$out"
[ -x "$ic" ] || { echo "no ic binary at $ic (cargo build --release --bin ic)"; exit 2; }
git -C "$root" rev-parse HEAD > "$out/commit.txt"
shasum -a 256 "$ic" > "$out/ic.sha256"
uptime > "$out/host_load_at_start.txt"
for n in "${rungs[@]}"; do
  tmp="$out/sweep_n$n.json"
  python3 - "$sweep" "$n" "$tmp" <<'EOF'
import json, sys
s = json.load(open(sys.argv[1])); s["instance"]["degree"] = int(sys.argv[2])
json.dump(s, open(sys.argv[3], "w"), indent=1)
EOF
  cmd=("$ic" --out "$out/report_n$n.json" bench --sweep "$tmp" --koblitz-a "${KOBLITZ_A:-0}" --rho-runs 16)
  printf '%q ' "${cmd[@]}" > "$out/cmd_n$n.txt"; echo >> "$out/cmd_n$n.txt"
  echo "== n=$n"; "${cmd[@]}" 2>&1 | tee "$out/stdout_n$n.txt" | tail -20
done
uptime > "$out/host_load_at_end.txt"
echo "reports in $out"
