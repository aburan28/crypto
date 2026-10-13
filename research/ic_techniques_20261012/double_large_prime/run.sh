#!/usr/bin/env bash
# Double-large-prime measurement runner.
#
#   research/ic_techniques_20261012/double_large_prime/run.sh [n ...]
#
# Runs every sweeps/dlp_n<n>_c<columns>.json through `ic bench` on K_0,
# writing report_n<n>_c<columns>.json(.gz), the command, stdout, the
# commit and the binary hash under runs/<stamp>/.  Operation counts are
# the metric; this host is shared and wall time is a practicality note.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../../.." && pwd)"
ic="${IC_BIN:-$root/target/release/ic}"
rungs=("$@")
[ ${#rungs[@]} -eq 0 ] && rungs=(23 31 41)
stamp="${RUN_STAMP:-$(date -u +%Y%m%dT%H%M%SZ)}"
out="$here/runs/$stamp"
mkdir -p "$out"
[ -x "$ic" ] || { echo "no ic binary at $ic (cargo build --release --bin ic)"; exit 2; }
git -C "$root" rev-parse HEAD > "$out/commit.txt"
shasum -a 256 "$ic" > "$out/ic.sha256"
uptime > "$out/host_load_at_start.txt"
for n in "${rungs[@]}"; do
  # SWEEP_PATTERN selects the sweep files of a rung (default: the
  # column cells; "control_envelope" runs the equal-envelope control).
  for sweep in "$here"/sweeps/dlp_n"$n"_"${SWEEP_PATTERN:-c[0-9]*}".json; do
    tag="$(basename "$sweep" .json)"
    cmd=("$ic" --out "$out/report_$tag.json" bench --sweep "$sweep" --koblitz-a 0 --rho-runs "${RHO_RUNS:-8}")
    printf '%q ' "${cmd[@]}" > "$out/cmd_$tag.txt"; echo >> "$out/cmd_$tag.txt"
    echo "== $tag"; "${cmd[@]}" 2>&1 | tee "$out/stdout_$tag.txt" | tail -30
    gzip -9 -n -f "$out/report_$tag.json" "$out/stdout_$tag.txt" 2>/dev/null || true
  done
done
uptime > "$out/host_load_at_end.txt"
echo "reports in $out"
