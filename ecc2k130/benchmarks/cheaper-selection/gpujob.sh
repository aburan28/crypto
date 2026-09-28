#!/bin/bash
# Build, verify and bench TABLE_PHASE_POPC against the 20 B/s reference
# (CHEAPER-SELECTION.md).  Runs inside nvidia/cuda:13.3.1-devel-ubuntu24.04
# with the ecc2k130 tree at /work and results at /results (or $RESULTS).
#
#   modal run --detach modal_job.py --job benchmarks/cheaper-selection/gpujob.sh \
#       --out DIR --gpu "RTX PRO 6000"
#   python3 runpod_job.py --job benchmarks/cheaper-selection/gpujob.sh --out DIR
#   python3 aws/bench_job.py --job benchmarks/cheaper-selection/gpujob.sh --out DIR
#
# Falsification (CHEAPER-SELECTION.md §2): paired median ≥ 20.50 B/s vs the
# same-session reference, every pair favouring the candidate, 300/300 re-walks,
# byte-identical DP set on a forced common walk count. Abandon below 20.20 B/s.
set -uo pipefail
cd /work || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq make g++ python3 >/dev/null 2>&1
R=${RESULTS:-/results}; mkdir -p "$R"
REPS=${REPS:-5}
{
  nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,power.limit,memory.total --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
} | tee "$R/host.txt"

build() {
  local name=$1; shift
  echo "=== build $name: $*"
  make -s "$@" > "$R/build-$name.log" 2>&1
  grep -A2 "eccPacked131.*walk" "$R/build-$name.log" | grep -E "registers|spill" | head -2 | tee "$R/build-$name.txt"
  if [ ! -x ecc2k130 ]; then echo "BUILD FAILED $name" | tee -a "$R/failures.txt"; tail -20 "$R/build-$name.log"; return 1; fi
  # Identity line must carry the phase knob so a silent default cannot fake a pair.
  ./ecc2k130 --curve 131 --packed --threads 256 --steps 1 --launches 1 --verify 0 --bench \
      > "$R/identity-$name.log" 2>&1 || true
  mv ecc2k130 "ecc2k130-$name"
}
fail=0
build ref     gpu-rtx-pro6000-20b TABLE_PHASE_POPC=0 || fail=1
build popc    gpu-rtx-pro6000-20b TABLE_PHASE_POPC=1 || fail=1
build ref-prof  gpu-rtx-pro6000-20b TABLE_PHASE_POPC=0 PHASE_PROFILE=1 || fail=1
build popc-prof gpu-rtx-pro6000-20b TABLE_PHASE_POPC=1 PHASE_PROFILE=1 || fail=1

# Assert the knob from each binary's own identity line before timing it.
for pair in "ref:0" "popc:1" "ref-prof:0" "popc-prof:1"; do
  b=${pair%%:*}; want=${pair##*:}
  [ -x "ecc2k130-$b" ] || continue
  # Re-run a one-shot report to print identity (bench above may have been empty).
  ./ecc2k130-$b --curve 131 --packed --threads 256 --steps 1 --launches 1 --verify 0 --bench \
      > "$R/identity-$b.log" 2>&1 || true
  if ! grep -q "packed table phase popc: $want" "$R/identity-$b.log"; then
    echo "KNOB MISMATCH: $b did not report phase popc=$want" | tee -a "$R/failures.txt"
    grep -E "phase popc|table pivot|table walk" "$R/identity-$b.log" | tee -a "$R/failures.txt" || true
    fail=1
  fi
done

verify() {
  local b=$1 status=0
  [ -x "ecc2k130-$b" ] || return 1
  echo "=== verify $b"
  ./ecc2k130-$b --curve 131 --packed --threads 96256 --dp-weight 48 --dp-cap 262144 \
      --steps 96 --launches 6 --verify 300 --run-id 7 --dp-file "$R/dp-$b.bin" \
      > "$R/verify-$b.log" 2>&1 || status=$?
  grep -E "MISMATCH|finished|registers|persist|OVERFLOW|dropped" "$R/verify-$b.log" | tee "$R/verify-$b.txt" || true
  return "$status"
}
verify ref || fail=1
verify popc || fail=1

python3 - "$R" <<'PY' | tee "$R/dp-identity.txt" || fail=1
import hashlib, os, sys
r = sys.argv[1]
ref = None
bad = 0
for name in ("ref", "popc"):
    path = os.path.join(r, "dp-%s.bin" % name)
    if not os.path.exists(path):
        print("%-8s missing" % name); bad = 1; continue
    data = open(path, "rb").read()
    recs = sorted(data[i:i+32] for i in range(0, len(data) - len(data) % 32, 32))
    digest = hashlib.sha256(b"".join(recs)).hexdigest()
    if ref is None: ref = digest
    if digest != ref: bad = 1
    print("%-8s %8d records  sha256 %s  %s" % (
        name, len(recs), digest[:16], "IDENTICAL to ref" if digest == ref else "DIFFERS from ref"))
sys.exit(bad)
PY

bench() {
  local b=$1 i=$2
  [ -x "ecc2k130-$b" ] || return 1
  ./ecc2k130-$b --curve 131 --packed --threads 96256 --bench --steps 1024 --launches 32 --verify 0 \
      > "$R/bench-$b-$i.log" 2>&1
  grep -E "finished:|registers|spill" "$R/bench-$b-$i.log" | tee -a "$R/bench-$b.txt" || true
}
echo "=== alternating benches (REPS=$REPS)" | tee "$R/bench-ref.txt" > "$R/bench-popc.txt"
for i in $(seq 1 "$REPS"); do
  bench ref "$i" || fail=1
  bench popc "$i" || fail=1
done
# Profile once each (not timed for the verdict).
for b in ref-prof popc-prof; do
  [ -x "ecc2k130-$b" ] || continue
  echo "=== profile $b"
  ./ecc2k130-$b --curve 131 --packed --threads 96256 --bench --steps 1024 --launches 8 --verify 0 \
      > "$R/profile-$b.log" 2>&1 || fail=1
  grep -E "phase|forward|invers|reverse|finished" "$R/profile-$b.log" | tee "$R/profile-$b.txt" || true
done

python3 - "$R" "$REPS" <<'PY' | tee "$R/summary.txt"
import re, statistics, sys
r, reps = sys.argv[1], int(sys.argv[2])
def rates(name):
    out = []
    for i in range(1, reps + 1):
        text = open("%s/bench-%s-%d.log" % (r, name, i)).read()
        # Client prints millions; ONE-BLOCK quotes billions (M/1000).
        m = re.search(r"finished:\s*([0-9.]+)\s*M it/s", text)
        if m: out.append(float(m.group(1)) / 1000.0)
    return out
ref, popc = rates("ref"), rates("popc")
print("ref samples B/s:", ["%.3f" % x for x in ref])
print("popc samples B/s:", ["%.3f" % x for x in popc])
if not ref or not popc or len(ref) != len(popc):
    print("INCOMPLETE"); sys.exit(1)
pairs = [p / a for a, p in zip(ref, popc)]
med_ref, med_popc = statistics.median(ref), statistics.median(popc)
print("median ref %.3f  popc %.3f  ratio %.4f" % (med_ref, med_popc, med_popc / med_ref))
print("pair ratios:", ["%.4f" % x for x in pairs], "all>1", all(x > 1 for x in pairs))
target, abandon = 20.50, 20.20
if med_popc >= target and all(x > 1 for x in pairs):
    print("VERDICT: SUCCESS (clears %.2f B/s and every pair)" % target)
elif med_popc < abandon:
    print("VERDICT: ABANDON (below %.2f B/s)" % abandon)
else:
    print("VERDICT: INCONCLUSIVE (between abandon and success bars)")
PY

exit "$fail"
