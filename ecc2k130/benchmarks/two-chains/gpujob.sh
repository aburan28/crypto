#!/bin/bash
# Build, verify and bench job behind benchmarks/two-chains/summary.json
# (TWO-CHAINS.md).  Runs inside nvidia/cuda:13.3.1-devel-ubuntu24.04 with the
# ecc2k130 tree mounted at /work and an output directory at /results (or
# $RESULTS).  Launchers, all running this same script on one RTX PRO 6000:
#   modal run modal_job.py --job benchmarks/two-chains/gpujob.sh --out DIR
#   python3 runpod_job.py --job benchmarks/two-chains/gpujob.sh --out DIR
#   python3 aws/bench_job.py --job benchmarks/two-chains/gpujob.sh --out DIR
#
# Every binary is the 20 B/s build of ONE-BLOCK-GEOMETRY.md plus what its name
# says.  Verification forces the same walk count on every binary (1,540,096
# walks, the automatic count of the reference) so the distinguished-point sets
# must be identical across binaries, not just re-walk correctly.
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
  mv ecc2k130 "ecc2k130-$name"
}
build ref                gpu-rtx-pro6000-20b
build ref-alusqr         gpu-rtx-pro6000-20b PACKED_ALU_SQR=1
build c2-256x32          gpu-rtx-pro6000-chains2
build c2-384x32          gpu-rtx-pro6000-chains2 CHAINS2_THREADS=384
build c2-512x16          gpu-rtx-pro6000-chains2 CHAINS2_THREADS=512 CHAINS2_BATCH=16
build c2-256x32-alusqr   gpu-rtx-pro6000-chains2 PACKED_ALU_SQR=1
build ref-prof           gpu-rtx-pro6000-20b PHASE_PROFILE=1
build c2-256x32-prof     gpu-rtx-pro6000-chains2 PHASE_PROFILE=1

# Verification: 300 device reports re-walked by the host reference, and the
# whole distinguished-point set written out for the cross-binary identity check.
verify() {
  local b=$1 threads=$2
  [ -x "ecc2k130-$b" ] || return
  echo "=== verify $b (threads $threads)"
  ./ecc2k130-$b --curve 131 --packed --threads "$threads" --dp-weight 48 --dp-cap 262144 \
      --steps 96 --launches 6 --verify 300 --run-id 7 --dp-file "$R/dp-$b.bin" \
      > "$R/verify-$b.log" 2>&1
  grep -E "MISMATCH|finished|packed chains|resident|registers|persist|OVERFLOW" "$R/verify-$b.log" | tee "$R/verify-$b.txt"
}
verify ref 96256
verify ref-alusqr 96256
verify c2-256x32 48128
verify c2-384x32 48128
verify c2-512x16 96256
verify c2-256x32-alusqr 48128
python3 - "$R" <<'PY' | tee "$R/dp-identity.txt"
import hashlib, os, sys
r = sys.argv[1]
ref = None
for name in ("ref", "ref-alusqr", "c2-256x32", "c2-384x32", "c2-512x16", "c2-256x32-alusqr"):
    path = os.path.join(r, "dp-%s.bin" % name)
    if not os.path.exists(path):
        print("%-20s missing" % name); continue
    data = open(path, "rb").read()
    recs = sorted(data[i:i+32] for i in range(0, len(data) - len(data) % 32, 32))
    digest = hashlib.sha256(b"".join(recs)).hexdigest()
    if ref is None: ref = digest
    print("%-20s %8d records  sha256 %s  %s" % (name, len(recs), digest[:16], "IDENTICAL to ref" if digest == ref else "DIFFERS from ref"))
PY

# Rate: alternating binaries, REPS rounds, automatic worker count.
sample() {
  local b=$1
  [ -x "ecc2k130-$b" ] || return
  ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 64 --verify 0 2>&1 \
      | grep -E "finished|resident" | tr '\n' ' '
  nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader
}
for rep in $(seq 1 "$REPS"); do
  for b in ref c2-256x32 c2-384x32 c2-512x16 ref-alusqr c2-256x32-alusqr; do
    echo "=== bench $b rep $rep"
    sample "$b" | tee -a "$R/bench.txt"
  done
done

# Phase profile of the reference and the two-chain kernel (cycles per warp-step).
for b in ref-prof c2-256x32-prof; do
  [ -x "ecc2k130-$b" ] || continue
  echo "=== phase profile $b"
  ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 16 --verify 0 2>&1 \
      | grep -E "finished|phase profile" | tee "$R/profile-$b.txt"
done
echo "=== done"
