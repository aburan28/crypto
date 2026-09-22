#!/bin/bash
# The table-walk kernel on a part whose carry-less unit is not 1/64-rate
# (TWO-CHAINS.md section 6; receipts in benchmarks/fast-clmad/).  Builds the
# RTX PRO 6000 20 B/s knob set for whatever GPU the container has, plus the
# ALU-to-CLMAD trades that lost on the 6000 because its unit was full, verifies
# every binary (300 host re-walks, distinguished-point set identity across
# binaries at a forced common walk count) and benches them alternating.
#   modal run modal_job.py --job benchmarks/fast-clmad/gpujob.sh --out DIR --gpu B200
set -uo pipefail
cd /work || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq make g++ python3 >/dev/null 2>&1
R=${RESULTS:-/results}; mkdir -p "$R"
REPS=${REPS:-5}
{
  nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-unknown}"
} | tee "$R/host.txt"
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCHFLAG="PRO6000_ARCH=-gencode arch=compute_${CAP},code=sm_${CAP}"
echo "arch: $ARCHFLAG" | tee -a "$R/host.txt"

build() {
  local name=$1; shift
  echo "=== build $name: $*"
  make -s "$@" "$ARCHFLAG" > "$R/build-$name.log" 2>&1
  grep -A2 "eccPacked131.*walk" "$R/build-$name.log" | grep -E "registers|spill" | head -2 | tee "$R/build-$name.txt"
  if [ ! -x ecc2k130 ]; then echo "BUILD FAILED $name" | tee -a "$R/failures.txt"; tail -20 "$R/build-$name.log"; return 1; fi
  mv ecc2k130 "ecc2k130-$name"
}
VARIANTS="ref clsq topclmad topclmad-clsq onbinv c2-256x32"
build ref            gpu-rtx-pro6000-20b
build clsq           gpu-rtx-pro6000-20b PACKED_ALU_SQUARE=0
build topclmad       gpu-rtx-pro6000-20b PACKED_TOP_CLMAD=1
build topclmad-clsq  gpu-rtx-pro6000-20b PACKED_TOP_CLMAD=1 PACKED_ALU_SQUARE=0
build onbinv         gpu-rtx-pro6000-20b PACKED_ONB_INV=1
build c2-256x32      gpu-rtx-pro6000-chains2
build ref-prof       gpu-rtx-pro6000-20b PHASE_PROFILE=1
build topclmad-clsq-prof gpu-rtx-pro6000-20b PACKED_TOP_CLMAD=1 PACKED_ALU_SQUARE=0 PHASE_PROFILE=1

# The automatic worker count of the 512-thread builds on this part, so every
# binary is verified on the same walks.
SMS=$(./ecc2k130-ref --curve 131 --packed --bench --steps 1 --launches 1 --verify 0 2>&1 | grep -oE '[0-9]+ SMs' | grep -oE '[0-9]+' | head -1)
T512=$((SMS * 512)); T256=$((SMS * 256))
echo "sms: $SMS, verify threads $T512 (512-thread builds) / $T256 (256-thread builds)" | tee -a "$R/host.txt"
verify() {
  local b=$1 threads=$2
  [ -x "ecc2k130-$b" ] || return
  echo "=== verify $b (threads $threads)"
  ./ecc2k130-$b --curve 131 --packed --threads "$threads" --dp-weight 48 --dp-cap 262144 \
      --steps 96 --launches 6 --verify 300 --run-id 7 --dp-file "$R/dp-$b.bin" \
      > "$R/verify-$b.log" 2>&1
  grep -E "MISMATCH|finished|packed chains|resident|registers|persist|OVERFLOW" "$R/verify-$b.log" | tee "$R/verify-$b.txt"
}
for b in ref clsq topclmad topclmad-clsq onbinv; do verify $b $T512; done
verify c2-256x32 $T256
python3 - "$R" $VARIANTS <<'PY' | tee "$R/dp-identity.txt"
import hashlib, os, sys
r = sys.argv[1]
ref = None
for name in sys.argv[2:]:
    path = os.path.join(r, "dp-%s.bin" % name)
    if not os.path.exists(path):
        print("%-20s missing" % name); continue
    data = open(path, "rb").read()
    recs = sorted(data[i:i+32] for i in range(0, len(data) - len(data) % 32, 32))
    digest = hashlib.sha256(b"".join(recs)).hexdigest()
    if ref is None: ref = digest
    print("%-20s %8d records  sha256 %s  %s" % (name, len(recs), digest[:16], "IDENTICAL to ref" if digest == ref else "DIFFERS from ref"))
PY

sample() {
  local b=$1
  [ -x "ecc2k130-$b" ] || return
  ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 64 --verify 0 2>&1 \
      | grep -E "finished|resident" | tr '\n' ' '
  nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader
}
for rep in $(seq 1 "$REPS"); do
  for b in $VARIANTS; do
    echo "=== bench $b rep $rep" | tee -a "$R/bench.txt"
    sample "$b" | tee -a "$R/bench.txt"
  done
done
for b in ref-prof topclmad-clsq-prof; do
  [ -x "ecc2k130-$b" ] || continue
  echo "=== phase profile $b"
  ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 16 --verify 0 2>&1 \
      | grep -E "finished|phase profile" | tee "$R/profile-$b.txt"
done
echo "=== done"
