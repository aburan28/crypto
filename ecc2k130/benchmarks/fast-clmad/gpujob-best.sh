#!/bin/bash
# The best B200 row: TOP_CLMAD + the CLMAD squaring + the ONB inversion, against
# the reference and the TOP_CLMAD + CLMAD-squaring row (TWO-CHAINS.md section 6.4).
#   modal run --detach modal_job.py --job benchmarks/fast-clmad/gpujob-best.sh --out DIR --gpu B200
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
fail=0
build() {
  local name=$1; shift
  echo "=== build $name: $*"
  make -s "$@" "$ARCHFLAG" > "$R/build-$name.log" 2>&1
  grep -A2 "eccPacked131.*walk" "$R/build-$name.log" | grep -E "registers|spill" | head -2 | tee "$R/build-$name.txt"
  if [ ! -x ecc2k130 ]; then echo "BUILD FAILED $name" | tee -a "$R/failures.txt"; tail -20 "$R/build-$name.log"; return 1; fi
  mv ecc2k130 "ecc2k130-$name"
}
VARIANTS="ref topclmad-clsq topclmad-clsq-onbinv"
build ref            gpu-rtx-pro6000-20b || fail=1


build topclmad-clsq  gpu-rtx-pro6000-20b PACKED_TOP_CLMAD=1 PACKED_ALU_SQUARE=0 || fail=1
build topclmad-clsq-onbinv gpu-rtx-pro6000-20b PACKED_TOP_CLMAD=1 PACKED_ALU_SQUARE=0 PACKED_ONB_INV=1 || fail=1

SMS=$(./ecc2k130-ref --curve 131 --packed --bench --steps 1 --launches 1 --verify 0 2>&1 | grep -oE '[0-9]+ SMs' | grep -oE '[0-9]+' | head -1)
T512=$((SMS * 512))
echo "sms: $SMS, verify threads $T512" | tee -a "$R/host.txt"
verify() {
  local b=$1 threads=$2 want=$3 status=0
  [ -x "ecc2k130-$b" ] || return 1
  echo "=== verify $b (threads $threads, expect alu square $want)"
  ./ecc2k130-$b --curve 131 --packed --threads "$threads" --dp-weight 48 --dp-cap 262144 \
      --steps 96 --launches 6 --verify 300 --run-id 7 --dp-file "$R/dp-$b.bin" \
      > "$R/verify-$b.log" 2>&1 || status=$?
  grep -E "MISMATCH|finished|packed alu square|packed top clmad|packed onb inv|resident|registers|OVERFLOW" "$R/verify-$b.log" | tee "$R/verify-$b.txt" || true
  if ! grep -q "packed alu square: $want" "$R/verify-$b.log"; then
    echo "KNOB MISMATCH: $b did not build with PACKED_ALU_SQUARE=$want" | tee -a "$R/failures.txt"; status=1
  fi
  return "$status"
}
verify ref $T512 1 || fail=1


verify topclmad-clsq $T512 0 || fail=1
verify topclmad-clsq-onbinv $T512 0 || fail=1
python3 - "$R" $VARIANTS <<'PY' | tee "$R/dp-identity.txt" || fail=1
import hashlib, os, sys
r = sys.argv[1]
ref = None
bad = 0
for name in sys.argv[2:]:
    path = os.path.join(r, "dp-%s.bin" % name)
    if not os.path.exists(path):
        print("%-20s missing" % name)
        bad = 1
        continue
    data = open(path, "rb").read()
    recs = sorted(data[i:i+32] for i in range(0, len(data) - len(data) % 32, 32))
    digest = hashlib.sha256(b"".join(recs)).hexdigest()
    if ref is None: ref = digest
    if digest != ref: bad = 1
    print("%-20s %8d records  sha256 %s  %s" % (name, len(recs), digest[:16], "IDENTICAL to ref" if digest == ref else "DIFFERS from ref"))
sys.exit(bad)
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
echo "=== done"
exit "$fail"
