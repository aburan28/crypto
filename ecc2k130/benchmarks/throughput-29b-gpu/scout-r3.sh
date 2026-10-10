#!/usr/bin/env bash
# Round-3 scouts on this RTX PRO 6000. Control is ALU_SQUARE (17.414).
# Not the frozen recipe; knobs stay off in Make defaults.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
OUT_DIR="benchmarks/throughput-29b-gpu"
mkdir -p "$OUT_DIR"

if [ -x /usr/local/cuda-13.3/bin/nvcc ]; then
  export PATH="/usr/local/cuda-13.3/bin:$PATH"
  export NVCC="/usr/local/cuda-13.3/bin/nvcc"
  export CUDA_HOME="/usr/local/cuda-13.3"
fi

export ARCH="-gencode arch=compute_120,code=sm_120"
preset() {
  make \
    "ARCH=${ARCH}" BATCH=16 THREADS=256 MINBLOCKS=2 \
    PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 \
    PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 \
    PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1 \
    PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256 \
    PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 \
    PACKED_TOP_CLMAD=0 PACKED_PAIR_ILP=1 PACKED_L2_PERSIST=1 UNROLL_SLOTS=2 \
    PACKED_SLOT_PREFETCH=0 TABLE_GLOBAL=0 TABLE_ADDEND_GLOBAL=0 \
    WALK_TABLE=1 TABLE_PIVOT_BYTES=1 PACKED_ALU_SQUARE=1 \
    "$@"
}

sudo -n nvidia-smi -lgc 2430,2430 >/tmp/ecc2k130-lgc.txt 2>&1 || true
sudo -n nvidia-smi -pl 600 >/tmp/ecc2k130-pl.txt 2>&1 || true

echo "=== device table-walk probe (union of new knobs) ==="
preset test-table-walk-cuda PACKED_ALU_SQR=1 PACKED_PAIR_CLMUL=1 PACKED_CLMUL_FLAT=1 \
  >"$OUT_DIR/probe-r3.log" 2>&1
probe_rc=$?
tail -20 "$OUT_DIR/probe-r3.log"
if [ "$probe_rc" -ne 0 ]; then
  echo "device probe failed" >&2
  exit 1
fi

build_one() {
  local name="$1"
  shift
  echo "=== build $name ==="
  preset -B ecc2k130 "$@" >"$OUT_DIR/build-${name}.log" 2>&1
  local rc=$?
  grep -E 'used [0-9]+ registers|stack frame|bytes spill' "$OUT_DIR/build-${name}.log" | tail -20
  if [ "$rc" -ne 0 ]; then
    echo "build $name failed" >&2
    tail -40 "$OUT_DIR/build-${name}.log" >&2
    return 1
  fi
  mv -f ecc2k130 "ecc2k130-${name}"
}

build_one alu_square_r3 || exit 1
build_one alu_sqr PACKED_ALU_SQR=1 || exit 1
build_one pair_clmul PACKED_PAIR_CLMUL=1 || exit 1
build_one clmul_flat PACKED_CLMUL_FLAT=1 || exit 1
build_one unroll4 UNROLL_SLOTS=4 || exit 1
build_one combo PACKED_ALU_SQR=1 PACKED_PAIR_CLMUL=1 PACKED_CLMUL_FLAT=1 || exit 1

python3 - "$OUT_DIR" <<'PY'
import os, subprocess, sys
sys.path.insert(0, ".")
from codegen.benchreport import benchResult, reportsVerified

out_dir = sys.argv[1]
names = ["alu_square_r3", "alu_sqr", "pair_clmul", "clmul_flat", "unroll4", "combo"]
for name in names:
    binp = f"./ecc2k130-{name}"
    cmd = f"{binp} --curve 131 --packed --dp-weight 50 --steps 16 --launches 6 --dp-cap 262144 --verify 300"
    print(f"=== verify {name} ===", flush=True)
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    sys.stdout.write(p.stdout[-1500:])
    ok = reportsVerified(p.returncode, p.stdout, required=300)
    print(f"verify {name}: {'OK' if ok else 'FAIL'}", flush=True)
    if not ok:
        sys.stderr.write(p.stderr[-2000:])
        sys.exit(1)

print("=== one-sample benches (1024 x 32) ===", flush=True)
rows = []
for name in names:
    binp = f"./ecc2k130-{name}"
    cmd = f"{binp} --curve 131 --packed --bench --steps 1024 --launches 32 --verify 0"
    print(f"--- bench {name} ---", flush=True)
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    smi = subprocess.run(
        "nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader",
        shell=True, capture_output=True, text=True)
    r = benchResult(cmd, p.returncode, p.stdout)
    bps = None if not r.get("valid") else r["rate"] / 1000.0
    ident = [ln for ln in p.stdout.splitlines() if ln.startswith("packed ")]
    print("\n".join(ident[-12:]))
    print(f"smi: {smi.stdout.strip()}")
    print(f"{name}: {bps} B/s valid={r.get('valid')}", flush=True)
    rows.append((name, bps, r.get("valid"), ident))

print("\nSCOUT TABLE")
for name, bps, valid, _ in rows:
    print(f"  {name:16s} {bps}  valid={valid}")
PY
