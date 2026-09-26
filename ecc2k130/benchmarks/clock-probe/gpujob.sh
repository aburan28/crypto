#!/usr/bin/env bash
# Sample SM clocks while the campaign packed bench runs.
# Usage (from ecc2k130/):
#   python3 modal_job.py --job benchmarks/clock-probe/gpujob.sh --out /tmp/clock-probe \
#     --gpu RTX-PRO-6000 --env ECC_PACKED_CLMAD=1 --env ECC_CUDA_VERSION=13.3.1 ...
set -euo pipefail
cd /work
RESULTS="${RESULTS:-/results}"
mkdir -p "$RESULTS"

# Campaign packed knobs (same as run.sh CURVE=131)
export ECC_CUDA_VERSION="${ECC_CUDA_VERSION:-13.3.1}"
export ECC_PACKED_SINGLE_PRODUCT=1 ECC_PACKED_CACHE_DENOM=1 ECC_PACKED_BY_VALUE=1
export ECC_PACKED_PERM_SIGMA=3 ECC_PACKED_POLY_CHAIN=1 ECC_PACKED_UNROLL_INV=1
export ECC_PACKED_PAIR_PRODUCTS=1 ECC_PACKED_POLY_STATE=1 ECC_PACKED_DIRECT_REDUCE=1
export ECC_PACKED_GENERATED_PRODUCT=1 ECC_PACKED_CLMAD=1 ECC_PACKED_STATE_TILE=256
export ECC_PACKED_WEIGHTED_PREFIX=2 ECC_PACKED_COMPACT_STATE=1 ECC_PACKED_SHARED_SIGMA=1
export ECC_PACKED_TOP_CLMAD=0 ECC_PACKED_CLMAD_SQUARE=0 ECC_PACKED_KARAT3=0
export ECC_WALK_TABLE=0

echo "=== pre-bench nvidia-smi ==="
nvidia-smi --query-gpu=name,driver_version,pstate,clocks.current.sm,clocks.max.sm,clocks.current.memory,power.limit,power.draw,temperature.gpu --format=csv | tee "$RESULTS/pre-smi.csv"

# Background sampler (~4 Hz)
(
  echo "timestamp_s,sm_mhz,mem_mhz,power_w,temp_c,pstate"
  while true; do
    nvidia-smi --query-gpu=timestamp,clocks.current.sm,clocks.current.memory,power.draw,temperature.gpu,pstate --format=csv,noheader,nounits \
      | awk -F', ' '{gsub(/ /,"",$0); print}' || true
    sleep 0.25
  done
) > "$RESULTS/smi-sample.csv" &
SAMPLER=$!
trap 'kill $SAMPLER 2>/dev/null || true' EXIT

echo "=== build packed campaign binary ==="
make -B gpu BATCH=16 THREADS=256 LEAF=0 MINBLOCKS=2 PACKED=1 \
  PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 \
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 \
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1 \
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256 \
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 \
  PACKED_TOP_CLMAD=0 PACKED_CLMAD_SQUARE=0 PACKED_KARAT3=0 \
  WALK_TABLE=0 2>&1 | tee "$RESULTS/build.log"

echo "=== bench automatic workers (survey geometry) ==="
./ecc2k130 --curve 131 --bench --steps 1024 --launches 32 --verify 0 --packed --threads 0 \
  2>&1 | tee "$RESULTS/bench-auto.log"

echo "=== bench 385024 workers (campaign geometry) ==="
./ecc2k130 --curve 131 --bench --steps 1024 --launches 16 --verify 0 --packed --threads 385024 \
  2>&1 | tee "$RESULTS/bench-385k.log"

echo "=== post-bench nvidia-smi ==="
nvidia-smi --query-gpu=name,pstate,clocks.current.sm,clocks.max.sm,power.draw,temperature.gpu --format=csv | tee "$RESULTS/post-smi.csv"

kill $SAMPLER 2>/dev/null || true
wait $SAMPLER 2>/dev/null || true

python3 - <<'PY' | tee "$RESULTS/summary.txt"
import csv, re, pathlib, statistics
root = pathlib.Path("/results")
samples = []
with open(root/"smi-sample.csv") as f:
    r = csv.DictReader(f)
    for row in r:
        try:
            samples.append({
                "sm": float(row["sm_mhz"]),
                "power": float(str(row["power_w"]).replace("W","") or 0),
                "temp": float(row["temp_c"]),
            })
        except Exception:
            pass
sms = [s["sm"] for s in samples if s["sm"] > 0]
powers = [s["power"] for s in samples if s["power"] > 0]
print(f"smi samples: {len(samples)}")
if sms:
    print(f"sm_mhz: min={min(sms):.0f} median={statistics.median(sms):.0f} max={max(sms):.0f} p90={sorted(sms)[int(0.9*len(sms))-1]:.0f}")
if powers:
    print(f"power_w: min={min(powers):.0f} median={statistics.median(powers):.0f} max={max(powers):.0f}")

def rate(path):
    text = path.read_text(errors="replace")
    m = re.findall(r"finished:\s+([0-9.]+)\s+M it/s", text)
    return [float(x)/1000.0 for x in m]

print("bench-auto B/s:", rate(root/"bench-auto.log"))
print("bench-385k B/s:", rate(root/"bench-385k.log"))
# device line
for p in (root/"bench-auto.log", root/"bench-385k.log"):
    for line in p.read_text(errors="replace").splitlines():
        if "device:" in line or "registers" in line and "packed kernel" in line:
            print(p.name+":", line)
PY
