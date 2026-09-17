#!/usr/bin/env bash
# Measure complete scalar iterations/s for GF(2^131) on a T4 (EC2 g4dn), sm_75.
#
# The Ada script (../ada/run.sh) refuses a T4 so an L4 number cannot be
# filed under the wrong part. This is the matching T4 script. PACKED_CLMAD
# stays 0: see ../../T4-G4DN.md.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
OUT="${OUT:-benchmarks/t4/result.json}"
ALLOW_CONTENTION="${ALLOW_CONTENTION:-0}"

if [ "$ALLOW_CONTENTION" != "1" ]; then
  BUSY="$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader 2>/dev/null)"
  if [ -n "$BUSY" ]; then
    echo "REFUSING: something is already using this GPU:" >&2
    echo "$BUSY" >&2
    exit 2
  fi
fi

NAME="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n1)"
case "$NAME" in
  *T4*) : ;;
  *) if [ "${ALLOW_OTHER_GPU:-0}" != "1" ]; then
       echo "REFUSING: this GPU is '$NAME', not a T4." >&2
       exit 2
     fi ;;
esac

export ARCH="-gencode arch=compute_75,code=sm_75"
COMMON="BATCH=16 THREADS=256 MINBLOCKS=2 \
PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 \
PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 \
PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 \
PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 \
PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 \
PACKED_STATE_TILE=256 PACKED_WEIGHTED_PREFIX=2 \
PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 \
PACKED_TOP_CLMAD=0 PACKED_CLMAD=0"
export COMMON
REPEATS="${REPEATS:-3}"

python3 - "$OUT" "$REPEATS" "$NAME" <<'PY'
import json, os, subprocess, sys
sys.path.insert(0, ".")
from codegen.benchreport import benchResult, summarizeSamples

out_path, repeats, gpuName = sys.argv[1], int(sys.argv[2]), sys.argv[3]
common = os.environ["COMMON"]
rec = {"kind": "Turing (EC2 g4dn / T4) complete-scalar-iteration throughput",
       "gpu": gpuName,
       "method": "codegen.benchreport parseRate/summarizeSamples; median of finished rates",
       "workers": "automatic",
       "comparison": {"rtx_pro_6000_bench_B_per_s": 15.115792},
       "preset": "software"}

def run(key, cmd):
    print(f"--- {key}: {cmd}", flush=True)
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    rec[key] = {"command": cmd, "returncode": p.returncode,
                "stdout": p.stdout, "stderr": p.stderr[-4000:]}
    print(p.stdout[-2000:], flush=True)
    return p

run("nvidiaSmi", "nvidia-smi --query-gpu=name,uuid,driver_version,"
                 "clocks.max.sm,clocks.max.memory,power.limit --format=csv")
run("nvcc", "nvcc --version")
BENCH = ("./ecc2k130 --curve 131 --packed --bench --steps 1024 "
         "--launches 32 --verify 0")
if run("build", f"make -B gpu -j4 {common}").returncode != 0:
    rec["summary"] = {"valid": False, "B_it_per_s": None, "reason": "build failed"}
else:
    samples = []
    for i in range(repeats):
        p = run(f"bench_{i}", BENCH)
        samples.append(benchResult(BENCH, p.returncode, p.stdout))
    summary = summarizeSamples(samples)
    summary["B_it_per_s"] = summary["rate"] / 1000.0 if summary["valid"] else None
    rec["summary"] = summary
json.dump(rec, open(out_path, "w"), indent=1)
print(f"\nwrote {out_path}")
print(rec["summary"])
PY
