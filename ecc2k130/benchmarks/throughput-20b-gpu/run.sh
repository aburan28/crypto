#!/usr/bin/env bash
# 20 B/s attempt on this RTX PRO 6000. Contract:
#   benchmarks/throughput-20b-gpu/README.md
# Run from ecc2k130/:
#   bash benchmarks/throughput-20b-gpu/run.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
OUT_DIR="benchmarks/throughput-20b-gpu"
ALLOW_CONTENTION="${ALLOW_CONTENTION:-0}"

if [ "$ALLOW_CONTENTION" != "1" ]; then
  BUSY="$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader 2>/dev/null)"
  if [ -n "$BUSY" ]; then
    echo "REFUSING: something is already using this GPU:" >&2
    echo "$BUSY" >&2
    exit 2
  fi
fi

if [ -x /usr/local/cuda-13.3/bin/nvcc ]; then
  export PATH="/usr/local/cuda-13.3/bin:$PATH"
  export NVCC="/usr/local/cuda-13.3/bin/nvcc"
  export CUDA_HOME="/usr/local/cuda-13.3"
fi

# Eighteen's geometry and the RTX PRO 6000 arithmetic/storage preset.
export ARCH="-gencode arch=compute_120,code=sm_120"
export BATCH=16 THREADS=256 MINBLOCKS=2 \
       PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 \
       PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 \
       PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 \
       PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 \
       PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 \
       PACKED_CLMAD=1 PACKED_STATE_TILE=256 \
       PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 \
       PACKED_SHARED_SIGMA=1 PACKED_TOP_CLMAD=0 \
       PACKED_PAIR_ILP=0 PACKED_L2_PERSIST=0 \
       UNROLL_SLOTS=1 PACKED_SLOT_PREFETCH=0 \
       WALK_TABLE=0 TABLE_PIVOT_BYTES=0

REPEATS="${REPEATS:-3}"
STEPS="${STEPS:-1024}"
LAUNCHES="${LAUNCHES:-32}"

# Best-effort lock. Failure is recorded, not fatal: the receipt still has
# the clocks.sm samples next to each bench.
sudo -n nvidia-smi -lgc 2430,2430 >/tmp/ecc2k130-lgc.txt 2>&1 || true
sudo -n nvidia-smi -pl 600 >/tmp/ecc2k130-pl.txt 2>&1 || true

python3 - "$OUT_DIR" "$REPEATS" "$STEPS" "$LAUNCHES" <<'PY'
import json, os, re, subprocess, sys, hashlib, pathlib

sys.path.insert(0, ".")
from codegen.benchreport import benchResult, summarizeSamples, reportsVerified

out_dir, repeats, steps, launches = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
out_path = os.path.join(out_dir, "result.json")
arch = os.environ["ARCH"]
preset = (
    f"ARCH='{arch}' BATCH=16 THREADS=256 MINBLOCKS=2 "
    "PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 "
    "PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 "
    "PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1 "
    "PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256 "
    "PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 "
    "PACKED_TOP_CLMAD=0 PACKED_L2_PERSIST=0 UNROLL_SLOTS=1 PACKED_SLOT_PREFETCH=0 TABLE_GLOBAL=0"
)
modes = [
    ("shipping", "WALK_TABLE=0 TABLE_PIVOT_BYTES=0 PACKED_PAIR_ILP=0 PACKED_L2_PERSIST=0 UNROLL_SLOTS=1"),
    ("table_pivot", "WALK_TABLE=1 TABLE_PIVOT_BYTES=1 PACKED_PAIR_ILP=0 PACKED_L2_PERSIST=0 UNROLL_SLOTS=1"),
    ("table_pivot_ilp", "WALK_TABLE=1 TABLE_PIVOT_BYTES=1 PACKED_PAIR_ILP=1 PACKED_L2_PERSIST=0 UNROLL_SLOTS=1"),
    ("table_pivot_persist", "WALK_TABLE=1 TABLE_PIVOT_BYTES=1 PACKED_PAIR_ILP=1 PACKED_L2_PERSIST=1 UNROLL_SLOTS=1"),
    ("table_pivot_unroll", "WALK_TABLE=1 TABLE_PIVOT_BYTES=1 PACKED_PAIR_ILP=1 PACKED_L2_PERSIST=1 UNROLL_SLOTS=2"),
]
rec = {
    "kind": "RTX PRO 6000 20 B/s attempt: complete scalar updates/s",
    "unit": "B scalar updates / s",
    "class": "engineering",
    "targetB": 20.0,
    "floorB": [22.0, 25.0],
    "eighteenReferenceB": 16.666638,
    "method": "codegen.benchreport parseRate/summarizeSamples; median of finished rates",
    "steps": int(steps),
    "launches": int(launches),
    "repeats": repeats,
    "workers": "automatic (multiProcessorCount-scaled)",
    "clockLock": "nvidia-smi -lgc 2430,2430 (best-effort)",
}

def run(key, cmd, cwd=None):
    print(f"--- {key}: {cmd}", flush=True)
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    rec[key] = {"command": cmd, "returncode": p.returncode,
                "stdout": p.stdout, "stderr": p.stderr[-8000:]}
    sys.stdout.write(p.stdout[-3000:])
    if p.returncode:
        sys.stdout.write(p.stderr[-2000:])
    sys.stdout.flush()
    return p

run("nvidiaSmi",
    "nvidia-smi --query-gpu=name,uuid,driver_version,clocks.sm,clocks.max.sm,"
    "power.limit,persistence_mode --format=csv")
run("clockLock", "cat /tmp/ecc2k130-lgc.txt /tmp/ecc2k130-pl.txt")
run("nvcc", "nvcc --version")

host = run("hostTests",
           "mkdir -p build && "
           "g++ -O2 -std=c++17 -DECC_PACKED_SINGLE_PRODUCT=1 -DECC_PACKED_BY_VALUE=1 "
           "-DECC_PACKED_PERM_SIGMA=3 -DECC_PACKED_DIRECT_REDUCE=1 "
           "-DECC_PACKED_GENERATED_PRODUCT=1 -DECC_PACKED_CLMAD=1 "
           "-DECC_PACKED_PAIR_ILP=1 -Wno-unknown-pragmas src/testpacked.cpp "
           "-o build/test-packed-pair-ilp-host && ./build/test-packed-pair-ilp-host && "
           "make test-table-walk-host")
if host.returncode != 0:
    rec["fatal"] = "host tests failed; no rate recorded"
    json.dump(rec, open(out_path, "w"), indent=1)
    sys.exit(1)

binaries = {}
for name, knobs in modes:
    log = os.path.join(out_dir, f"build-{name}.log")
    cmd = f"make -B ecc2k130 {preset} {knobs} >{log} 2>&1 && mv ecc2k130 ecc2k130-{name}"
    p = run(f"build_{name}", cmd)
    rec[f"build_{name}"]["logTail"] = open(log).read()[-4000:]
    if p.returncode != 0:
        rec["fatal"] = f"build {name} failed"
        json.dump(rec, open(out_path, "w"), indent=1)
        sys.exit(1)
    binaries[name] = f"./ecc2k130-{name}"

# Device selection primitives, both pivot layouts, with the ILP binary's flags
# on the last build (PAIR_ILP does not enter those kernels).
dev = run("tableWalkCuda",
          f"make test-table-walk-cuda {preset} WALK_TABLE=1 TABLE_PIVOT_BYTES=1 PACKED_PAIR_ILP=1")
if dev.returncode != 0:
    rec["fatal"] = "device table-walk probe failed"
    json.dump(rec, open(out_path, "w"), indent=1)
    sys.exit(1)

verify_cmd = ("{bin} --curve 131 --packed --dp-weight 50 --steps 16 "
              "--launches 6 --dp-cap 262144 --verify 300")
for name, binp in binaries.items():
    p = run(f"verify_{name}", verify_cmd.format(bin=binp))
    rec[f"verify_{name}"]["reportsOk"] = reportsVerified(
        p.returncode, p.stdout, required=300)
    if not rec[f"verify_{name}"]["reportsOk"]:
        rec["fatal"] = f"verify {name}: need 300 replayed reports, none dropped"
        json.dump(rec, open(out_path, "w"), indent=1)
        sys.exit(1)

BENCH = ("{bin} --curve 131 --packed --bench --steps " + steps +
         " --launches " + launches + " --verify 0")
order = []
for rep in range(repeats):
    # Alternate so thermal drift cannot favour the last binary.
    seq = ["shipping", "table_pivot", "table_pivot_ilp", "table_pivot_persist", "table_pivot_unroll"]
    if rep % 2:
        seq = list(reversed(seq))
    order.extend((rep, n) for n in seq)
rec["benchOrder"] = [n for _, n in order]

samples = {name: [] for name, _ in modes}
for rep, name in order:
    binp = binaries[name]
    cmd = BENCH.format(bin=binp)
    p = run(f"bench_{name}_{rep}", cmd)
    smi = subprocess.run(
        "nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader",
        shell=True, capture_output=True, text=True)
    rec[f"bench_{name}_{rep}"]["smi"] = smi.stdout.strip()
    samples[name].append(benchResult(cmd, p.returncode, p.stdout))

def bps(summary):
    return None if not summary.get("valid") else summary["rate"] / 1000.0

table = {}
any_over = False
for name, _ in modes:
    summary = summarizeSamples(samples[name])
    summary["B_it_per_s"] = bps(summary)
    rec[name] = summary
    table[name] = summary["B_it_per_s"]
    if summary["B_it_per_s"] is not None and summary["B_it_per_s"] > 20.0:
        any_over = True

rec["anyOver20"] = any_over
rec["tableB"] = table
for name, val in table.items():
    if val is None:
        rec[f"ratio_{name}_to_20"] = None
        rec[f"ratio_{name}_to_floor22"] = None
        continue
    rec[f"ratio_{name}_to_20"] = val / 20.0
    rec[f"ratio_{name}_to_floor22"] = val / 22.0
    rec[f"ratio_{name}_to_eighteen"] = val / 16.666638

json.dump(rec, open(out_path, "w"), indent=1)
print(f"\nwrote {out_path}", flush=True)
print(json.dumps({
    "anyOver20": any_over,
    "tableB": table,
    "ratios_to_20": {n: rec[f"ratio_{n}_to_20"] for n, _ in modes},
}, indent=1), flush=True)
if not all(rec[n].get("valid") for n, _ in modes):
    sys.exit(1)
PY
