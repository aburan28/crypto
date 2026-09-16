#!/usr/bin/env bash
# Measure complete scalar iterations/s for GF(2^131) on an Ada part: the L4
# (EC2 g6) or the L40S (EC2 g6e), both sm_89.
#
# Why this exists: every rate in this tree is for the RTX PRO 6000 Blackwell,
# sm_120 (14.637530 B/s benchmarking, 14.1 B/s collecting -- THROUGHPUT-30B.md),
# and the g6/g6e purchasing decision has no measurement behind it at all. See
# ../../ADA-L4-L40S.md for the break-even this feeds.
#
# It also settles the one preset knob that does NOT transfer off sm_120.
# PACKED_CLMAD buys a carryless op costing ~37.7 ALU ops with a pipe balance
# measured only on Blackwell -- 87.3% ALU against 51.4% carryless. Ada divides
# its units differently, so this runs BOTH arms on the same device allocation,
# alternating, and lets the receipt pick. Nothing here assumes an answer.
#
# Run it ON a g6 or g6e instance, from the ecc2k130 directory:
#   bash benchmarks/ada/run.sh
#
# Writes benchmarks/ada/result.json with every command's verbatim output,
# matching benchmarks/rtx-pro4500/result.json.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
OUT="${OUT:-benchmarks/ada/result.json}"
ALLOW_CONTENTION="${ALLOW_CONTENTION:-0}"

# A collecting worker on the same GPU steals SMs from the benchmark and the
# benchmark steals them back. Both numbers come out low and neither measures
# anything. Refuse rather than record a quietly wrong rate.
if [ "$ALLOW_CONTENTION" != "1" ]; then
  BUSY="$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader 2>/dev/null)"
  if [ -n "$BUSY" ]; then
    echo "REFUSING: something is already using this GPU:" >&2
    echo "$BUSY" >&2
    echo >&2
    echo "A bench sharing the GPU with a collecting worker measures neither." >&2
    echo "Stop the worker first (it checkpoints; see --checkpoint), then re-run." >&2
    echo "Override with ALLOW_CONTENTION=1 only if you intend a contention test," >&2
    echo "in which case the result is NOT comparable to any other figure here." >&2
    exit 2
  fi
fi

# Refuse a part this script does not describe. An A10G or a T4 would build and
# run and produce a number filed under "Ada", which is worse than no number.
NAME="$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -n1)"
case "$NAME" in
  *L4*|*L40S*|*L40*) : ;;
  *) if [ "${ALLOW_OTHER_GPU:-0}" != "1" ]; then
       echo "REFUSING: this GPU is '$NAME', not an L4 or L40S." >&2
       echo "Set ALLOW_OTHER_GPU=1 to record it anyway; the result is then not" >&2
       echo "an Ada g6/g6e figure whatever this file is named." >&2
       exit 2
     fi ;;
esac

# Build only sm_89. The Makefile default is a five-architecture fat binary and
# ptxas takes minutes per architecture on this kernel; four of those passes
# would be for parts this instance does not have.
export ARCH="-gencode arch=compute_89,code=sm_89"

# The RTX PRO 6000 preset's arithmetic and layout options, in the names the
# MAKEFILE consumes. These are NOT the ECC_PACKED_* names: those are the Modal
# image interface (modal_app.py reads them from the environment), while a local
# `make` build reads BATCH/THREADS/MINBLOCKS/PACKED_* and turns them into the
# -DECC_PACKED_* defines itself. Using the Modal spelling here would leave every
# option at its Makefile default -- batch 32, all PACKED_* zero -- and quietly
# measure a different kernel than the figures this compares to.
#
# PACKED_CLMAD is deliberately absent: it is the variable, set per arm below.
COMMON="BATCH=16 THREADS=256 MINBLOCKS=2 \
PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 \
PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 \
PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 \
PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 \
PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 \
PACKED_STATE_TILE=256 PACKED_WEIGHTED_PREFIX=2 \
PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 PACKED_TOP_CLMAD=0"
export COMMON

REPEATS="${REPEATS:-3}"

python3 - "$OUT" "$REPEATS" "$NAME" <<'PY'
import json, os, subprocess, sys

# The tree's own rate reader and summariser, not a local reimplementation.
# parseRate accepts ONLY the single `finished:` line -- printed after the
# pending GPU work is synchronised -- and rejects a run with none. The periodic
# progress lines also carry "M it/s" and read high during boost-clock warmup;
# the L4 is a 72 W part whose sustained clock is well under its boost, so
# scraping those and taking a maximum would overstate it by more than the
# g6-vs-g6e break-even is wide. summarizeSamples takes the MEDIAN across
# repeats, as the RTX PRO 6000 figures do.
sys.path.insert(0, ".")
from codegen.benchreport import benchResult, summarizeSamples

out_path, repeats, gpuName = sys.argv[1], int(sys.argv[2]), sys.argv[3]
common = os.environ["COMMON"]
rec = {"kind": "Ada (EC2 g6 / g6e) complete-scalar-iteration throughput",
       "gpu": gpuName,
       "method": "codegen.benchreport parseRate/summarizeSamples, identical to "
                 "the RTX PRO 6000 figures; median of finished rates",
       "workers": "automatic (multiProcessorCount-scaled, occupancy-limited)",
       "comparison": {"rtx_pro_6000_bench_B_per_s": 14.637530,
                      "rtx_pro_6000_collecting_B_per_s": 14.1},
       "note": "--bench excludes distinguished-point handling; the collecting "
               "rate is lower. The two arms differ ONLY in PACKED_CLMAD."}

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
ARMS = [("software", "0"), ("clmad", "1")]

rec["arms"] = {}
for arm, clmad in ARMS:
    knobs = f"{common} PACKED_CLMAD={clmad}"
    a = rec["arms"][arm] = {"packedClmad": int(clmad)}
    # Record the defines actually compiled in, so a future reader can confirm
    # the arm reached the build rather than trusting that it did.
    d = run(f"{arm}_buildDefines",
            f"make -n gpu {knobs} | tr ' ' '\\n' | grep -E '^-DECC_' | sort -u")
    a["defines"] = sorted(d.stdout.split())
    if run(f"{arm}_build", f"make -B gpu -j4 {knobs}").returncode != 0:
        # clmad is PTX 9.3 / CUDA 13.3; an older toolkit or a target that
        # rejects it fails HERE, and that is a result, not a crash. The other
        # arm still has a rate and the JSON still gets written.
        a["built"] = False
        a["summary"] = {"valid": False, "B_it_per_s": None,
                        "reason": "build failed; see the arm's build output"}
        continue
    a["built"] = True
    samples = []
    for i in range(repeats):
        p = run(f"{arm}_bench_{i}", BENCH)
        samples.append(benchResult(BENCH, p.returncode, p.stdout))
    summary = summarizeSamples(samples)
    # Never a number from an invalid run.
    summary["B_it_per_s"] = summary["rate"] / 1000.0 if summary["valid"] else None
    a["summary"] = summary

sw = rec["arms"]["software"]["summary"].get("B_it_per_s")
cl = rec["arms"]["clmad"]["summary"].get("B_it_per_s")
if sw and cl:
    rec["clmadGain"] = cl / sw - 1.0
    rec["preset"] = "clmad" if cl > sw else "software"
else:
    rec["clmadGain"] = None
    rec["preset"] = "software" if sw else None

json.dump(rec, open(out_path, "w"), indent=1)
print(f"\nwrote {out_path}")
for arm, _ in ARMS:
    s = rec["arms"][arm]["summary"]
    print(f"  {arm:9s} {s.get('B_it_per_s')} B it/s  (valid={s.get('valid')})")
print(f"  clmad gain: {rec['clmadGain']}   preset: {rec['preset']}")
PY
