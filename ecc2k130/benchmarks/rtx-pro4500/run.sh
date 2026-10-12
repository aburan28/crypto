#!/usr/bin/env bash
# Measure complete scalar iterations/s for GF(2^131) on an RTX PRO 4500 (EC2 g7).
#
# Why this exists: every throughput figure in this tree is for the RTX PRO 6000
# (14.637530 B/s benchmarking, 14.1 B/s collecting -- THROUGHPUT-30B.md). The
# g7 fleet runs the 4500 and has NO recorded rate, so the g7-vs-g7e purchasing
# decision currently rests on recollection. See ../../RTX-PRO4500.md for the
# break-even this feeds.
#
# Run it ON a g7 instance, from the ecc2k130 directory:
#   bash benchmarks/rtx-pro4500/run.sh
#
# It writes benchmarks/rtx-pro4500/result.json with every command's verbatim
# output, matching benchmarks/hardware-limits/result.json.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
OUT="benchmarks/rtx-pro4500/result.json"
ALLOW_CONTENTION="${ALLOW_CONTENTION:-0}"

# A collecting worker on the same GPU steals SMs from the benchmark and the
# benchmark steals them back. Both numbers come out low and neither is a
# measurement of anything. Refuse rather than record a quietly wrong rate.
if [ "$ALLOW_CONTENTION" != "1" ]; then
  BUSY="$(nvidia-smi --query-compute-apps=pid,process_name --format=csv,noheader 2>/dev/null)"
  if [ -n "$BUSY" ]; then
    echo "REFUSING: something is already using this GPU:" >&2
    echo "$BUSY" >&2
    echo >&2
    echo "A bench sharing the GPU with a collecting worker measures neither." >&2
    echo "Stop the worker first (it checkpoints; see --checkpoint), then re-run." >&2
    echo "Override with ALLOW_CONTENTION=1 only if you intend a contention test," >&2
    echo "in which case the result is NOT comparable to the RTX PRO 6000 figure." >&2
    exit 2
  fi
fi

# The RTX PRO 6000 preset's arithmetic and layout options, in the names the
# MAKEFILE consumes. These are NOT the ECC_PACKED_* names: those are the Modal
# image interface (modal_app.py reads them from the environment), while a local
# `make` build reads BATCH/THREADS/MINBLOCKS/PACKED_* and turns them into the
# -DECC_PACKED_* defines itself. Using the Modal spelling here would leave every
# option at its Makefile default -- batch 32, all PACKED_* zero -- and quietly
# measure a different kernel than the 14.6/14.1 B/s figures this compares to.
#
# These options are field-arithmetic and storage choices, so they transfer
# across GPU models. The 6000 preset's 385,024 WORKERS deliberately do not:
# that is 4x automatic on a 188-SM part and has no claim on a smaller one.
# Automatic scales with multiProcessorCount (src/main.cu), so start there.
export BATCH=16 THREADS=256 MINBLOCKS=2 \
       PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 \
       PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 \
       PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 \
       PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 \
       PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 \
       PACKED_CLMAD=1 PACKED_STATE_TILE=256 \
       PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 \
       PACKED_SHARED_SIGMA=1 PACKED_TOP_CLMAD=0

REPEATS="${REPEATS:-3}"

python3 - "$OUT" "$REPEATS" <<'PY'
import json, subprocess, sys

# The tree's own rate reader and summariser, not a local reimplementation.
# parseRate accepts ONLY the single `finished:` line -- which is printed after
# the pending GPU work is synchronised -- and rejects a run with none. The
# periodic progress lines also carry "M it/s" and read high during boost-clock
# warmup, so scraping them and taking a maximum would sit the 4500 above a
# same-method comparison and could flip a break-even only a few percent wide.
# summarizeSamples takes the MEDIAN across repeats, as the 6000 figure does.
sys.path.insert(0, ".")
from codegen.benchreport import benchResult, summarizeSamples

out_path, repeats = sys.argv[1], int(sys.argv[2])
rec = {"kind": "RTX PRO 4500 (EC2 g7) complete-scalar-iteration throughput",
       "method": "codegen.benchreport parseRate/summarizeSamples, identical to "
                 "the RTX PRO 6000 figures; median of finished rates",
       "comparison": {"rtx_pro_6000_bench_B_per_s": 14.637530,
                      "rtx_pro_6000_collecting_B_per_s": 14.1},
       "note": "--bench excludes distinguished-point handling; the collecting "
               "rate is lower."}

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
# Record the defines actually compiled in, so a future reader can confirm the
# preset reached the build rather than trusting this comment.
run("buildDefines", "make -n gpu | tr ' ' '\\n' | grep -E '^-DECC_' | sort -u")
if run("build", "make gpu -j4").returncode != 0:
    rec["fatal"] = "build failed; no rate recorded"
    json.dump(rec, open(out_path, "w"), indent=1)
    sys.exit(1)

BENCH = ("./ecc2k130 --curve 131 --packed --bench --steps 1024 "
         "--launches 32 --verify 0")
samples = []
for i in range(repeats):
    p = run(f"bench_{i}", BENCH)
    samples.append(benchResult(BENCH, p.returncode, p.stdout))

summary = summarizeSamples(samples)
rec["workers"] = "automatic (multiProcessorCount-scaled)"
rec["summary"] = summary
if summary["valid"]:
    rec["summary"]["B_it_per_s"] = summary["rate"] / 1000.0
else:
    rec["summary"]["B_it_per_s"] = None   # never a number from an invalid run

json.dump(rec, open(out_path, "w"), indent=1)
print(f"\nwrote {out_path}")
print(json.dumps({k: v for k, v in rec["summary"].items() if k != "samples"},
                 indent=1))
PY
