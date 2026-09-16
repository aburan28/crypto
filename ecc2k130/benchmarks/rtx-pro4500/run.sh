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
# output, matching the convention of benchmarks/hardware-limits/result.json.
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

# Same arithmetic/layout options as the RTX PRO 6000 preset (Makefile
# RTX_PRO6000_ENV). These are field-arithmetic and storage choices, not
# GPU-model choices, so they transfer; the WORKER COUNT deliberately does not
# -- 385,024 was tuned as 4x automatic on a 188-SM part and has no claim to be
# right on a smaller one. Hence the sweep below.
export ECC_PACKED_SINGLE_PRODUCT=1 ECC_PACKED_CACHE_DENOM=1 \
       ECC_PACKED_BY_VALUE=1 ECC_PACKED_PERM_SIGMA=3 \
       ECC_PACKED_POLY_CHAIN=1 ECC_PACKED_UNROLL_INV=1 \
       ECC_PACKED_PAIR_PRODUCTS=1 ECC_PACKED_POLY_STATE=1 \
       ECC_PACKED_DIRECT_REDUCE=1 ECC_PACKED_GENERATED_PRODUCT=1 \
       ECC_PACKED_CLMAD=1 ECC_PACKED_STATE_TILE=256 \
       ECC_PACKED_WEIGHTED_PREFIX=2 ECC_PACKED_COMPACT_STATE=1 \
       ECC_PACKED_SHARED_SIGMA=1 ECC_PACKED_TOP_CLMAD=0 \
       ECC_BATCH=16 ECC_THREADS=256

python3 - "$OUT" <<'PY'
import json, subprocess, sys, os, re

out_path = sys.argv[1]
rec = {"kind": "RTX PRO 4500 (EC2 g7) complete-scalar-iteration throughput",
       "note": "Rates are from --bench (no distinguished-point handling); the "
               "collecting rate is lower. The RTX PRO 6000 comparison figures "
               "are 14.637530 B/s benchmarking and 14.1 B/s collecting."}

def run(key, cmd, **kw):
    print(f"--- {key}: {cmd}", flush=True)
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)
    rec[key] = {"command": cmd, "returncode": p.returncode,
                "stdout": p.stdout, "stderr": p.stderr[-4000:]}
    print(p.stdout[-2000:], flush=True)
    return p

run("nvidiaSmi", "nvidia-smi --query-gpu=name,uuid,driver_version,"
                 "clocks.max.sm,clocks.max.memory,power.limit --format=csv")
run("smCount", "nvidia-smi --query-gpu=name,count --format=csv,noheader")
run("nvcc", "nvcc --version")
build = run("build", "make gpu -j4")
if build.returncode != 0:
    rec["fatal"] = "build failed; no rate recorded"
    json.dump(rec, open(out_path, "w"), indent=1)
    sys.exit(1)

# Automatic workers first (scales with multiProcessorCount), then multiples of
# it, because the 6000's 4x may or may not be this part's optimum.
rec["runs"] = []
for workers in (0, 0, 0):        # three repeats at automatic
    p = run(f"bench_auto_{len(rec['runs'])}",
            "./ecc2k130 --curve 131 --packed --bench --steps 1024 "
            "--launches 32 --verify 0")
    m = re.findall(r"([0-9.]+)\s*M it/s", p.stdout)
    rec["runs"].append({"workers": "automatic",
                        "m_it_per_s": [float(x) for x in m] or None,
                        "returncode": p.returncode})

rates = [r for run_ in rec["runs"] for r in (run_["m_it_per_s"] or [])]
if rates:
    best = max(rates)
    rec["summary"] = {
        "best_M_it_per_s": best,
        "best_B_it_per_s": best / 1000.0,
        "all_M_it_per_s": rates,
        "median_M_it_per_s": sorted(rates)[len(rates) // 2],
    }
else:
    rec["summary"] = {"error": "no rate parsed from output; nothing recorded"}

json.dump(rec, open(out_path, "w"), indent=1)
print(f"\nwrote {out_path}")
print(json.dumps(rec.get("summary"), indent=1))
PY
