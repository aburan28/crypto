#!/usr/bin/env bash
# Paired WITNESS=0 vs WITNESS=1 on the campaign packed walk.
# Hypothesis: WITNESS=1 (default since cairn) explains 9.3 vs survey 14.5 B/s
# (~36% drop matches CAIRN-WITNESS.md's predicted hot-state increase).
set -euo pipefail
cd /work
RESULTS="${RESULTS:-/results}"
mkdir -p "$RESULTS"

common=(
  BATCH=16 THREADS=256 LEAF=0 MINBLOCKS=2 PACKED=1
  PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
  PACKED_TOP_CLMAD=0 PACKED_CLMAD_SQUARE=0 PACKED_KARAT3=0
  WALK_TABLE=0
)

nvidia-smi --query-gpu=name,driver_version,power.limit --format=csv | tee "$RESULTS/gpu.csv"

run_one() {
  local tag="$1" witness="$2"
  echo "=== build WITNESS=$witness ($tag) ==="
  make -B gpu "${common[@]}" WITNESS="$witness" 2>&1 | tee "$RESULTS/build-$tag.log" | tail -20
  # Confirm the binary identity line
  strings ecc2k130 | grep -E 'witness|WITNESS' | head -5 || true
  echo "=== bench auto workers WITNESS=$witness ==="
  ./ecc2k130 --curve 131 --bench --steps 1024 --launches 32 --verify 0 --packed --threads 0 \
    2>&1 | tee "$RESULTS/bench-$tag.log"
}

# Control first (survey-equivalent), then treatment (fleet default)
run_one control 0
run_one witness 1

python3 - <<'PY' | tee "$RESULTS/summary.txt"
import re, pathlib, statistics
root = pathlib.Path("/results")

def rates(tag):
    text = (root / f"bench-{tag}.log").read_text(errors="replace")
    finished = [float(x)/1000 for x in re.findall(r"finished:\s+([0-9.]+)\s+M it/s", text)]
    # also grab register line
    regs = re.search(r"packed kernel:\s+(\d+) registers", text)
    device = next((ln for ln in text.splitlines() if ln.startswith("device:")), "")
    return finished, int(regs.group(1)) if regs else None, device

c, cr, cd = rates("control")
w, wr, wd = rates("witness")
print("device control:", cd)
print("device witness:", wd)
print("regs control/witness:", cr, wr)
print("control B/s samples:", ["%.3f" % x for x in c])
print("witness B/s samples:", ["%.3f" % x for x in w])
if c and w:
    cm, wm = statistics.median(c), statistics.median(w)
    print(f"median control={cm:.3f} B/s  witness={wm:.3f} B/s  ratio={wm/cm:.4f}  drop={(1-wm/cm)*100:.1f}%")
    print(f"survey reference=14.470 B/s  control/survey={cm/14.470:.4f}  witness/survey={wm/14.470:.4f}")
    if cm >= 13.5 and wm < 11.0:
        print("VERDICT: WITNESS=1 is the cause of the fleet regression")
    elif cm < 11.0 and wm < 11.0:
        print("VERDICT: both slow — not just WITNESS; look at clocks/SKU")
    else:
        print("VERDICT: inconclusive — inspect samples")
PY
