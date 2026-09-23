#!/bin/bash
# Energy per update of the RTX PRO 6000 builds, as the predictor of their rate
# on a power-capped part of the same SM (the RTX PRO 4500: sm_120, 82 SMs,
# 165 W; RTX-PRO4500.md, POWER-BOUND.md).  For each binary: 300 device reports
# re-walked, then REPS alternating benches with the board power sampled every
# 100 ms *during* the run; energy per update is mean power over the timed
# window divided by the rate.  If the container may set the power limit, the
# benches are repeated with the 6000 capped to the 4500's power per SM
# (165 W x 188/82 = 378 W, or the card's minimum if that is higher).
#   modal run --detach modal_job.py --job benchmarks/power-bound/gpujob.sh --out DIR --gpu RTX-PRO-6000 [--env SET=fused]
#   python3 aws/bench_job.py --job benchmarks/power-bound/gpujob.sh --out DIR --region us-west-2 \
#       --instance-type g7.2xlarge --env SET=rtx4500,REPS=5     # on a real RTX PRO 4500
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq make g++ python3 >/dev/null 2>&1
R=${RESULTS:-/results}; mkdir -p "$R"
REPS=${REPS:-3}; CAP_W=${CAP_W:-378}
{
  nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,power.limit,power.min_limit,power.max_limit,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-unknown}"
} | tee "$R/host.txt"
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCH="-gencode arch=compute_${CAP},code=sm_${CAP}"
# The audited shipping/table preset of benchmarks/table-walk/gpujob.sh: the geometry the 4500 was measured at.
K="BATCH=16 THREADS=256 MINBLOCKS=2 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256 PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1"
fail=0
build() {  # build NAME make-args...
  local name=$1; shift
  echo "=== build $name: $*"
  make -s "$@" > "$R/build-$name.log" 2>&1
  grep -A2 "eccPacked131.*walk" "$R/build-$name.log" | grep -E "registers|spill" | head -2 | tee "$R/build-$name.txt"
  if [ ! -x ecc2k130 ]; then echo "BUILD FAILED $name" | tee -a "$R/failures.txt"; return 1; fi
  mv ecc2k130 "ecc2k130-$name"
}
PRE="PRO6000_ARCH=$ARCH"
# SET=survey: the tree's builds at both geometries; SET=fused: the fused pass
# of the 20 B/s build against it, with the distinguished-point identity check.
SET=${SET:-survey}
if [ "$SET" = fused ]; then
  build ref20b             gpu-rtx-pro6000-20b "$PRE" || fail=1
  build ref20b-fused       gpu-preset "$PRE" KNOBS="TABLE_FUSED=1 TABLE_PIPE_SELECT=0" || fail=1
  build ref20b-fused-pipe  gpu-preset "$PRE" KNOBS="TABLE_FUSED=1 TABLE_FUSED_PIPE=1 TABLE_PIPE_SELECT=0" || fail=1
  build ref20b-fused-384x1 gpu-preset "$PRE" KNOBS="TABLE_FUSED=1 TABLE_PIPE_SELECT=0 THREADS=384 MINBLOCKS=1" || fail=1
  build ref20b-fused-sqclmad gpu-preset "$PRE" KNOBS="TABLE_FUSED=1 TABLE_PIPE_SELECT=0 PACKED_ALU_SQUARE=0" || fail=1
  VARIANTS="ref20b ref20b-fused ref20b-fused-pipe ref20b-fused-384x1 ref20b-fused-sqclmad"
  FORCE=96256   # the 512-thread builds' automatic count: 1,540,096 walks for every binary
elif [ "$SET" = rtx4500 ]; then
  # On a real RTX PRO 4500 (g7.2xlarge): the two receipts' builds as controls,
  # the 20 B/s build, and the two fused builds POWER-BOUND.md predicts for it.
  build ship-256x2          -B ecc2k130 "ARCH=$ARCH" $K WALK_TABLE=0 WITNESS=0 || fail=1
  build table-256x2         -B ecc2k130 "ARCH=$ARCH" $K WALK_TABLE=1 || fail=1
  build ref20b              gpu-rtx-pro6000-20b "$PRE" || fail=1
  build rtx4500             gpu-rtx-pro4500 "$PRE" || fail=1
  build rtx4500-sqclmad     gpu-rtx-pro4500 "$PRE" RTX4500_EXTRA="PACKED_ALU_SQUARE=0" || fail=1
  VARIANTS="ship-256x2 table-256x2 ref20b rtx4500 rtx4500-sqclmad"
  SMS=$(nvidia-smi --query-gpu=name --format=csv,noheader >/dev/null; ./ecc2k130-ref20b --curve 131 --packed --bench --steps 1 --launches 1 --verify 0 2>&1 | grep -oE '[0-9]+ SMs' | grep -oE '[0-9]+' | head -1)
  FORCE=$((SMS * 512))
else
  build ship-256x2   -B ecc2k130 "ARCH=$ARCH" $K WALK_TABLE=0 WITNESS=0 || fail=1
  build table-256x2  -B ecc2k130 "ARCH=$ARCH" $K WALK_TABLE=1 || fail=1
  build ref20b       gpu-rtx-pro6000-20b "$PRE" || fail=1
  build ref20b-384x1 gpu-preset "$PRE" KNOBS="THREADS=384 MINBLOCKS=1" || fail=1
  build ref20b-640x1 gpu-preset "$PRE" KNOBS="THREADS=640 MINBLOCKS=1" || fail=1
  build ref20b-fused gpu-preset "$PRE" KNOBS="TABLE_FUSED=1 TABLE_PIPE_SELECT=0" || fail=1
  build ref20b-alusqr gpu-preset "$PRE" KNOBS="PACKED_ALU_SQR=1" || fail=1
  build c2-256x32    gpu-rtx-pro6000-chains2 "$PRE" || fail=1
  VARIANTS="ship-256x2 table-256x2 ref20b ref20b-384x1 ref20b-640x1 ref20b-fused ref20b-alusqr c2-256x32"
  FORCE=""
fi

for b in $VARIANTS; do
  [ -x "ecc2k130-$b" ] || continue
  echo "=== verify $b"
  ./ecc2k130-$b --curve 131 --packed ${FORCE:+--threads $FORCE --dp-file $R/dp-$b.bin} --dp-weight 48 --dp-cap 262144 --steps 96 --launches 6 --verify 300 --run-id 7 \
      > "$R/verify-$b.log" 2>&1 || { echo "VERIFY FAILED $b" | tee -a "$R/failures.txt"; fail=1; }
  grep -E "MISMATCH|finished|resident|registers|OVERFLOW|packed (table walk|chains|alu onb square|alu square)" "$R/verify-$b.log" | tee "$R/verify-$b.txt" || true
  grep -q "300 verified against the reference, 0 dropped" "$R/verify-$b.log" || { echo "NOT 300/300: $b" | tee -a "$R/failures.txt"; fail=1; }
done
if [ -n "$FORCE" ]; then
# The shipping walk is a different walk, so it is held to 300/300 but not to the table walk's points.
python3 - "$R" $(echo $VARIANTS | tr ' ' '\n' | grep -v '^ship-') <<'PY' | tee "$R/dp-identity.txt" || fail=1
import hashlib, os, sys
r = sys.argv[1]; ref = None; bad = 0
for name in sys.argv[2:]:
    path = os.path.join(r, "dp-%s.bin" % name)
    if not os.path.exists(path): print("%-22s missing" % name); bad = 1; continue
    data = open(path, "rb").read()
    recs = sorted(data[i:i+32] for i in range(0, len(data) - len(data) % 32, 32))
    digest = hashlib.sha256(b"".join(recs)).hexdigest()
    if ref is None: ref = digest
    if digest != ref: bad = 1
    print("%-22s %8d records  sha256 %s  %s" % (name, len(recs), digest[:16], "IDENTICAL to ref" if digest == ref else "DIFFERS from ref"))
sys.exit(bad)
PY
fi

# One timed run with the board sampled every 100 ms while it runs.
timed() {  # timed NAME TAG
  local b=$1 tag=$2 log="$R/power-$2-$1-$3.csv"
  nvidia-smi --query-gpu=timestamp,power.draw,clocks.sm,temperature.gpu --format=csv,noheader,nounits -lms 100 > "$log" 2>/dev/null &
  local smi=$!
  local t0; t0=$(date +%s.%N)
  local out; out=$(./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 64 --verify 0 2>&1 | grep -E "finished" )
  local t1; t1=$(date +%s.%N)
  kill $smi 2>/dev/null; wait $smi 2>/dev/null
  echo "$tag $b rep $3: $out window $t0 $t1" | tee -a "$R/bench.txt"
}
for rep in $(seq 1 "$REPS"); do for b in $VARIANTS; do [ -x "ecc2k130-$b" ] && timed "$b" uncapped "$rep"; done; done

# The power-capped pass, if the container may set the limit.
[ "${SKIP_CAP:-0}" = 1 ] && { echo "=== done"; exit "$fail"; }
echo "=== power limit: try ${CAP_W} W" | tee -a "$R/powercap.txt"
MINW=$(nvidia-smi --query-gpu=power.min_limit --format=csv,noheader,nounits | head -1 | cut -d. -f1)
[ -n "$MINW" ] && [ "$MINW" -gt "$CAP_W" ] && CAP_W=$MINW
if nvidia-smi -pl "$CAP_W" >> "$R/powercap.txt" 2>&1; then
  nvidia-smi --query-gpu=power.limit --format=csv,noheader | tee -a "$R/powercap.txt"
  for rep in $(seq 1 "$REPS"); do for b in $VARIANTS; do [ -x "ecc2k130-$b" ] && timed "$b" capped "$rep"; done; done
  nvidia-smi -pl "$(nvidia-smi --query-gpu=power.default_limit --format=csv,noheader,nounits | head -1 | cut -d. -f1)" >> "$R/powercap.txt" 2>&1
else
  echo "power limit not settable in this container; capped pass skipped" | tee -a "$R/powercap.txt"
fi
echo "=== done"
exit "$fail"
