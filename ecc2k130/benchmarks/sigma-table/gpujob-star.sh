#!/bin/bash
# One-knob star on top of the Frobenius-table build (SIGMA-TABLE.md): the
# table at B16/T512/min1 is the baseline; each arm changes one knob.  Gates
# (arithmetic, integration, replay, sorted-corpus identity) run for every
# binary; then one warm-up per binary and three alternating pairs per arm at
# 32 launches.  Container contract of modal_job.py.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
R=${RESULTS:-/results}
mkdir -p "$R"
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCH="-gencode arch=compute_${CAP},code=sm_${CAP}"
VERIFY_THREADS=96256
BENCH_THREADS=385024
LAUNCHES=32
PAIRS=${PAIRS:-3}
fail=0
{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
} | tee "$R/host.txt"
sha256sum Makefile include/packedkernels.cuh include/packedengine.cuh include/packed131.h \
  src/main.cu src/testpackedcuda.cu benchmarks/sigma-table/gpujob-star.sh > "$R/source-files.sha256"

COMMON=(
  PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
  PACKED_INLINE_POLY=3 SIGMA_FUSED=1 WITNESS=0 WALK_TABLE=0
  THREADS=512 MINBLOCKS=1 PACKED_SIGMA_TABLE=1
)
ARMS="baseline b32 invpoly2 fromreduced"
knobs() {
  case "$1" in
    baseline)    echo "BATCH=16" ;;
    b32)         echo "BATCH=32" ;;
    invpoly2)    echo "BATCH=16 PACKED_INV_POLY=2" ;;
    fromreduced) echo "BATCH=16 PACKED_FROM_REDUCED=1" ;;
  esac
}
batchOf() { case "$1" in b32) echo 32 ;; *) echo 16 ;; esac; }

build() {
  local name=$1
  echo "=== build $name $(knobs "$name")"
  # shellcheck disable=SC2046
  make -s -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" $(knobs "$name") > "$R/build-$name.log" 2>&1 || return 1
  grep -E "registers|spill|stack frame" "$R/build-$name.log" | tail -4 | tee "$R/build-$name.txt" || true
  [ -x ecc2k130 ] || return 1
  mv ecc2k130 "ecc2k130-$name"
  echo "=== arithmetic $name"
  # shellcheck disable=SC2046
  make -s test-packed-cuda ARCH="$ARCH" "${COMMON[@]}" $(knobs "$name") > "$R/arithmetic-$name.log" 2>&1 || return 1
  grep -E "PASS|mismatch" "$R/arithmetic-$name.log" | tee "$R/arithmetic-$name.txt"
  grep -qx "packed arithmetic sigma table: 1" "$R/arithmetic-$name.log" || return 1
  echo "=== integration $name"
  python3 codegen/testpackedclient.py "./ecc2k130-$name" > "$R/integration-$name.log" 2>&1 || return 1
  grep -E "PASS|Error|assert" "$R/integration-$name.log" | tee "$R/integration-$name.txt"
}
for arm in $ARMS; do build "$arm" || fail=1; done
if [ "$fail" != 0 ]; then tail -40 "$R"/build-*.log "$R"/arithmetic-*.log "$R"/integration-*.log 2>/dev/null; exit "$fail"; fi
sha256sum ecc2k130-* > "$R/binary-sha256.txt"

# The replay keeps the walk population fixed (threads x batch = 96256 x 16), so
# every arm seeds the same lanes and must reproduce the same sorted corpus.
verify() {
  local name=$1 status=0 b threads
  b=$(batchOf "$name"); threads=$((VERIFY_THREADS * 16 / b))
  ./ecc2k130-$name --curve 131 --packed --threads "$threads" \
    --dp-weight 48 --dp-cap 262144 --steps 95 --launches 7 \
    --verify 300 --run-id 23 --dp-file "$R/dp-$name.bin" > "$R/verify-$name.log" 2>&1 || status=$?
  grep -E "MISMATCH|OVERFLOW|finished|packed sigma|registers" "$R/verify-$name.log" | tee "$R/verify-$name.txt" || true
  grep -q "^packed sigma table: 1 " "$R/verify-$name.log" || status=1
  grep -qx "packed witness: 0" "$R/verify-$name.log" || status=1
  grep -Eq "^backend cuda-packed131: $threads threads x $b slots x 1 lanes = $((threads*b)) walks, dp weight 48, 95 steps per launch$" "$R/verify-$name.log" || status=1
  grep -Eq "\(300 verified against the reference, 0 dropped\)" "$R/verify-$name.log" || status=1
  ! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-$name.log" || status=1
  return "$status"
}
for arm in $ARMS; do verify "$arm" || fail=1; done
# Same lanes, same seeds, same steps: every arm must reproduce the baseline corpus.
python3 - "$R" <<'PY' | tee "$R/corpus-identity.txt" || fail=1
import hashlib, sys, os
R = sys.argv[1]
def records(name):
    data = open(os.path.join(R, f'dp-{name}.bin'), 'rb').read()
    n = len(data) // 32
    return sorted(data[32 * i:32 * i + 32] for i in range(n))
ref = records('baseline'); ok = len(ref) > 0
print(f'baseline {len(ref)} records sha256 {hashlib.sha256(b"".join(ref)).hexdigest()}')
for name in ('b32', 'invpoly2', 'fromreduced'):
    r = records(name); same = r == ref; ok = ok and same
    print(f'{name} {len(r)} records {"IDENTICAL" if same else "DIFFERENT"}')
print('PASS: sorted corpus identity for every arm' if ok else 'FAIL: corpora differ or empty')
sys.exit(0 if ok else 1)
PY
if [ "$fail" != 0 ]; then echo "correctness preflight failed; timing suppressed" | tee "$R/preflight.txt"; exit "$fail"; fi
echo "PASS: arithmetic, integration, replay and corpus gates" | tee "$R/preflight.txt"

printf 'phase\tarm\tpair\torder\tvariant\trateMps\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 arm=$2 pair=$3 order=$4 variant=$5 log rc count rate state b
  b=$(batchOf "$variant")
  log="$R/${phase}-${arm}-${pair}-${order}-${variant}.log"
  ./ecc2k130-$variant --curve 131 --packed --threads "$BENCH_THREADS" --bench --steps 1024 --launches "$LAUNCHES" --verify 0 > "$log" 2>&1
  rc=$?
  count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  grep -q "^packed sigma table: 1 " "$log" || rc=1
  grep -Eq "^backend cuda-packed131: $BENCH_THREADS threads x $b slots x 1 lanes = $((BENCH_THREADS*b)) walks, dp weight 0, 1024 steps per launch$" "$log" || rc=1
  if [ "$count" != 1 ] || ! awk -v r="$rate" 'BEGIN{exit !(r+0>0)}'; then rc=1; fi
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$phase" "$arm" "$pair" "$order" "$variant" "$rate" "$state" | tee -a "$R/samples.tsv"
  return "$rc"
}
for arm in $ARMS; do sample warmup "$arm" 0 1 "$arm" || fail=1; done
for arm in b32 invpoly2 fromreduced; do
  for pair in $(seq 1 "$PAIRS"); do
    if [ $((pair % 2)) = 1 ]; then
      sample ab "$arm" "$pair" 1 baseline || fail=1; sample ab "$arm" "$pair" 2 "$arm" || fail=1
    else
      sample ab "$arm" "$pair" 1 "$arm" || fail=1; sample ab "$arm" "$pair" 2 baseline || fail=1
    fi
  done
done
if [ "$fail" != 0 ]; then echo "timing row failed" | tee -a "$R/failures.txt"; exit "$fail"; fi
python3 - "$R/samples.tsv" "$R/result.json" <<'PY' | tee "$R/summary.txt"
import csv, json, statistics, sys
rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter='\t') if r['phase'] == 'ab']
out = {}
for arm in ('b32', 'invpoly2', 'fromreduced'):
    pairs = {}
    for r in rows:
        if r['arm'] == arm:
            pairs.setdefault(int(r['pair']), {})[r['variant']] = float(r['rateMps'])
    ratios = [pairs[p][arm] / pairs[p]['baseline'] for p in sorted(pairs)]
    out[arm] = dict(baseline=[pairs[p]['baseline'] for p in sorted(pairs)], candidate=[pairs[p][arm] for p in sorted(pairs)],
                    ratios=ratios, median=statistics.median(ratios), minimum=min(ratios))
    print('%-12s ratios %s median %.6f min %.6f' % (arm, ['%.6f' % x for x in ratios], out[arm]['median'], out[arm]['minimum']))
json.dump(out, open(sys.argv[2], 'w'), indent=1)
PY
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo "=== done"
