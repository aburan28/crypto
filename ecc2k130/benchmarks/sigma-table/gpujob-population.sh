#!/bin/bash
# Worker population, L2 persisting window and geometry attribution for the
# Frobenius-table build (SIGMA-TABLE.md), after ONE-BLOCK-GEOMETRY.md: the
# table walk's 20 B/s needed one wave of workers (the state in L2), the L2
# persist window and one block per SM.  Arms: the fused preset at 256x2 and at
# 512x1 (no table), the table, and the table with PACKED_L2_PERSIST=1; each at
# one wave (96,256 workers) and the first and third also at four waves
# (385,024).  Gates for every binary, then rotating rounds at 32 launches.
# Also runs roofline.py on the fused and table builds.  modal_job.py contract.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
R=${RESULTS:-/results}
mkdir -p "$R"
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCH="-gencode arch=compute_${CAP},code=sm_${CAP}"
VERIFY_THREADS=96256
ROUNDS=${ROUNDS:-3}
fail=0
{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
} | tee "$R/host.txt"
sha256sum Makefile include/packedkernels.cuh include/packedengine.cuh include/packed131.h \
  src/main.cu benchmarks/sigma-table/gpujob-population.sh > "$R/source-files.sha256"

COMMON=(
  BATCH=16 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
  PACKED_INLINE_POLY=3 SIGMA_FUSED=1 WITNESS=0 WALK_TABLE=0
)
ARMS="ctrl256 ctrl512 table tablep"
knobs() {
  case "$1" in
    ctrl256) echo "THREADS=256 MINBLOCKS=2 PACKED_SIGMA_TABLE=0 PACKED_L2_PERSIST=0" ;;
    ctrl512) echo "THREADS=512 MINBLOCKS=1 PACKED_SIGMA_TABLE=0 PACKED_L2_PERSIST=0" ;;
    table)   echo "THREADS=512 MINBLOCKS=1 PACKED_SIGMA_TABLE=1 PACKED_L2_PERSIST=0" ;;
    tablep)  echo "THREADS=512 MINBLOCKS=1 PACKED_SIGMA_TABLE=1 PACKED_L2_PERSIST=1" ;;
  esac
}
tableOf() { case "$1" in table|tablep) echo 1 ;; *) echo 0 ;; esac; }
persistOf() { case "$1" in tablep) echo 1 ;; *) echo 0 ;; esac; }

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
  grep -qx "packed arithmetic sigma table: $(tableOf "$name")" "$R/arithmetic-$name.log" || return 1
  echo "=== integration $name"
  python3 codegen/testpackedclient.py "./ecc2k130-$name" > "$R/integration-$name.log" 2>&1 || return 1
  grep -E "PASS|Error|assert" "$R/integration-$name.log" | tee "$R/integration-$name.txt"
}
for arm in $ARMS; do build "$arm" || fail=1; done
if [ "$fail" != 0 ]; then tail -40 "$R"/build-*.log "$R"/arithmetic-*.log "$R"/integration-*.log 2>/dev/null; exit "$fail"; fi
sha256sum ecc2k130-* > "$R/binary-sha256.txt"

echo "=== roofline"
for t in sigma-fused:15.51 sigma-table:15.88; do
  target=${t%%:*}; rate=${t##*:}
  python3 roofline.py --target "gpu-rtx-pro6000-$target" --measured "$rate" --ops 25 --functions 12 \
    --json "$R/roofline-$target.json" > "$R/roofline-$target.txt" 2>&1 || echo "roofline $target failed (non-fatal)" | tee -a "$R/roofline-$target.txt"
  head -40 "$R/roofline-$target.txt"
done

verify() {
  local name=$1 status=0
  ./ecc2k130-$name --curve 131 --packed --threads "$VERIFY_THREADS" \
    --dp-weight 48 --dp-cap 262144 --steps 95 --launches 7 \
    --verify 300 --run-id 23 --dp-file "$R/dp-$name.bin" > "$R/verify-$name.log" 2>&1 || status=$?
  grep -E "MISMATCH|OVERFLOW|finished|packed sigma|packed L2|registers" "$R/verify-$name.log" | tee "$R/verify-$name.txt" || true
  grep -q "^packed sigma table: $(tableOf "$name") " "$R/verify-$name.log" || status=1
  grep -qx "packed L2 persist: $(persistOf "$name")" "$R/verify-$name.log" || status=1
  grep -qx "packed witness: 0" "$R/verify-$name.log" || status=1
  grep -Eq "\(300 verified against the reference, 0 dropped\)" "$R/verify-$name.log" || status=1
  ! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-$name.log" || status=1
  return "$status"
}
for arm in $ARMS; do verify "$arm" || fail=1; done
python3 - "$R" <<'PY' | tee "$R/corpus-identity.txt" || fail=1
import hashlib, sys, os
R = sys.argv[1]
def records(name):
    data = open(os.path.join(R, f'dp-{name}.bin'), 'rb').read()
    return sorted(data[32 * i:32 * i + 32] for i in range(len(data) // 32))
ref = records('ctrl256'); ok = len(ref) > 0
print(f'ctrl256 {len(ref)} records sha256 {hashlib.sha256(b"".join(ref)).hexdigest()}')
for name in ('ctrl512', 'table', 'tablep'):
    r = records(name); same = r == ref; ok = ok and same
    print(f'{name} {len(r)} records {"IDENTICAL" if same else "DIFFERENT"}')
print('PASS: sorted corpus identity for every arm' if ok else 'FAIL: corpora differ or empty')
sys.exit(0 if ok else 1)
PY
if [ "$fail" != 0 ]; then echo "correctness preflight failed; timing suppressed" | tee "$R/preflight.txt"; exit "$fail"; fi
echo "PASS: arithmetic, integration, replay and corpus gates" | tee "$R/preflight.txt"

printf 'phase\tround\tworkers\tvariant\trateMps\tpersistWindow\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 round=$2 workers=$3 variant=$4 log rc count rate state window
  log="$R/${phase}-${round}-${workers}-${variant}.log"
  ./ecc2k130-$variant --curve 131 --packed --threads "$workers" --bench --steps 1024 --launches 32 --verify 0 > "$log" 2>&1
  rc=$?
  count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  window=$(grep -o "^packed L2 persist window: .*" "$log" | head -1 | cut -c26- | tr '\t' ' ')
  grep -q "^packed sigma table: $(tableOf "$variant") " "$log" || rc=1
  grep -Eq "^backend cuda-packed131: $workers threads x 16 slots x 1 lanes = $((workers*16)) walks, dp weight 0, 1024 steps per launch$" "$log" || rc=1
  if [ "$count" != 1 ] || ! awk -v r="$rate" 'BEGIN{exit !(r+0>0)}'; then rc=1; fi
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$phase" "$round" "$workers" "$variant" "$rate" "${window:-none}" "$state" | tee -a "$R/samples.tsv"
  return "$rc"
}
for arm in $ARMS; do sample warmup 0 96256 "$arm" || fail=1; done
for round in $(seq 1 "$ROUNDS"); do
  case $((round % 4)) in
    1) order="ctrl256 ctrl512 table tablep" ;;
    2) order="tablep table ctrl512 ctrl256" ;;
    3) order="table ctrl256 tablep ctrl512" ;;
    0) order="ctrl512 tablep ctrl256 table" ;;
  esac
  for arm in $order; do sample one "$round" 96256 "$arm" || fail=1; done
  if [ "$round" -le 2 ]; then
    if [ "$round" = 1 ]; then sample four "$round" 385024 ctrl256 || fail=1; sample four "$round" 385024 table || fail=1
    else sample four "$round" 385024 table || fail=1; sample four "$round" 385024 ctrl256 || fail=1; fi
  fi
done
if [ "$fail" != 0 ]; then echo "timing row failed" | tee -a "$R/failures.txt"; exit "$fail"; fi
python3 - "$R/samples.tsv" "$R/result.json" <<'PY' | tee "$R/summary.txt"
import csv, json, statistics, sys
rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter='\t') if r['phase'] in ('one', 'four')]
out = {}
for pop in ('96256', '385024'):
    for v in ('ctrl256', 'ctrl512', 'table', 'tablep'):
        rates = [float(r['rateMps']) for r in rows if r['workers'] == pop and r['variant'] == v]
        if rates:
            out[f'{v}@{pop}'] = dict(rates=rates, median=statistics.median(rates))
            print('%-8s @%-6s rates %s median %.3f' % (v, pop, ['%.3f' % x for x in rates], statistics.median(rates)))
def ratio(a, b):
    return out[a]['median'] / out[b]['median'] if a in out and b in out else None
out['ratios'] = {
    'ctrl512/ctrl256 @96256 (geometry alone)': ratio('ctrl512@96256', 'ctrl256@96256'),
    'table/ctrl256 @96256': ratio('table@96256', 'ctrl256@96256'),
    'table/ctrl512 @96256 (table alone)': ratio('table@96256', 'ctrl512@96256'),
    'tablep/table @96256 (L2 persist)': ratio('tablep@96256', 'table@96256'),
    'table@96256 / table@385024 (one wave vs four)': ratio('table@96256', 'table@385024'),
    'ctrl256@96256 / ctrl256@385024': ratio('ctrl256@96256', 'ctrl256@385024'),
}
for k, v in out['ratios'].items():
    print('%-48s %s' % (k, '%.6f' % v if v else 'n/a'))
json.dump(out, open(sys.argv[2], 'w'), indent=1)
PY
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo "=== done"
