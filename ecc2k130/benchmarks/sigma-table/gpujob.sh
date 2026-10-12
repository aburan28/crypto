#!/bin/bash
# Frobenius nibble tables (PACKED_SIGMA_TABLE) against the confirmed fused sigma
# preset on one RTX PRO 6000: build both, gate on the GPU arithmetic suite, the
# client integration test, a 300-report replay and sorted-corpus identity, then
# time five alternating A/B pairs after two warmups per binary.  Container
# contract of modal_job.py: tree at /work, results under $RESULTS.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
R=${RESULTS:-/results}
mkdir -p "$R"
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCH="-gencode arch=compute_${CAP},code=sm_${CAP}"
VERIFY_THREADS=96256
BENCH_THREADS=${BENCH_THREADS:-385024}
LAUNCHES=${LAUNCHES:-32}
PAIRS=${PAIRS:-5}
fail=0
{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
  echo "verify threads: $VERIFY_THREADS; benchmark threads: $BENCH_THREADS; launches: $LAUNCHES; pairs: $PAIRS"
} | tee "$R/host.txt"
sha256sum Makefile include/packedkernels.cuh include/packedengine.cuh include/packed131.h \
  include/packedsigma131.h src/main.cu src/testpackedcuda.cu benchmarks/sigma-table/gpujob.sh \
  > "$R/source-files.sha256"

COMMON=(
  BATCH=16 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
  PACKED_INLINE_POLY=3 SIGMA_FUSED=1 WITNESS=0 WALK_TABLE=0
)
knobs() {  # name -> geometry and table knobs
  case "$1" in
    control)   echo "THREADS=256 MINBLOCKS=2 PACKED_SIGMA_TABLE=0" ;;
    candidate) echo "THREADS=512 MINBLOCKS=1 PACKED_SIGMA_TABLE=1" ;;
  esac
}
tableOf() { case "$1" in control) echo 0 ;; candidate) echo 1 ;; esac; }

build() {
  local name=$1
  echo "=== build $name $(knobs "$name")"
  # shellcheck disable=SC2046
  make -s -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" $(knobs "$name") > "$R/build-$name.log" 2>&1 || return 1
  grep -E "Function properties for.*walk|registers|spill|stack frame" "$R/build-$name.log" | tail -12 | tee "$R/build-$name.txt" || true
  [ -x ecc2k130 ] || return 1
  mv ecc2k130 "ecc2k130-$name"
  cuobjdump -sass "ecc2k130-$name" > "$R/sass-$name.txt" 2>/dev/null || true
  echo "=== arithmetic $name"
  # shellcheck disable=SC2046
  make -s test-packed-cuda ARCH="$ARCH" "${COMMON[@]}" $(knobs "$name") > "$R/arithmetic-$name.log" 2>&1 || return 1
  grep -E "PASS|mismatch|sigma table" "$R/arithmetic-$name.log" | tee "$R/arithmetic-$name.txt"
  grep -qx "packed arithmetic sigma table: $(tableOf "$name")" "$R/arithmetic-$name.log" || return 1
  echo "=== integration $name"
  python3 codegen/testpackedclient.py "./ecc2k130-$name" > "$R/integration-$name.log" 2>&1 || return 1
  grep -E "PASS|Error|assert" "$R/integration-$name.log" | tee "$R/integration-$name.txt"
}
build control || fail=1
build candidate || fail=1
if [ "$fail" != 0 ]; then
  tail -60 "$R"/build-*.log "$R"/arithmetic-*.log "$R"/integration-*.log 2>/dev/null
  exit "$fail"
fi
sha256sum ecc2k130-control ecc2k130-candidate > "$R/binary-sha256.txt"

verify() {
  local name=$1 status=0
  ./ecc2k130-$name --curve 131 --packed --threads "$VERIFY_THREADS" \
    --dp-weight 48 --dp-cap 262144 --steps 95 --launches 7 \
    --verify 300 --run-id 23 --dp-file "$R/dp-$name.bin" > "$R/verify-$name.log" 2>&1 || status=$?
  grep -E "MISMATCH|OVERFLOW|finished|packed sigma|resident|registers" "$R/verify-$name.log" | tee "$R/verify-$name.txt" || true
  grep -q "^packed sigma table: $(tableOf "$name") " "$R/verify-$name.log" || status=1
  grep -qx "packed witness: 0" "$R/verify-$name.log" || status=1
  grep -Eq "\(300 verified against the reference, 0 dropped\)" "$R/verify-$name.log" || status=1
  ! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-$name.log" || status=1
  return "$status"
}
verify control || fail=1
verify candidate || fail=1
python3 - "$R/dp-control.bin" "$R/dp-candidate.bin" <<'PY' | tee "$R/corpus-identity.txt" || fail=1
import hashlib, sys
def records(path):
    data = open(path, 'rb').read()
    n = len(data) // 32
    return sorted(data[32 * i:32 * i + 32] for i in range(n)), len(data)
a, la = records(sys.argv[1]); b, lb = records(sys.argv[2])
ha = hashlib.sha256(b''.join(a)).hexdigest(); hb = hashlib.sha256(b''.join(b)).hexdigest()
print(f'control {len(a)} records {la} bytes sha256 {ha}')
print(f'candidate {len(b)} records {lb} bytes sha256 {hb}')
print('PASS: sorted corpus identity' if a == b and len(a) > 0 else 'FAIL: corpora differ or empty')
sys.exit(0 if a == b and len(a) > 0 else 1)
PY
if [ "$fail" != 0 ]; then
  echo "correctness preflight failed; timing suppressed" | tee "$R/preflight.txt"
  exit "$fail"
fi
echo "PASS: arithmetic, integration, replay and sorted corpus identity" | tee "$R/preflight.txt"

printf 'phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 pair=$2 order=$3 variant=$4 log rc count rate digest state
  log="$R/${phase}-${pair}-${order}-${variant}.log"
  ./ecc2k130-$variant --curve 131 --packed --threads "$BENCH_THREADS" \
    --bench --steps 1024 --launches "$LAUNCHES" --verify 0 > "$log" 2>&1
  rc=$?
  count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  grep -q "^packed sigma table: $(tableOf "$variant") " "$log" || rc=1
  grep -qx "packed witness: 0" "$log" || rc=1
  grep -Eq "^backend cuda-packed131: $BENCH_THREADS threads x 16 slots x 1 lanes = $((BENCH_THREADS*16)) walks, dp weight 0, 1024 steps per launch$" "$log" || rc=1
  if [ "$count" != 1 ] || ! awk -v r="$rate" 'BEGIN{exit !(r+0>0)}'; then rc=1; fi
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$phase" "$pair" "$order" "$variant" "$rate" "$digest" "$state" | tee -a "$R/samples.tsv"
  return "$rc"
}
sample warmup 0 1 control || fail=1
sample warmup 0 2 candidate || fail=1
sample warmup 0 3 candidate || fail=1
sample warmup 0 4 control || fail=1
for pair in $(seq 1 "$PAIRS"); do
  if [ $((pair % 2)) = 1 ]; then
    sample ab "$pair" 1 control || fail=1
    sample ab "$pair" 2 candidate || fail=1
  else
    sample ab "$pair" 1 candidate || fail=1
    sample ab "$pair" 2 control || fail=1
  fi
done
if [ "$fail" != 0 ]; then
  echo "timing row failed" | tee -a "$R/failures.txt"
  exit "$fail"
fi
python3 - "$R/samples.tsv" "$R/result.json" <<'PY' | tee "$R/summary.txt"
import csv, json, statistics, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
ab = [r for r in rows if r['phase'] == 'ab']
pairs = {}
for r in ab:
    pairs.setdefault(int(r['pair']), {})[r['variant']] = float(r['rateMps'])
ratios = [pairs[p]['candidate'] / pairs[p]['control'] for p in sorted(pairs)]
c = [pairs[p]['control'] for p in sorted(pairs)]; k = [pairs[p]['candidate'] for p in sorted(pairs)]
out = dict(controlRatesMps=c, candidateRatesMps=k, controlMedianMps=statistics.median(c),
           candidateMedianMps=statistics.median(k), pairedRatios=ratios,
           pairedMedianRatio=statistics.median(ratios), pairedMinRatio=min(ratios), pairedMaxRatio=max(ratios),
           warmups=[(r['variant'], float(r['rateMps'])) for r in rows if r['phase'] == 'warmup'])
json.dump(out, open(sys.argv[2], 'w'), indent=1)
print('control median %.3f M/s, candidate median %.3f M/s' % (out['controlMedianMps'], out['candidateMedianMps']))
print('paired ratios', ['%.6f' % r for r in ratios], 'median %.6f min %.6f' % (out['pairedMedianRatio'], out['pairedMinRatio']))
PY
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo "=== done"
