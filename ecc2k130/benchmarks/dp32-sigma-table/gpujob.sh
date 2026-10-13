#!/bin/bash
# Screen the exact polynomial sigma table in the ordinary witness-carrying
# DP32 walk. The three arms separate live geometry from the candidate's
# one-block geometry and the table change. Run under runpod_job.py.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
R=${RESULTS:-/results}
mkdir -p "$R"
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCH="-gencode arch=compute_${CAP},code=sm_${CAP}"
VERIFY_THREADS=96256
VERIFY_REPORTS=5
VERIFY_LAUNCHES=2
BENCH_THREADS=${BENCH_THREADS:-385024}
LAUNCHES=${LAUNCHES:-16}
PAIRS=${PAIRS:-3}
fail=0

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
  echo "verify threads: $VERIFY_THREADS; verify reports: $VERIFY_REPORTS; verify launches: $VERIFY_LAUNCHES; benchmark threads: $BENCH_THREADS; launches: $LAUNCHES; pairs: $PAIRS"
} | tee "$R/host.txt"

sha256sum Makefile include/packedkernels.cuh include/packedengine.cuh include/packed131.h \
  include/packedsigma131.h src/main.cu src/testpackedcuda.cu \
  benchmarks/dp32-sigma-table/gpujob.sh benchmarks/sigma-fused/corpus_identity.cpp \
  > "$R/source-files.sha256"

g++ -O2 -std=c++17 -Wall -Wextra -Werror \
  benchmarks/sigma-fused/corpus_identity.cpp -o /tmp/corpus-identity || exit 1

COMMON=(
  BATCH=16 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
  PACKED_INLINE_POLY=3 SIGMA_FUSED=0 WITNESS=1 WALK_TABLE=0
)
knobs() {
  case "$1" in
    live)      echo "THREADS=256 MINBLOCKS=2 PACKED_SIGMA_TABLE=0" ;;
    geometry)  echo "THREADS=512 MINBLOCKS=1 PACKED_SIGMA_TABLE=0" ;;
    candidate) echo "THREADS=512 MINBLOCKS=1 PACKED_SIGMA_TABLE=1" ;;
  esac
}
tableOf() { case "$1" in live|geometry) echo 0 ;; candidate) echo 1 ;; esac; }

build() {
  local name=$1
  echo "=== build $name $(knobs "$name")"
  # shellcheck disable=SC2046
  make -s -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" $(knobs "$name") \
    > "$R/build-$name.log" 2>&1 || return 1
  grep -E "Function properties for.*walk|registers|spill|stack frame" \
    "$R/build-$name.log" | tail -12 | tee "$R/build-$name.txt" || true
  [ -x ecc2k130 ] || return 1
  mv ecc2k130 "ecc2k130-$name"
  echo "=== arithmetic $name"
  # shellcheck disable=SC2046
  make -s test-packed-cuda ARCH="$ARCH" "${COMMON[@]}" $(knobs "$name") \
    > "$R/arithmetic-$name.log" 2>&1 || return 1
  grep -E "PASS|mismatch|sigma table" "$R/arithmetic-$name.log" | tee "$R/arithmetic-$name.txt"
  grep -qx "packed arithmetic sigma table: $(tableOf "$name")" "$R/arithmetic-$name.log" || return 1
  # testpackedclient.py predates the witness v2 framing used by DP32. The
  # device arithmetic suite above and the v2 replay/corpus gate below cover
  # this experiment without silently parsing v2 output as v1.
}

for arm in live geometry candidate; do
  build "$arm" || fail=1
done
if [ "$fail" != 0 ]; then
  tail -80 "$R"/build-*.log "$R"/arithmetic-*.log 2>/dev/null
  exit "$fail"
fi
sha256sum ecc2k130-live ecc2k130-geometry ecc2k130-candidate > "$R/binary-sha256.txt"

verify() {
  local name=$1 status=0
  ./ecc2k130-$name --curve 131 --packed --threads "$VERIFY_THREADS" \
    --dp-weight 32 --dp-cap 262144 --steps 1024 --launches "$VERIFY_LAUNCHES" \
    --verify "$VERIFY_REPORTS" --run-id 32 --dp-file "$R/dp-$name.bin" \
    > "$R/verify-$name.log" 2>&1 || status=$?
  grep -E "MISMATCH|OVERFLOW|finished|packed sigma|resident|registers" \
    "$R/verify-$name.log" | tee "$R/verify-$name.txt" || true
  grep -q "^packed sigma table: $(tableOf "$name") " "$R/verify-$name.log" || status=1
  grep -qx "packed witness: 1" "$R/verify-$name.log" || status=1
  grep -Eq "\($VERIFY_REPORTS verified against the reference, 0 dropped\)" "$R/verify-$name.log" || status=1
  ! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-$name.log" || status=1
  return "$status"
}
for arm in live geometry candidate; do verify "$arm" || fail=1; done
/tmp/corpus-identity "$R/dp-live.bin" "$R/dp-geometry.bin" "$R/dp-candidate.bin" \
  | tee "$R/corpus-identity.txt" || fail=1
if [ "$fail" != 0 ]; then
  echo "correctness preflight failed; timing suppressed" | tee "$R/preflight.txt"
  exit "$fail"
fi
echo "PASS: arithmetic, integration, replay and sorted v2 corpus identity" | tee "$R/preflight.txt"

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
  grep -qx "packed witness: 1" "$log" || rc=1
  grep -Eq "^backend cuda-packed131: $BENCH_THREADS threads x 16 slots x 1 lanes = $((BENCH_THREADS*16)) walks, dp weight 0, 1024 steps per launch$" "$log" || rc=1
  if [ "$count" != 1 ] || ! awk -v r="$rate" 'BEGIN{exit !(r+0>0)}'; then rc=1; fi
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$phase" "$pair" "$order" "$variant" "$rate" "$digest" "$state" | tee -a "$R/samples.tsv"
  return "$rc"
}

# Equal-work warmups, then three rotating three-arm panels.
sample warmup 0 1 live || fail=1
sample warmup 0 2 geometry || fail=1
sample warmup 0 3 candidate || fail=1
for pair in $(seq 1 "$PAIRS"); do
  case $((pair % 3)) in
    1) order="live geometry candidate" ;;
    2) order="geometry candidate live" ;;
    0) order="candidate live geometry" ;;
  esac
  pos=1
  for arm in $order; do
    sample ab "$pair" "$pos" "$arm" || fail=1
    pos=$((pos + 1))
  done
done
if [ "$fail" != 0 ]; then
  echo "timing row failed" | tee -a "$R/failures.txt"
  exit "$fail"
fi

python3 - "$R/samples.tsv" "$R/result.json" <<'PY' | tee "$R/summary.txt"
import csv, json, statistics, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
ab = [r for r in rows if r['phase'] == 'ab']
by = {}
for r in ab:
    by.setdefault(r['variant'], []).append(float(r['rateMps']))
med = {k: statistics.median(v) for k, v in by.items()}
live = med['live']; geometry = med['geometry']; candidate = med['candidate']
out = {
    'ratesMps': by, 'mediansMps': med,
    'candidateOverLive': candidate / live,
    'candidateOverGeometry': candidate / geometry,
    'pairs': len(set(r['pair'] for r in ab)),
    'launches': 16,
}
json.dump(out, open(sys.argv[2], 'w'), indent=1)
for k in ('live', 'geometry', 'candidate'):
    print('%s median %.3f M/s' % (k, med[k]))
print('candidate/live %.6f candidate/geometry %.6f' % (out['candidateOverLive'], out['candidateOverGeometry']))
PY
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo "=== done"
