#!/bin/bash
# Build and screen the opt-in r-adding table-walk redesign against the
# shipping sigma walk. This is a throughput/correctness screen only: the
# redesign has a different v3 iteration function and corpus identity.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
R=${RESULTS:-/results}
mkdir -p "$R"
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCH="-gencode arch=compute_${CAP},code=sm_${CAP}"
VERIFY_REPORTS=${VERIFY_REPORTS:-300}
VERIFY_STEPS=${VERIFY_STEPS:-16}
VERIFY_LAUNCHES=${VERIFY_LAUNCHES:-6}
VERIFY_DP_CAP=${VERIFY_DP_CAP:-262144}
LAUNCHES=${LAUNCHES:-32}
REPEATS=${REPEATS:-3}
REDESIGN_PIVOT=${REDESIGN_PIVOT:-1}
fail=0

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
  echo "verify reports: $VERIFY_REPORTS; verify dp cap: $VERIFY_DP_CAP; verify steps: $VERIFY_STEPS; verify launches: $VERIFY_LAUNCHES; benchmark launches: $LAUNCHES; repeats: $REPEATS; redesign pivot bytes: $REDESIGN_PIVOT"
} | tee "$R/host.txt"

sha256sum Makefile include/packedkernels.cuh include/packedtablewalk.cuh include/tablewalk.h \
  include/packedengine.cuh include/packed131.h src/main.cu src/testpackedcuda.cu \
  src/testtablewalkcuda.cu benchmarks/iteration-redesign/gpujob.sh \
  > "$R/source-files.sha256"

COMMON=(
  BATCH=16 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1
  PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1
  PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 PACKED_DIRECT_REDUCE=1
  PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=1 PACKED_STATE_TILE=256
  PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1
)
knobs() {
  case "$1" in
    shipping) echo "THREADS=256 MINBLOCKS=2 WITNESS=1 WALK_TABLE=0" ;;
    redesign) echo "THREADS=256 MINBLOCKS=2 WITNESS=0 WALK_TABLE=1 TABLE_PIVOT_BYTES=$REDESIGN_PIVOT" ;;
  esac
}

build() {
  local name=$1
  # shellcheck disable=SC2046
  make -s -B ecc2k130 ARCH="$ARCH" "${COMMON[@]}" $(knobs "$name") \
    > "$R/build-$name.log" 2>&1 || return 1
  grep -E "Function properties for.*walk|registers|spill|stack frame" \
    "$R/build-$name.log" | tail -16 | tee "$R/build-$name.txt" || true
  [ -x ecc2k130 ] || return 1
  mv ecc2k130 "ecc2k130-$name"
}

for arm in shipping redesign; do
  build "$arm" || fail=1
done
if [ "$fail" != 0 ]; then
  tail -100 "$R"/build-*.log 2>/dev/null
  exit "$fail"
fi
sha256sum ecc2k130-shipping ecc2k130-redesign > "$R/binary-sha256.txt"

# The device table-walk suite exercises the redesigned selector and addend
# tables against the independent reference implementation.
# shellcheck disable=SC2046
make -s test-table-walk-cuda ARCH="$ARCH" "${COMMON[@]}" $(knobs redesign) \
  > "$R/arithmetic-redesign.log" 2>&1 || fail=1
grep -E "PASS|mismatch|table walk" "$R/arithmetic-redesign.log" \
  | tee "$R/arithmetic-redesign.txt" || true

verify() {
  local name=$1 status=0
  ./ecc2k130-$name --curve 131 --packed --dp-weight 50 \
    --dp-cap "$VERIFY_DP_CAP" --steps "$VERIFY_STEPS" --launches "$VERIFY_LAUNCHES" --verify "$VERIFY_REPORTS" \
    --run-id 35 --dp-file "$R/dp-$name.bin" > "$R/verify-$name.log" 2>&1 || status=$?
  grep -E "MISMATCH|OVERFLOW|finished|table walk|resident|registers" \
    "$R/verify-$name.log" | tee "$R/verify-$name.txt" || true
  if [ "$name" = redesign ]; then
    grep -q '^packed table walk: 1 ' "$R/verify-$name.log" || status=1
    grep -Eq "\($VERIFY_REPORTS verified against the reference, 0 dropped\)" \
      "$R/verify-$name.log" || status=1
  else
    grep -q '^packed table walk: 0 ' "$R/verify-$name.log" || status=1
  fi
  ! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-$name.log" || status=1
  return "$status"
}
for arm in shipping redesign; do
  verify "$arm" || fail=1
done

if [ "$fail" != 0 ]; then
  echo "correctness preflight failed; timing suppressed" | tee "$R/preflight.txt"
  exit "$fail"
fi
echo "PASS: arithmetic, table-walk replay and v3 verification" | tee "$R/preflight.txt"

printf 'arm\trepeat\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local arm=$1 repeat=$2 log rc rate digest state
  log="$R/bench-${arm}-${repeat}.log"
  ./ecc2k130-$arm --curve 131 --packed --bench --steps 1024 --launches "$LAUNCHES" \
    > "$log" 2>&1
  rc=$?
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  if [ "$(printf '%s\n' "$rate" | wc -l)" != 1 ] || ! awk -v r="$rate" 'BEGIN{exit !(r+0>0)}'; then
    rc=1
  fi
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\n' "$arm" "$repeat" "$rate" "$digest" "$state" | tee -a "$R/samples.tsv"
  return "$rc"
}

for repeat in $(seq 1 "$REPEATS"); do
  sample shipping "$repeat" || fail=1
  sample redesign "$repeat" || fail=1
done
if [ "$fail" != 0 ]; then
  echo "timing row failed" | tee "$R/failures.txt"
  exit "$fail"
fi

awk -F '\t' -v expected="$REPEATS" '
  NR > 1 { n[$1]++; vals[$1,n[$1]] = $3 }
  END {
    for (a in n) {
      if (n[a] == expected) {
        if (vals[a,1] > vals[a,2]) { t=vals[a,1]; vals[a,1]=vals[a,2]; vals[a,2]=t }
        if (vals[a,2] > vals[a,3]) { t=vals[a,2]; vals[a,2]=vals[a,3]; vals[a,3]=t }
        if (vals[a,1] > vals[a,2]) { t=vals[a,1]; vals[a,1]=vals[a,2]; vals[a,2]=t }
        printf "%s median %.3f M/s range %.3f-%.3f M/s\n", a, vals[a,2], vals[a,1], vals[a,3]
      } else printf "%s incomplete (%d samples)\n", a, n[a]
    }
  }
' "$R/samples.tsv" | tee "$R/summary.txt"
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo "=== done"
