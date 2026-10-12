#!/bin/bash
# Frozen five-pair both2 confirmation plus matched A/A control. See README.md.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1
apt-get install -y -qq make g++ >/dev/null 2>&1

R=${RESULTS:-/results}
mkdir -p "$R"
ARCH="PRO6000_ARCH=-gencode arch=compute_120,code=sm_120"
THREADS=${VERIFY_THREADS:-96256}
COMMON="TABLE_SPLIT_FORWARD=1 TABLE_BATCH_HINTS=1 TABLE_BLOCK_HINTS=1 TABLE_HINT_QUEUE=512 CYCLE_FAST2=1 BATCH=16 THREADS=512 MINBLOCKS=1"
fail=0

g++ -O2 -std=c++17 benchmarks/block-both2-confirm5/corpus_identity.cpp -o /tmp/corpus-identity || exit 1
g++ -O2 -std=c++17 benchmarks/block-both2-confirm5/summarize.cpp -o /tmp/both2-summarize || exit 1
{
  /tmp/corpus-identity --self-test
  /tmp/both2-summarize --self-test
} | tee "$R/helper-self-test.txt" || exit 1
cp /tmp/corpus-identity /tmp/both2-summarize "$R"/

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
  echo "verify threads: $THREADS; batch: 16; live slots: $((THREADS*16))"
  echo "protocol: five alternating A/B pairs plus five matched A/A pairs; 64 launches/sample"
} | tee "$R/host.txt"

sha256sum Makefile include/packedkernels.cuh include/packedengine.cuh \
  include/packedtablewalk.cuh include/tablewalk.h include/cycleanchor_body.h \
  include/ref.h src/main.cu benchmarks/block-both2-confirm5/README.md \
  benchmarks/block-both2-confirm5/gpujob.sh \
  benchmarks/block-both2-confirm5/corpus_identity.cpp \
  benchmarks/block-both2-confirm5/summarize.cpp > "$R/source-files.sha256"

build() {
  local name=$1 mode=$2 extra="PACKED_SQUARE_TABLE=0 PACKED_INV_POLY=0"
  [ "$mode" = 1 ] && extra="PACKED_SQUARE_TABLE=1 PACKED_INV_POLY=2"
  echo "=== build $name: $COMMON $extra"
  make -s gpu-preset "$ARCH" KNOBS="$COMMON $extra" \
    > "$R/build-$name.log" 2>&1 || return 1
  [ -x ecc2k130 ] || return 1
  mv ecc2k130 "ecc2k130-$name"
}

build control 0 || fail=1
build candidate 1 || fail=1
if [ "$fail" != 0 ]; then
  echo "build failure" | tee "$R/failures.txt"
  exit "$fail"
fi
cp ecc2k130-control ecc2k130-candidate "$R"/
sha256sum "$R/ecc2k130-control" "$R/ecc2k130-candidate" > "$R/binary-sha256.txt"

check_markers() {
  local log=$1 mode=$2 square=0 inv=0 shared=48732
  if [ "$mode" = 1 ]; then square=1; inv=2; shared=57052; fi
  grep -qx "packed table split forward: 1" "$log" || return 1
  grep -qx "packed table batch hints: 1" "$log" || return 1
  grep -qx "packed cycle fast2: 1" "$log" || return 1
  grep -qx "packed table block hints: 1, queue 512" "$log" || return 1
  grep -qx "packed table global: 0" "$log" || return 1
  grep -qx "packed table addend global: 0" "$log" || return 1
  grep -qx "packed square table: $square" "$log" || return 1
  grep -qx "packed polynomial inversion: $inv" "$log" || return 1
  grep -qx "packed table walk: 1 (8 branches, $shared shared bytes)" "$log" || return 1
  grep -qx "packed kernel: 128 registers/thread, 400 local bytes/thread, 1040 shared bytes/block, single-product multiplier" "$log" || return 1
  grep -Eq "^backend cuda-packed131: $THREADS threads x 16 slots x 1 lanes = $((THREADS*16)) walks," "$log" || return 1
  ! grep -Eq "MISMATCH|OVERFLOW" "$log"
}

verify() {
  local name=$1 mode=$2 status=0
  ./ecc2k130-$name --curve 131 --packed --threads "$THREADS" \
    --dp-weight 48 --dp-cap 262144 --steps 96 --launches 6 \
    --verify 300 --run-id 7 --dp-file "$R/dp-$name.bin" \
    > "$R/verify-$name.log" 2>&1 || status=$?
  check_markers "$R/verify-$name.log" "$mode" || status=1
  grep -Eq "[[:space:]]887095296 iterations[[:space:]]+1480278 dp[[:space:]]+1480278 stored[[:space:]]+0 dropped" "$R/verify-$name.log" || status=1
  grep -Eq "finished: [0-9.]+ M it/s, 1480278 distinguished points \(300 verified against the reference, 0 dropped\)" "$R/verify-$name.log" || status=1
  return "$status"
}

verify control 0 || fail=1
verify candidate 1 || fail=1
/tmp/corpus-identity "$R/dp-control.bin" "$R/dp-candidate.bin" \
  /tmp/sorted-control.bin /tmp/sorted-candidate.bin > "$R/dp-identity.txt" || fail=1
{
  sha256sum /tmp/sorted-control.bin /tmp/sorted-candidate.bin
  cmp /tmp/sorted-control.bin /tmp/sorted-candidate.bin && echo "sorted payloads byte-identical"
} >> "$R/dp-identity.txt" || fail=1
rm -f /tmp/sorted-control.bin /tmp/sorted-candidate.bin

if [ "$fail" != 0 ]; then
  echo "correctness/corpus preflight failed; timing suppressed" | tee -a "$R/failures.txt"
  exit "$fail"
fi
echo "PASS exact v3 replay and corpus preflight" > "$R/preflight.txt"

printf 'phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 pair=$2 order=$3 variant=$4 mode=0 binary=control launches=$5
  local log rc=0 rate count digest state final expected
  if [ "$variant" = candidate ]; then mode=1; binary=candidate; fi
  log="$R/${phase}-${pair}-${order}-${variant}.log"
  ./ecc2k130-$binary --curve 131 --packed --threads "$THREADS" \
    --bench --steps 1024 --launches "$launches" --verify 0 > "$log" 2>&1 || rc=$?
  check_markers "$log" "$mode" || rc=1
  count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  final=$(sed -nE 's/.*[[:space:]]([0-9]+) iterations.*/\1/p' "$log" | tail -1)
  expected=$((THREADS*16*1024*launches))
  [ "$count" = 1 ] || rc=1
  [ "$final" = "$expected" ] || rc=1
  awk -v r="$rate" 'BEGIN{exit !(r+0>0)}' || rc=1
  grep -Eq "finished: $rate M it/s, 0 distinguished points \(0 verified against the reference, 0 dropped\)" "$log" || rc=1
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu \
    --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$phase" "$pair" "$order" "$variant" "$rate" "$digest" "$state" \
    | tee -a "$R/samples.tsv"
  return "$rc"
}

# Two excluded warmups.
sample warmup 0 1 control 16 || fail=1
sample warmup 0 2 candidate 16 || fail=1

# Fixed interleaved five-round schedule from README.md.
sample ab 1 1 control 64 || fail=1
sample ab 1 2 candidate 64 || fail=1
sample aa 1 1 aa1 64 || fail=1
sample aa 1 2 aa2 64 || fail=1

sample aa 2 1 aa2 64 || fail=1
sample aa 2 2 aa1 64 || fail=1
sample ab 2 1 candidate 64 || fail=1
sample ab 2 2 control 64 || fail=1

sample ab 3 1 control 64 || fail=1
sample ab 3 2 candidate 64 || fail=1
sample aa 3 1 aa1 64 || fail=1
sample aa 3 2 aa2 64 || fail=1

sample aa 4 1 aa2 64 || fail=1
sample aa 4 2 aa1 64 || fail=1
sample ab 4 1 candidate 64 || fail=1
sample ab 4 2 control 64 || fail=1

sample ab 5 1 control 64 || fail=1
sample ab 5 2 candidate 64 || fail=1
sample aa 5 1 aa1 64 || fail=1
sample aa 5 2 aa2 64 || fail=1

sha256sum "$R"/*.log > "$R/log-sha256.txt"
if [ "$fail" != 0 ]; then
  echo "one or more timed samples failed exact-work/marker validation" | tee -a "$R/failures.txt"
  exit "$fail"
fi
/tmp/both2-summarize "$R/samples.tsv" "$R/result.json" || exit 1
sha256sum "$R/result.json" "$R/samples.tsv" "$R/dp-identity.txt" > "$R/result-sha256.txt"
echo "=== done"
