#!/bin/bash
# Bounded per-lane vs block-compacted hint screen. See README.md.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1
apt-get install -y -qq make g++ python3 >/dev/null 2>&1

R=${RESULTS:-/results}
mkdir -p "$R"
ARCH="PRO6000_ARCH=-gencode arch=compute_120,code=sm_120"
THREADS=${VERIFY_THREADS:-96256}
COMMON="TABLE_SPLIT_FORWARD=1 TABLE_BATCH_HINTS=1 CYCLE_FAST2=1 BATCH=16 THREADS=512 MINBLOCKS=1"
fail=0

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
  echo "verify threads: $THREADS; batch: 16; live slots: $((THREADS*16))"
} | tee "$R/host.txt"

sha256sum Makefile include/packedkernels.cuh include/packedengine.cuh \
  include/packedtablewalk.cuh include/tablewalk.h include/cycleanchor_body.h \
  include/ref.h src/main.cu > "$R/source-files.sha256"

build() {
  local name=$1 block=$2
  echo "=== build $name: $COMMON TABLE_BLOCK_HINTS=$block TABLE_HINT_QUEUE=512"
  make -s gpu-preset "$ARCH" KNOBS="$COMMON TABLE_BLOCK_HINTS=$block TABLE_HINT_QUEUE=512" \
    > "$R/build-$name.log" 2>&1 || return 1
  if [ ! -x ecc2k130 ]; then
    echo "BUILD FAILED $name" | tee -a "$R/failures.txt"
    return 1
  fi
  mv ecc2k130 "ecc2k130-$name"
}

build control 0 || fail=1
build candidate 1 || fail=1
if [ "$fail" != 0 ]; then
  tail -40 "$R"/build-*.log
  exit "$fail"
fi
sha256sum ecc2k130-control ecc2k130-candidate > "$R/binary-sha256.txt"

verify() {
  local name=$1 block=$2 status=0
  ./ecc2k130-$name --curve 131 --packed --threads "$THREADS" \
    --dp-weight 48 --dp-cap 262144 --steps 96 --launches 6 \
    --verify 300 --run-id 7 --dp-file "$R/dp-$name.bin" \
    > "$R/verify-$name.log" 2>&1 || status=$?
  grep -E "MISMATCH|OVERFLOW|finished|packed table (split forward|batch hints|block hints)|packed cycle fast2|resident" \
    "$R/verify-$name.log" | tee "$R/verify-$name.txt" || true
  grep -qx "packed table split forward: 1" "$R/verify-$name.log" || status=1
  grep -qx "packed table batch hints: 1" "$R/verify-$name.log" || status=1
  grep -qx "packed cycle fast2: 1" "$R/verify-$name.log" || status=1
  grep -qx "packed table block hints: $block, queue 512" "$R/verify-$name.log" || status=1
  grep -Eq "^backend cuda-packed131: $THREADS threads x 16 slots x 1 lanes = $((THREADS*16)) walks," \
    "$R/verify-$name.log" || status=1
  grep -Eq "\(300 verified against the reference, 0 dropped\)" "$R/verify-$name.log" || status=1
  ! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-$name.log" || status=1
  return "$status"
}

verify control 0 || fail=1
verify candidate 1 || fail=1
python3 benchmarks/dp_identity.py "$R" control candidate \
  | tee "$R/dp-identity.txt" || fail=1

python3 benchmarks/block-hints/summarize.py "$R" \
  --out "$R/preflight.json" --preflight || fail=1
if [ "$fail" != 0 ]; then
  echo "correctness/corpus preflight failed; timing suppressed" | tee -a "$R/failures.txt"
  exit "$fail"
fi

printf 'phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 pair=$2 order=$3 name=$4 launches=$5 block log rc rate count digest state
  block=0; [ "$name" = candidate ] && block=1
  log="$R/${phase}-${pair}-${order}-${name}.log"
  ./ecc2k130-$name --curve 131 --packed --threads "$THREADS" \
    --bench --steps 1024 --launches "$launches" --verify 0 > "$log" 2>&1
  rc=$?
  count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  grep -qx "packed table split forward: 1" "$log" || rc=1
  grep -qx "packed table batch hints: 1" "$log" || rc=1
  grep -qx "packed cycle fast2: 1" "$log" || rc=1
  grep -qx "packed table block hints: $block, queue 512" "$log" || rc=1
  grep -Eq "^backend cuda-packed131: $THREADS threads x 16 slots x 1 lanes = $((THREADS*16)) walks," \
    "$log" || rc=1
  if [ "$count" != 1 ] || ! awk -v r="$rate" 'BEGIN{exit !(r+0>0)}'; then rc=1; fi
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu \
    --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$phase" "$pair" "$order" "$name" "$rate" "$digest" "$state" \
    | tee -a "$R/samples.tsv"
  return "$rc"
}

# Exactly two excluded warmups.
sample warmup 0 1 control 16 || fail=1
sample warmup 0 2 candidate 16 || fail=1

# Candidate must clear the faster bracketing control by at least 0.5%.
sample screen 0 1 control 32 || fail=1
sample screen 0 2 candidate 32 || fail=1
sample screen 0 3 control 32 || fail=1

python3 benchmarks/block-hints/summarize.py "$R" --out "$R/screen.json" || fail=1
qualified=$(python3 -c 'import json,sys; print(int(json.load(open(sys.argv[1]))["timing"]["qualified"]))' "$R/screen.json")

if [ "$fail" = 0 ] && [ "$qualified" = 1 ]; then
  sample confirm 1 1 control 64 || fail=1
  sample confirm 1 2 candidate 64 || fail=1
  sample confirm 2 1 candidate 64 || fail=1
  sample confirm 2 2 control 64 || fail=1
  sample confirm 3 1 control 64 || fail=1
  sample confirm 3 2 candidate 64 || fail=1
fi

python3 benchmarks/block-hints/summarize.py "$R" --out "$R/result.json" || fail=1
echo "=== done"
exit "$fail"
