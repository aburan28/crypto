#!/bin/bash
# B16/B32/B64 geometry sweep for the reconverged-hints table-v3 walk.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1
apt-get install -y -qq make g++ python3 >/dev/null 2>&1

R=${RESULTS:-/results}; mkdir -p "$R"
ARCH="PRO6000_ARCH=-gencode arch=compute_120,code=sm_120"
ARMS=(b16 b32 b64)
fail=0

batch() { case "$1" in b16) echo 16;; b32) echo 32;; b64) echo 64;; esac; }
block_threads() { case "$1" in b16) echo 512;; b32) echo 256;; b64) echo 128;; esac; }
runtime_threads() { case "$1" in b16) echo 96256;; b32) echo 48128;; b64) echo 24064;; esac; }

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
  echo "common live slots: 1540096; 32-launch slot populations: 49283072"
} | tee "$R/host.txt"
printf 'arm\tbatch\tblockThreads\truntimeThreads\tliveSlots\tscreenSlotPopulations\n' > "$R/geometry.tsv"
for arm in "${ARMS[@]}"; do
  b=$(batch "$arm"); bt=$(block_threads "$arm"); rt=$(runtime_threads "$arm")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$arm" "$b" "$bt" "$rt" "$((b*rt))" "$((b*rt*32))" >> "$R/geometry.tsv"
done
sha256sum Makefile include/packedkernels.cuh include/packedengine.cuh \
  include/packedtablewalk.cuh include/cycleanchor_body.h src/main.cu > "$R/source-files.sha256"

build() {
  local arm=$1 b bt
  b=$(batch "$arm"); bt=$(block_threads "$arm")
  echo "=== build $arm: BATCH=$b THREADS=$bt MINBLOCKS=1"
  make -s gpu-preset "$ARCH" KNOBS="TABLE_SPLIT_FORWARD=1 TABLE_BATCH_HINTS=1 BATCH=$b THREADS=$bt MINBLOCKS=1" \
    > "$R/build-$arm.log" 2>&1 || return 1
  [ -x ecc2k130 ] || return 1
  mv ecc2k130 "ecc2k130-$arm"
}
for arm in "${ARMS[@]}"; do build "$arm" || fail=1; done
if [ "$fail" != 0 ]; then tail -40 "$R"/build-*.log; exit "$fail"; fi
sha256sum ecc2k130-b16 ecc2k130-b32 ecc2k130-b64 > "$R/binary-sha256.txt"

verify() {
  local arm=$1 b rt status=0
  b=$(batch "$arm"); rt=$(runtime_threads "$arm")
  ./ecc2k130-$arm --curve 131 --packed --threads "$rt" --dp-weight 48 --dp-cap 262144 \
    --steps 96 --launches 6 --verify 300 --run-id 7 --dp-file "$R/dp-$arm.bin" \
    > "$R/verify-$arm.log" 2>&1 || status=$?
  grep -E "MISMATCH|OVERFLOW|finished|packed table (split forward|batch hints)|backend cuda-packed131" \
    "$R/verify-$arm.log" | tee "$R/verify-$arm.txt" || true
  grep -qx "packed table split forward: 1" "$R/verify-$arm.log" || status=1
  grep -qx "packed table batch hints: 1" "$R/verify-$arm.log" || status=1
  grep -Eq "^backend cuda-packed131: $rt threads x $b slots x 1 lanes = 1540096 walks," "$R/verify-$arm.log" || status=1
  grep -Eq "\(300 verified against the reference, 0 dropped\)" "$R/verify-$arm.log" || status=1
  ! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-$arm.log" || status=1
  return "$status"
}
for arm in "${ARMS[@]}"; do verify "$arm" || fail=1; done
python3 benchmarks/dp_identity.py "$R" b16 b32 b64 | tee "$R/dp-identity.txt" || fail=1
python3 benchmarks/hint-geometry/summarize.py "$R" --out "$R/preflight.json" --preflight || fail=1
if [ "$fail" != 0 ]; then echo "correctness preflight failed; timing suppressed"; exit "$fail"; fi

printf 'phase\tpair\torder\tvariant\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
sample() {
  local phase=$1 pair=$2 order=$3 arm=$4 launches=$5 b rt log rc rate count digest state
  b=$(batch "$arm"); rt=$(runtime_threads "$arm"); log="$R/${phase}-${pair}-${order}-${arm}.log"
  ./ecc2k130-$arm --curve 131 --packed --threads "$rt" --bench --steps 1024 --launches "$launches" --verify 0 > "$log" 2>&1
  rc=$?; count=$(grep -c '^[[:space:]]*finished: [0-9.][0-9.]* M it/s' "$log" || true)
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  grep -qx "packed table split forward: 1" "$log" || rc=1
  grep -qx "packed table batch hints: 1" "$log" || rc=1
  grep -Eq "^backend cuda-packed131: $rt threads x $b slots x 1 lanes = 1540096 walks," "$log" || rc=1
  if [ "$count" != 1 ] || ! awk -v r="$rate" 'BEGIN{exit !(r+0>0)}'; then rc=1; fi
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$phase" "$pair" "$order" "$arm" "$rate" "$digest" "$state" | tee -a "$R/samples.tsv"
  return "$rc"
}

# One excluded warmup per arm.
sample warmup 0 1 b16 16 || fail=1
sample warmup 0 2 b32 16 || fail=1
sample warmup 0 3 b64 16 || fail=1

# Equal-launch two-pass screen, reversing order on pass two.
sample screen 1 1 b16 32 || fail=1
sample screen 1 2 b32 32 || fail=1
sample screen 1 3 b64 32 || fail=1
sample screen 2 1 b64 32 || fail=1
sample screen 2 2 b32 32 || fail=1
sample screen 2 3 b16 32 || fail=1

python3 benchmarks/hint-geometry/summarize.py "$R" --out "$R/screen.json" || fail=1
top=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["timing"]["topCandidate"])' "$R/screen.json")
qualified=$(python3 -c 'import json,sys; print(int(json.load(open(sys.argv[1]))["timing"]["qualified"]))' "$R/screen.json")
if [ "$fail" = 0 ] && [ "$qualified" = 1 ]; then
  sample confirm 1 1 b16 64 || fail=1; sample confirm 1 2 "$top" 64 || fail=1
  sample confirm 2 1 "$top" 64 || fail=1; sample confirm 2 2 b16 64 || fail=1
  sample confirm 3 1 b16 64 || fail=1; sample confirm 3 2 "$top" 64 || fail=1
fi
python3 benchmarks/hint-geometry/summarize.py "$R" --out "$R/result.json" || fail=1
echo "=== done"
exit "$fail"
