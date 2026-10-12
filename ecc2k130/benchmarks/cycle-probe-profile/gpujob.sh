#!/bin/bash
# Diagnostic-only exact table-v3 cold-probe counts on the B16 fast2 build.
set -euo pipefail
cd "${WORK:-/work}"

R=${RESULTS:-/results}
mkdir -p "$R"
ARCH="PRO6000_ARCH=-gencode arch=compute_120,code=sm_120"
THREADS=96256
BATCH=16
STEPS=64
DP_WEIGHT=0
DP_CAP=262144
COMMON="BATCH=$BATCH THREADS=512 MINBLOCKS=1 TABLE_SPLIT_FORWARD=1 TABLE_BATCH_HINTS=1 CYCLE_FAST2=1 CYCLE_PROFILE=1"

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-$(cat /work/SOURCE_REV 2>/dev/null || echo unknown)}"
  echo "diagnostic only: instrumented timing is not performance evidence"
  echo "profile geometry: $THREADS threads x $BATCH slots; $STEPS steps/launch; dp weight $DP_WEIGHT"
  echo "launch budgets: early=1 medium=8 long=32"
} | tee "$R/host.txt"

sha256sum Makefile include/kernel.h include/packedengine.cuh include/packedkernels.cuh \
  include/packedtablewalk.cuh include/tablewalk.h include/cycleanchor_body.h \
  src/main.cu src/testcycleprofile.cpp \
  benchmarks/cycle-probe-profile/gpujob.sh \
  benchmarks/cycle-probe-profile/README.md > "$R/source-files.sha256"

make test-cycle-profile > "$R/host-controls.log" 2>&1
grep -E '^PASS: hints ' "$R/host-controls.log" | tee "$R/host-controls.txt"
test "$(grep -c '^PASS: hints ' "$R/host-controls.log")" -eq 2

echo "=== build profile: $COMMON"
make -s gpu-preset "$ARCH" KNOBS="$COMMON" > "$R/build.log" 2>&1
test -x ecc2k130
mv ecc2k130 ecc2k130-cycle-profile
sha256sum ecc2k130-cycle-profile > "$R/binary-sha256.txt"

printf 'phase\tlaunches\tstepsPerLaunch\texactUpdates\tlogSha256\tcountersJson\n' > "$R/counts.tsv"
run_profile() {
  local phase=$1 launches=$2 log json digest updates
  log="$R/$phase.log"
  ./ecc2k130-cycle-profile --curve 131 --packed --threads "$THREADS" \
    --bench --steps "$STEPS" --launches "$launches" --verify 0 \
    --dp-weight "$DP_WEIGHT" --dp-cap "$DP_CAP" --run-id 7 > "$log" 2>&1
  grep -qx "packed table split forward: 1" "$log"
  grep -qx "packed table batch hints: 1" "$log"
  grep -qx "packed cycle fast2: 1" "$log"
  grep -qx "packed cycle profile: 1" "$log"
  test "$(grep -c '^cycle profile: ' "$log")" -eq 1
  json=$(sed -n 's/^cycle profile: //p' "$log")
  printf '%s\n' "$json" | grep -q '"reconciled":true'
  digest=$(sha256sum "$log" | awk '{print $1}')
  updates=$((THREADS * BATCH * STEPS * launches))
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$phase" "$launches" "$STEPS" "$updates" "$digest" "$json" \
    | tee -a "$R/counts.tsv"
}

# Each invocation starts from the same run id and initial state. Medium and
# long therefore contain the exact early prefix; their only changed input is
# the fixed launch budget.
run_profile early 1
run_profile medium 8
run_profile long 32

echo "=== diagnostic profile complete; no timing claim is admissible"
