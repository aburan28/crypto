#!/bin/bash
# Run the full tuned v3 redesign recipe, rather than the deliberately small
# arm used by gpujob.sh for quick selector screening. This is a benchmark
# harness only: it never writes campaign output or DP-ingester state.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
R=${RESULTS:-/results}
mkdir -p "$R"
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCH="-gencode arch=compute_${CAP},code=sm_${CAP}"
VERIFY_REPORTS=${VERIFY_REPORTS:-1}
VERIFY_STEPS=${VERIFY_STEPS:-1}
VERIFY_LAUNCHES=${VERIFY_LAUNCHES:-1}
VERIFY_DP_CAP=${VERIFY_DP_CAP:-16384}
BENCH_STEPS=${BENCH_STEPS:-8}
BENCH_LAUNCHES=${BENCH_LAUNCHES:-1}
REPEATS=${REPEATS:-3}
fail=0

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-unknown}"
  echo "verify reports: $VERIFY_REPORTS; verify dp cap: $VERIFY_DP_CAP; verify steps: $VERIFY_STEPS; verify launches: $VERIFY_LAUNCHES; benchmark steps: $BENCH_STEPS; benchmark launches: $BENCH_LAUNCHES; repeats: $REPEATS"
} | tee "$R/host.txt"

sha256sum Makefile include/packedkernels.cuh include/packedtablewalk.cuh \
  include/tablewalk.h include/packedengine.cuh include/packed131.h src/main.cu \
  src/testpackedcuda.cu src/testtablewalkcuda.cu \
  benchmarks/iteration-redesign/tunedjob.sh > "$R/source-files.sha256"

make -s -B gpu-rtx-pro6000-iteration-redesign PRO6000_ARCH="$ARCH" \
  > "$R/build-tuned.log" 2>&1 || fail=1
grep -E "Function properties for.*walk|registers|spill|stack frame|packed table walk" \
  "$R/build-tuned.log" | tail -24 | tee "$R/build-tuned.txt" || true
[ -x ecc2k130 ] && mv ecc2k130 ecc2k130-tuned || fail=1

if [ "$fail" != 0 ]; then
  tail -120 "$R/build-tuned.log"
  exit "$fail"
fi

./ecc2k130-tuned --curve 131 --packed --dp-weight 50 \
  --dp-cap "$VERIFY_DP_CAP" --steps "$VERIFY_STEPS" \
  --launches "$VERIFY_LAUNCHES" --verify "$VERIFY_REPORTS" \
  --run-id 35 --dp-file "$R/dp-tuned.bin" > "$R/verify-tuned.log" 2>&1 || fail=1
grep -E "MISMATCH|OVERFLOW|finished|table walk|resident|registers" \
  "$R/verify-tuned.log" | tee "$R/verify-tuned.txt" || true
grep -q '^packed table walk: 1 ' "$R/verify-tuned.log" || fail=1
grep -Eq "\($VERIFY_REPORTS verified against the reference, 0 dropped\)" \
  "$R/verify-tuned.log" || fail=1
! grep -Eq "MISMATCH|OVERFLOW" "$R/verify-tuned.log" || fail=1

if [ "$fail" != 0 ]; then
  echo "correctness preflight failed; timing suppressed" | tee "$R/preflight.txt"
  exit "$fail"
fi
echo "PASS: tuned v3 table-walk verification" | tee "$R/preflight.txt"

printf 'repeat\trateMps\tlogSha256\tgpuState\n' > "$R/samples.tsv"
for repeat in $(seq 1 "$REPEATS"); do
  log="$R/bench-$repeat.log"
  ./ecc2k130-tuned --curve 131 --packed --bench \
    --steps "$BENCH_STEPS" --launches "$BENCH_LAUNCHES" > "$log" 2>&1 || fail=1
  rate=$(sed -nE 's/^[[:space:]]*finished: ([0-9.]+) M it\/s.*/\1/p' "$log")
  if [ "$(printf '%s\n' "$rate" | wc -l)" != 1 ] || ! awk -v r="$rate" 'BEGIN{exit !(r+0>0)}'; then
    fail=1
  fi
  digest=$(sha256sum "$log" | awk '{print $1}')
  state=$(nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader,nounits | head -1 | tr '\t' ' ')
  printf '%s\t%s\t%s\t%s\n' "$repeat" "$rate" "$digest" "$state" | tee -a "$R/samples.tsv"
done

if [ "$fail" != 0 ]; then
  echo "timing row failed" | tee "$R/failures.txt"
  exit "$fail"
fi

awk -F '\t' -v expected="$REPEATS" '
  NR > 1 { n++; vals[n] = $2 }
  END {
    if (n == expected) {
      for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++)
        if (vals[j] < vals[i]) { t=vals[i]; vals[i]=vals[j]; vals[j]=t }
      printf "tuned median %.3f M/s range %.3f-%.3f M/s\n", vals[int((n+1)/2)], vals[1], vals[n]
    } else printf "tuned incomplete (%d samples)\n", n
  }
' "$R/samples.tsv" | tee "$R/summary.txt"
sha256sum "$R"/* > "$R/artifact-files.sha256"
echo "=== done"
