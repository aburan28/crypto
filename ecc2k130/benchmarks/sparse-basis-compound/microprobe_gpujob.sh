#!/bin/bash
set -euo pipefail

cd /work
R=${RESULTS:-/results}
mkdir -p "$R"
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1
apt-get install -y -qq build-essential libssl-dev >/dev/null 2>&1

GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
GPU_COUNT=$(nvidia-smi --query-gpu=count --format=csv,noheader 2>/dev/null | head -1 || true)
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
if [ "$GPU_NAME" != "NVIDIA RTX PRO 6000 Blackwell Server Edition" ] || [ "$CAP" != 120 ]; then
  echo "frozen GPU mismatch: $GPU_NAME sm_$CAP" >&2
  exit 2
fi
if [ -n "$GPU_COUNT" ] && [ "$GPU_COUNT" != 1 ]; then
  echo "frozen GPU count mismatch: $GPU_COUNT" >&2
  exit 2
fi
nvcc --version | grep -q 'release 13\.3, V13\.3\.73' || {
  nvcc --version >&2
  echo "frozen CUDA compiler mismatch" >&2
  exit 2
}
case "${SOURCE_REV:-}" in
  [0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f][0-9a-f]) ;;
  *) echo "SOURCE_REV must be an exact clean 40-hex revision" >&2; exit 2 ;;
esac

{
  nvidia-smi --query-gpu=name,uuid,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: $SOURCE_REV"
  echo "probe: 188 blocks x 512 threads x 65536 rounds"
  echo "dynamic table bytes: 77440; driver reservation expected: 1024"
} | tee "$R/host.txt"

DIR=benchmarks/sparse-basis-compound
SOURCE_SHA=$(sha256sum "$DIR/microprobe.cu" | awk '{print $1}')
sha256sum \
  "$DIR/MICROPROBE-PROTOCOL.md" "$DIR/PROTOCOL.md" "$DIR/RESULTS.md" \
  "$DIR/audit.cpp" "$DIR/microprobe.cu" "$DIR/summarize_microprobe.cpp" \
  "$DIR/audit_microprobe_sass.cpp" "$DIR/microprobe_gpujob.sh" \
  include/bitslice.h include/packed131.h include/packedtransform131.h include/packedsigma131.h \
  include/packeddirectreduce131.h > "$R/source-files.sha256"

g++ -O2 -std=c++17 -Wall -Wextra -Werror "$DIR/summarize_microprobe.cpp" \
  -o /tmp/summarize-microprobe
g++ -O2 -std=c++17 -Wall -Wextra -Werror "$DIR/audit_microprobe_sass.cpp" \
  -o /tmp/audit-microprobe-sass

nvcc -O3 -std=c++17 -gencode arch=compute_120,code=sm_120 -lineinfo \
  -Xptxas -v -DECC_PACKED_PERM_SIGMA=1 -DECC_PACKED_FROM_REDUCED=1 \
  -DECC_PACKED_DIRECT_REDUCE=1 -DECC_PACKED_INLINE_POLY=3 \
  "-DMICROPROBE_SOURCE_SHA=\"$SOURCE_SHA\"" \
  "$DIR/microprobe.cu" -o "$R/microprobe" -lcrypto \
  > "$R/build.log" 2>&1
test -x "$R/microprobe"
grep -q "$SOURCE_SHA" < <(strings "$R/microprobe")
sha256sum "$R/microprobe" > "$R/binary.sha256"
cuobjdump --dump-resource-usage "$R/microprobe" > "$R/cuobjdump-resources.txt" 2>&1

mkdir -p "$R/code"
(cd "$R/code" && cuobjdump --extract-elf all "$R/microprobe" > extract.log 2>&1)
CUBIN=$(find "$R/code" -maxdepth 1 -type f -name '*.cubin' -print)
if [ "$(printf '%s\n' "$CUBIN" | sed '/^$/d' | wc -l)" -ne 1 ]; then
  echo "expected exactly one native cubin" >&2
  exit 1
fi
nvdisasm -c "$CUBIN" > "$R/microprobe.sass"
/tmp/audit-microprobe-sass "$R/microprobe.sass" "$R/sass-audit.tsv" \
  | tee "$R/sass-audit.txt"
sha256sum "$CUBIN" "$R/microprobe.sass" "$R/sass-audit.tsv" \
  > "$R/native-code.sha256"

"$R/microprobe" "$R" | tee "$R/probe.log"
grep -qx "source_sha256=$SOURCE_SHA" "$R/probe.log"
grep -qx "PASS: exact-layout sparse table multicast probe" "$R/probe.log"
grep -q '^PASS:' "$R/correctness.txt"
test "$(tail -n +2 "$R/key-histogram.tsv" | wc -l)" -eq 64
test "$(tail -n +2 "$R/samples.tsv" | wc -l)" -eq 42
test "$(tail -n +2 "$R/resources.tsv" | wc -l)" -eq 7

cat > "$R/gates.txt" <<EOF
correctness=PASS
key_coverage=PASS
resources=PASS
sass=PASS
source_binding=PASS
binary_binding=PASS
EOF

cat > "$R/bindings.txt" <<EOF
sourceRevision=$SOURCE_REV
sourceSha256=$SOURCE_SHA
binarySha256=$(sha256sum "$R/microprobe" | awk '{print $1}')
cubinSha256=$(sha256sum "$CUBIN" | awk '{print $1}')
sassSha256=$(sha256sum "$R/microprobe.sass" | awk '{print $1}')
sassAuditSha256=$(sha256sum "$R/sass-audit.tsv" | awk '{print $1}')
tablesSha256=$(sha256sum "$R/tables.bin" | awk '{print $1}')
samplesSha256=$(sha256sum "$R/samples.tsv" | awk '{print $1}')
EOF

/tmp/summarize-microprobe "$R/samples.tsv" "$R/gates.txt" "$R/bindings.txt" "$R/result.json"
cat "$R/result.json"

nvidia-smi --query-gpu=name,uuid,clocks.sm,power.draw,temperature.gpu --format=csv,noheader \
  > "$R/gpu-after.txt"
(
  cd "$R"
  find . -type f ! -name job.log ! -name exit-code ! -name scientific-files.sha256 \
    -print0 | sort -z | xargs -0 sha256sum > /tmp/sparse-microprobe-files.sha256
  mv /tmp/sparse-microprobe-files.sha256 scientific-files.sha256
)
echo "PASS: bounded native exact-layout microprobe complete"
