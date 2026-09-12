#!/usr/bin/env bash
# test.sh -- exercise ecc2k130-fpga through the worker.py contract on the
# software model: selftest, a fresh run at a loose dp weight with the exact
# flags worker.py passes, resume from its checkpoint, a refused checkpoint
# (exit 6), a corpus reload, and the dp file's acceptance by merge.py.
set -euo pipefail
BIN=${1:-./ecc2k130-fpga}
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/ecc2k130-fpga-test.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

echo "== selftest"
"$BIN" --selftest

# worker.py's command line, plus the model and a stop condition
COMMON=(--curve 131 --steps 1024 --launches 0 --dp-file "$WORK/dp.bin" --checkpoint "$WORK/walk.ck"
        --checkpoint-every 600 --packed --device 0 --threads 192512 --dp-weight 56
        --max-iters 1073741824 --load-max 2000000 --dp-cap 65536
        --sim --sim-engines 2 --sim-idw 4 --sim-steps 4)

echo "== fresh run"
"$BIN" "${COMMON[@]}" --run-id 7 --verify 8 --dps 500 | tee "$WORK/run1.log"
grep -q "verified against the reference" "$WORK/run1.log"
n1=$(( $(stat -c %s "$WORK/dp.bin" 2>/dev/null || stat -f %z "$WORK/dp.bin") / 32 ))
test "$n1" -gt 0
test "$(head -c 8 "$WORK/walk.ck")" = "ECC2K130"
echo "   $n1 records"

echo "== resume"
"$BIN" "${COMMON[@]}" --run-id 7 --verify 2 --dps 100 | tee "$WORK/run2.log"
grep -q "^resumed from" "$WORK/run2.log"
grep -q "^reloaded $n1 points" "$WORK/run2.log"
n2=$(( $(stat -c %s "$WORK/dp.bin" 2>/dev/null || stat -f %z "$WORK/dp.bin") / 32 ))
test "$n2" -gt "$n1"

echo "== seeds never repeat across the resume"
python3 - "$WORK/dp.bin" <<'EOF'
import struct, sys
d = open(sys.argv[1], "rb").read()
seeds = [struct.unpack_from("<Q", d, i)[0] for i in range(0, len(d), 32)]
assert len(seeds) == len(set(seeds)), "duplicate seed in the corpus"
assert all(s >> 48 == 7 for s in seeds), "run id not in the seed"
print("   %d distinct seeds" % len(seeds))
EOF

echo "== refused checkpoint (other run id) -> exit 6"
set +e
"$BIN" "${COMMON[@]}" --run-id 8 --dps 1 >/dev/null
rc=$?
set -e
test "$rc" -eq 6

echo "== wrong dp weight for the bitstream -> exit 7"
set +e
"$BIN" --sim --sim-dp-weight 30 --launches 1 >/dev/null
rc=$?
set -e
test "$rc" -eq 7

echo "== merge.py accepts the corpus"
mkdir -p "$WORK/corpus" "$WORK/merge"
cp "$WORK/dp.bin" "$WORK/corpus/1-0000000000000000.bin"
if python3 -c "import numpy" 2>/dev/null; then
    python3 "$HERE/../../../ecc2k130/aws/merge.py" --work "$WORK/merge" --local "$WORK/corpus" --detect-only \
        | grep -q "\"corpus\": $n2"
else
    echo "   (numpy not installed, skipped)"
fi

echo "test.sh: PASS"
