#!/bin/bash
# Automatic sweep of the table-walk kernel's knobs and geometry on whatever
# GPU the container has (AUTOSWEEP.md).  Four phases, all on the card:
#
#   1. star     build the base point and every arm (base + one change), and
#               screen each with a short bench, twice, alternating;
#   2. greedy   take the arms that beat the base by MIN_GAIN, best first, and
#               add them to the base one at a time, keeping each that still
#               improves the screened rate -- the combination is measured, not
#               assumed, because knobs that share a pipe do not add;
#   3. verify   the base, the best single arm and the greedy result: 300
#               device reports re-walked by the host reference, and the whole
#               distinguished-point set on a forced common walk count compared
#               by hash across the three (a knob may not change the walk);
#   4. final    the same three binaries alternating for REPS repetitions of
#               the full bench command, SM clock and power sampled after each.
#
# The base point and the arm list are the tree's knowledge of the kernel; the
# ranking is the card's.  Static counts choose nothing here: on a logic-pipe-
# bound part a knob's sign is not the sign it has on a carry-less-bound one
# (TWO-CHAINS.md section 6), and only the card knows which part it is.
#
#   modal run --detach modal_job.py --job benchmarks/autosweep/gpujob.sh --out DIR --gpu B200
#   BASE_KNOBS / ARMS_FILE / MIN_GAIN / QUICK_LAUNCHES / REPS in the environment.
set -uo pipefail
cd "${WORK:-/work}" || exit 1
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq make g++ python3 >/dev/null 2>&1
R=${RESULTS:-/results}; mkdir -p "$R"
REPS=${REPS:-5}; QUICK=${QUICK_LAUNCHES:-32}; GAIN=${MIN_GAIN:-1.005}; MAX_GREEDY=${MAX_GREEDY:-8}
{
  nvidia-smi --query-gpu=name,driver_version,clocks.max.sm,power.limit,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -2
  echo "source: ${SOURCE_REV:-unknown}"
} | tee "$R/host.txt"
CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '. ')
ARCHFLAG="PRO6000_ARCH=-gencode arch=compute_${CAP},code=sm_${CAP}"
echo "arch: $ARCHFLAG" | tee -a "$R/host.txt"

# The base point: the 20 B/s knob set plus the B200 trades of TWO-CHAINS.md 6.4.
# On an RTX PRO 6000 pass BASE_KNOBS="" to sweep around the 20 B/s build itself.
BASE_KNOBS=${BASE_KNOBS-"PACKED_TOP_CLMAD=1 PACKED_ALU_SQUARE=0 PACKED_ONB_INV=1"}
# Arms: name|knob overrides appended to the base (last assignment wins).
if [ -n "${ARMS_FILE:-}" ] && [ -f "$ARMS_FILE" ]; then
  mapfile -t ARMS < <(grep -vE '^\s*(#|$)' "$ARMS_FILE")
else
  ARMS=(
    "geo-256x2|THREADS=256 MINBLOCKS=2"
    "geo-384x1|THREADS=384 MINBLOCKS=1"
    "geo-640x1|THREADS=640 MINBLOCKS=1"
    "geo-768x1|THREADS=768 MINBLOCKS=1"
    "batch32-256x1|BATCH=32 THREADS=256 MINBLOCKS=1"
    "batch32-512x1|BATCH=32 THREADS=512 MINBLOCKS=1"
    "karat3-not-topclmad|PACKED_KARAT3=1 PACKED_TOP_CLMAD=0"
    "no-pipe-select|TABLE_PIPE_SELECT=0"
    "no-chain-first|PACKED_CHAIN_FIRST=0"
    "inline-poly-0|PACKED_INLINE_POLY=0"
    "inline-poly-1|PACKED_INLINE_POLY=1"
    "unroll2|UNROLL_SLOTS=2"
    "no-from-reduced|PACKED_FROM_REDUCED=0"
    "no-l2-persist|PACKED_L2_PERSIST=0"
    "pivot-nibble|TABLE_PIVOT_BYTES=0"
    "clmul-flat|PACKED_CLMUL_FLAT=1"
    "fused|TABLE_FUSED=1 TABLE_PIPE_SELECT=0"
    "fused-pipe|TABLE_FUSED=1 TABLE_FUSED_PIPE=1 TABLE_PIPE_SELECT=0"
    "alu-sqr|PACKED_ALU_SQR=1"
    "no-pair-ilp|PACKED_PAIR_ILP=0"
  )
fi
echo "base: $BASE_KNOBS" | tee -a "$R/host.txt"
printf '%s\n' "${ARMS[@]}" > "$R/arms.txt"

declare -A KNOBS REGS SPILL RATE1 RATE2
fail=0
build() {  # build NAME KNOBS...
  local name=$1; shift
  KNOBS[$name]="$*"
  echo "=== build $name: $*"
  make -s gpu-preset "$ARCHFLAG" KNOBS="$BASE_KNOBS $*" > "$R/build-$name.log" 2>&1
  local info; info=$(grep -A2 "eccPacked131.*walk" "$R/build-$name.log" | grep -E "registers|spill" | head -2 | tr '\n' ' ')
  echo "    $info" | tee "$R/build-$name.txt"
  REGS[$name]=$(echo "$info" | grep -oE "Used [0-9]+ registers" | grep -oE "[0-9]+" | head -1)
  SPILL[$name]=$(echo "$info" | grep -oE "[0-9]+ bytes spill stores" | grep -oE "[0-9]+" | head -1)
  if [ ! -x ecc2k130 ]; then
    echo "BUILD FAILED $name" | tee -a "$R/failures.txt"; grep -E "error" "$R/build-$name.log" | head -3; return 1
  fi
  mv ecc2k130 "ecc2k130-$name"
}
rate() {  # rate NAME LAUNCHES -> B/s on stdout, clock/power line appended to bench.txt
  local b=$1 launches=$2
  ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches "$launches" --verify 0 2>&1 \
      | grep -oE "finished: [0-9.]+ M it/s" | grep -oE "[0-9.]+" | awk '{printf "%.4f", $1/1000}'
}
tick() { nvidia-smi --query-gpu=clocks.sm,power.draw,temperature.gpu --format=csv,noheader; }

echo "##### phase 1: star"
build base || { echo "base build failed; nothing to sweep"; exit 1; }
NAMES=(base)
for arm in "${ARMS[@]}"; do
  name=${arm%%|*}; knobs=${arm#*|}
  build "$name" $knobs && NAMES+=("$name")
done
for pass in 1 2; do
  for b in "${NAMES[@]}"; do
    r=$(rate "$b" "$QUICK"); [ -z "$r" ] && r=0
    if [ "$pass" = 1 ]; then RATE1[$b]=$r; else RATE2[$b]=$r; fi
    echo "screen $b pass $pass: $r B/s  $(tick)" | tee -a "$R/screen.txt"
  done
done
# Rank arms by their screened median against the base's, same passes.
python3 - "$R" "${NAMES[@]}" <<'PY' > "$R/star.txt"
import sys, re, statistics
r = sys.argv[1]; names = sys.argv[2:]
rates = {}
for line in open(r + "/screen.txt"):
    m = re.match(r"screen (\S+) pass (\d): ([0-9.]+) B/s", line)
    if m: rates.setdefault(m.group(1), []).append(float(m.group(3)))
base = statistics.median(rates["base"])
rows = []
for n in names:
    if n == "base" or n not in rates: continue
    med = statistics.median(rates[n]); rows.append((med / base, n, med))
rows.sort(reverse=True)
print("base %.4f" % base)
for ratio, n, med in rows: print("%.4f %s %.4f" % (ratio, n, med))
PY
cat "$R/star.txt"

echo "##### phase 2: greedy combination"
current="$BASE_KNOBS"; currentRate=$(python3 -c "import statistics;print(statistics.median([${RATE1[base]},${RATE2[base]}]))")
bestSingle=$(awk 'NR==2{print $2}' "$R/star.txt")
combo=0; comboName=base
while read -r ratio name med; do
  [ "$name" = base ] && continue
  awk -v r="$ratio" -v g="$GAIN" 'BEGIN{exit !(r+0 >= g+0)}' || break
  [ "$combo" -ge "$MAX_GREEDY" ] && break
  combo=$((combo + 1)); cand="combo$combo"
  knobs="${KNOBS[$name]}"
  build "$cand" "$(echo "${KNOBS[$comboName]:-} $knobs" | sed 's/^ *//')" || continue
  r1=$(rate "$cand" "$QUICK"); r2=$(rate "$cand" "$QUICK"); rc=$(python3 -c "import statistics;print(statistics.median([${r1:-0},${r2:-0}]))")
  echo "greedy $cand = $comboName + $name: $rc B/s vs $currentRate  $(tick)" | tee -a "$R/greedy.txt"
  if awk -v a="$rc" -v b="$currentRate" 'BEGIN{exit !(a+0 > b*1.002)}'; then
    echo "  keep" | tee -a "$R/greedy.txt"; comboName=$cand; currentRate=$rc
  else
    echo "  drop" | tee -a "$R/greedy.txt"
  fi
done < <(tail -n +2 "$R/star.txt")
echo "greedy result: $comboName (${KNOBS[$comboName]:-base}) at $currentRate B/s" | tee -a "$R/greedy.txt"

echo "##### phase 3: verify"
FINAL=(base)
[ -n "$bestSingle" ] && [ "$bestSingle" != base ] && FINAL+=("$bestSingle")
[ "$comboName" != base ] && [ "$comboName" != "$bestSingle" ] && FINAL+=("$comboName")
SMS=$(./ecc2k130-base --curve 131 --packed --bench --steps 1 --launches 1 --verify 0 2>&1 | grep -oE '[0-9]+ SMs' | grep -oE '[0-9]+' | head -1)
W=$((SMS * 512 * 16))
batchOf() { local k="$BASE_KNOBS ${KNOBS[$1]:-}"; echo "$k" | grep -oE "BATCH=[0-9]+" | tail -1 | grep -oE "[0-9]+" || echo 16; }
verify() {
  local b=$1 batch status=0; batch=$(batchOf "$b"); [ -z "$batch" ] && batch=16
  local threads=$((W / batch))
  if [ $((W % (batch * 256))) -ne 0 ]; then echo "=== verify $b: batch $batch does not divide the common walk count; identity skipped"; threads=0; fi
  echo "=== verify $b (batch $batch, threads $threads)"
  ./ecc2k130-$b --curve 131 --packed ${threads:+--threads $threads} --dp-weight 48 --dp-cap 262144 \
      --steps 96 --launches 6 --verify 300 --run-id 7 --dp-file "$R/dp-$b.bin" > "$R/verify-$b.log" 2>&1 || status=$?
  grep -E "MISMATCH|finished|resident|registers|OVERFLOW|packed (alu square|top clmad|onb inv|chains|table pivot)" "$R/verify-$b.log" | tee "$R/verify-$b.txt" || true
  return "$status"
}
for b in "${FINAL[@]}"; do verify "$b" || fail=1; done
python3 - "$R" "${FINAL[@]}" <<'PY' | tee "$R/dp-identity.txt" || fail=1
import hashlib, os, sys
r = sys.argv[1]; ref = None; bad = 0
for name in sys.argv[2:]:
    path = os.path.join(r, "dp-%s.bin" % name)
    if not os.path.exists(path): print("%-20s missing" % name); bad = 1; continue
    data = open(path, "rb").read()
    recs = sorted(data[i:i+32] for i in range(0, len(data) - len(data) % 32, 32))
    digest = hashlib.sha256(b"".join(recs)).hexdigest()
    if ref is None: ref = digest
    if digest != ref: bad = 1
    print("%-20s %8d records  sha256 %s  %s" % (name, len(recs), digest[:16], "IDENTICAL to base" if digest == ref else "DIFFERS from base"))
sys.exit(bad)
PY

echo "##### phase 4: final alternating bench"
for rep in $(seq 1 "$REPS"); do
  for b in "${FINAL[@]}"; do
    echo "=== bench $b rep $rep" | tee -a "$R/bench.txt"
    { ./ecc2k130-$b --curve 131 --packed --bench --steps 1024 --launches 64 --verify 0 2>&1 | grep -E "finished|resident" | tr '\n' ' '; tick; } | tee -a "$R/bench.txt"
  done
done
# Everything in one JSON for the summariser.
python3 - "$R" "${NAMES[@]}" <<'PY'
import json, os, re, sys, statistics
r = sys.argv[1]; names = sys.argv[2:]
knobs = {}
for n in names + [l.split()[1] for l in open(r + "/greedy.txt") if l.startswith("greedy combo")]:
    log = os.path.join(r, "build-%s.log" % n)
    txt = open(r + "/build-%s.txt" % n).read() if os.path.exists(r + "/build-%s.txt" % n) else ""
    regs = re.search(r"Used (\d+) registers", txt); spill = re.search(r"(\d+) bytes spill stores", txt)
    knobs[n] = {"registers": int(regs.group(1)) if regs else None, "spillBytes": int(spill.group(1)) if spill else None}
screen = {}
for line in open(r + "/screen.txt"):
    m = re.match(r"screen (\S+) pass (\d): ([0-9.]+) B/s\s+(\d+) MHz, ([0-9.]+) W, (\d+)", line)
    if m: screen.setdefault(m.group(1), []).append({"rateBps": float(m.group(3)), "smClockMHz": int(m.group(4)), "powerW": float(m.group(5)), "tempC": int(m.group(6))})
greedy = [l.rstrip() for l in open(r + "/greedy.txt")]
json.dump({"arms": open(r + "/arms.txt").read().splitlines(), "build": knobs, "screen": screen,
           "star": open(r + "/star.txt").read().splitlines(), "greedy": greedy}, open(r + "/sweep.json", "w"), indent=2)
PY
echo "=== done"
exit "$fail"
