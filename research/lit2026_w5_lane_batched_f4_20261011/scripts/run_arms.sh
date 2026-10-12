#!/bin/sh
# lit2026 W5 interleaved arms on the four f6_ic_probe production cells.
#   OUT     output directory (raw/ is created under it)
#   ROUNDS  rounds (default 5); the arm order rotates every round
# Arms: orig (4b6e09b78 binary), prev (candidate binary, switches off),
# gate (KIC_F6_GATE_FIRST=1), lanes (KIC_F4_LANES=<targets of the cell>
# with KIC_F4_INHERIT_ROOT=rref), rref (KIC_F4_INHERIT_ROOT=rref, the
# sequential control the lane roots reproduce row for row), and gaterref
# (both switches).  ARMS selects and orders them.
BIN=/Volumes/SSD990-2/crypto/targets/lit2026-w5/bin
OUT=${OUT:?set OUT}
ROUNDS=${ROUNDS:-5}
mkdir -p "$OUT/raw"
cd /Volumes/SSD990-2/crypto/targets/lit2026-w5/scratch || exit 1

cell_flags() {
  case $1 in
    K1_17) echo "--a 1 --n 17 --dim 6 --targets 8" ;;
    K1_23) echo "--a 1 --n 23 --dim 7 --targets 4" ;;
    K0_23) echo "--a 0 --n 23 --dim 8 --targets 4" ;;
    K0_31) echo "--a 0 --n 31 --dim 9 --targets 1 --budget 2000" ;;
  esac
}
cell_lanes() {
  case $1 in
    K1_17) echo 8 ;;
    K1_23) echo 4 ;;
    K0_23) echo 4 ;;
    K0_31) echo 1 ;;
  esac
}

run_arm() {
  arm=$1; cell=$2; round=$3
  flags=$(cell_flags "$cell")
  file="$OUT/raw/r${round}-${cell}-${arm}.jsonl"
  echo "{\"uptime\":\"$(uptime | sed 's/.*load averages*: //')\",\"utc\":\"$(date -u +%Y-%m-%dT%H:%M:%SZ)\",\"arm\":\"$arm\",\"cell\":\"$cell\",\"round\":$round}" > "$file"
  case $arm in
    orig)  $BIN/f6_ic_probe.orig $flags --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    prev)  $BIN/f6_ic_probe.cand $flags --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    gate)  KIC_F6_GATE_FIRST=1 $BIN/f6_ic_probe.cand $flags --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    lanes) KIC_F4_INHERIT_ROOT=rref KIC_F4_LANES=$(cell_lanes "$cell") $BIN/f6_ic_probe.cand $flags --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    rref)  KIC_F4_INHERIT_ROOT=rref $BIN/f6_ic_probe.cand $flags --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    gaterref) KIC_F6_GATE_FIRST=1 KIC_F4_INHERIT_ROOT=rref $BIN/f6_ic_probe.cand $flags --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
  esac
  echo "exit=$?" >> "$OUT/raw/r${round}-${cell}-${arm}.status"
}

ARMS=${ARMS:-"orig prev gate lanes rref"}
round=1
while [ "$round" -le "$ROUNDS" ]; do
  for cell in K1_17 K1_23 K0_23 K0_31; do
    # rotate the arm order by round
    order=$(echo $ARMS | awk -v r="$round" '{n=NF; s=""; for(i=0;i<n;i++){j=((i+r-1)%n)+1; s=s" "$j}; print s}')
    for arm in $order; do
      run_arm "$arm" "$cell" "$round"
    done
  done
  round=$((round + 1))
done
echo done
