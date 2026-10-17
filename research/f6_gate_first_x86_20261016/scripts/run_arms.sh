#!/bin/sh
# Gate-first / RREF-root confirmation on x86 (BMI2, eager rewrite default).
# One binary (f6_ic_probe built from the commit in the report), arms set by
# environment only: off (no switch), gate (KIC_F6_GATE_FIRST=1), rref
# (KIC_F4_INHERIT_ROOT=rref), gaterref (both). The four cells of
# research/lit2026_w5_lane_batched_f4_20261011 plus five new ones.
#   PROBE  path to f6_ic_probe    OUT  output directory    ROUNDS (default 5)
PROBE=${PROBE:?set PROBE}
OUT=${OUT:?set OUT}
ROUNDS=${ROUNDS:-5}
mkdir -p "$OUT/raw"
cell_flags() {
  case $1 in
    K1_17) echo "--a 1 --n 17 --dim 6 --targets 8" ;;
    K1_23) echo "--a 1 --n 23 --dim 7 --targets 4" ;;
    K0_23) echo "--a 0 --n 23 --dim 8 --targets 4" ;;
    K0_31) echo "--a 0 --n 31 --dim 9 --targets 1 --budget 2000" ;;
    K0_19) echo "--a 0 --n 19 --dim 7 --targets 4" ;;
    K1_19) echo "--a 1 --n 19 --dim 6 --targets 8" ;;
    K1_17d7) echo "--a 1 --n 17 --dim 7 --targets 8" ;;
    K1_29) echo "--a 1 --n 29 --dim 8 --targets 2" ;;
    K0_37) echo "--a 0 --n 37 --dim 9 --targets 1 --budget 2000" ;;
  esac
}
run_arm() {
  arm=$1; cell=$2; round=$3
  file="$OUT/raw/r${round}-${cell}-${arm}.jsonl"
  echo "{\"uptime\":\"$(uptime | sed 's/.*load average: //')\",\"utc\":\"$(date -u +%Y-%m-%dT%H:%M:%SZ)\",\"arm\":\"$arm\",\"cell\":\"$cell\",\"round\":$round}" > "$file"
  case $arm in
    off)      env -u KIC_F6_GATE_FIRST -u KIC_F4_INHERIT_ROOT "$PROBE" $(cell_flags "$cell") --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    gate)     KIC_F6_GATE_FIRST=1 "$PROBE" $(cell_flags "$cell") --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    rref)     KIC_F4_INHERIT_ROOT=rref "$PROBE" $(cell_flags "$cell") --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    gaterref) KIC_F6_GATE_FIRST=1 KIC_F4_INHERIT_ROOT=rref "$PROBE" $(cell_flags "$cell") --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
  esac
  echo "exit=$?" > "$OUT/raw/r${round}-${cell}-${arm}.status"
}
ARMS="off gate rref gaterref"
round=1
while [ "$round" -le "$ROUNDS" ]; do
  for cell in K1_17 K1_23 K0_23 K0_31 K0_19 K1_19 K1_17d7 K1_29 K0_37; do
    order=$(echo $ARMS | awk -v r="$round" '{n=NF; s=""; for(i=0;i<n;i++){j=((i+r-1)%n)+1; s=s" "$j}; print s}')
    for arm in $order; do run_arm "$arm" "$cell" "$round"; done
  done
  round=$((round + 1))
done
echo done
