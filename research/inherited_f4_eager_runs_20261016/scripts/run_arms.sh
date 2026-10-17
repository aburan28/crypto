#!/bin/sh
# Portable eager rewrite (shift runs) against BMI2 eager, the lazy rewrite and the bit scatter (the current non-BMI2 path).
# One binary, arms set by environment only: bmi2 (default on this host), runs (KIC_F4_EAGER=runs),
# lazy (KIC_F4_EAGER=0), scatter (KIC_F4_REWRITE=scatter: lazy and bit scatter, what an aarch64 host ran).
# Cells as in research/f6_gate_first_x86_20261016.
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
    bmi2)    env -u KIC_F4_EAGER -u KIC_F4_REWRITE "$PROBE" $(cell_flags "$cell") --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    runs)    KIC_F4_EAGER=runs "$PROBE" $(cell_flags "$cell") --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    lazy)    KIC_F4_EAGER=0 "$PROBE" $(cell_flags "$cell") --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
    scatter) KIC_F4_REWRITE=scatter "$PROBE" $(cell_flags "$cell") --engines f4,f6p --reps 3 >> "$file" 2>/dev/null ;;
  esac
  echo "exit=$?" > "$OUT/raw/r${round}-${cell}-${arm}.status"
}
ARMS="bmi2 runs lazy scatter"
round=1
while [ "$round" -le "$ROUNDS" ]; do
  for cell in K1_17 K1_23 K0_23 K0_31 K0_19 K1_19 K1_17d7 K1_29 K0_37; do
    order=$(echo $ARMS | awk -v r="$round" '{n=NF; s=""; for(i=0;i<n;i++){j=((i+r-1)%n)+1; s=s" "$j}; print s}')
    for arm in $order; do run_arm "$arm" "$cell" "$round"; done
  done
  round=$((round + 1))
done
echo done
