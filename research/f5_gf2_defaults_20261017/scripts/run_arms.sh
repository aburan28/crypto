#!/bin/sh
# Confirmation of the matrix-F5 / GF(2) defaults (direct pack, direct
# unpack, table reuse), one binary per tool, arms by environment only:
#   on        no switch set (the new defaults)
#   nopack    KIC_F5_DIRECT_PACK=0
#   nounpack  KIC_F5_UNPACK_DIRECT=0
#   noreuse   KIC_GF2_REUSE_TABLE=0
#   off       all three =0 (the previous defaults)
# f5: f4_f2_bench f5 cases (n <= 24), Reduced and SelectiveEchelon output,
#     RAYON_NUM_THREADS 1 and 4, all five arms.
# f4: f4_f2_bench f4 cases (n <= 20), 1 thread, arms on / noreuse.
# inh: f6_ic_probe inherited F4 (f4 engine) on three cells, arms on / noreuse.
#   BENCH PROBE OUT ROUNDS(default 3)
BENCH=${BENCH:?}; PROBE=${PROBE:?}; OUT=${OUT:?}; ROUNDS=${ROUNDS:-3}
mkdir -p "$OUT/raw"
armenv() {
  case $1 in
    on) echo "" ;;
    nopack) echo "KIC_F5_DIRECT_PACK=0" ;;
    nounpack) echo "KIC_F5_UNPACK_DIRECT=0" ;;
    noreuse) echo "KIC_GF2_REUSE_TABLE=0" ;;
    off) echo "KIC_F5_DIRECT_PACK=0 KIC_F5_UNPACK_DIRECT=0 KIC_GF2_REUSE_TABLE=0" ;;
  esac
}
stamp() { echo "{\"uptime\":\"$(uptime | sed 's/.*load average: //')\",\"utc\":\"$(date -u +%Y-%m-%dT%H:%M:%SZ)\",\"tool\":\"$1\",\"arm\":\"$2\",\"form\":\"$3\",\"threads\":$4,\"round\":$5}"; }
rot() { echo $1 | awk -v r="$2" '{n=NF; s=""; for(i=0;i<n;i++){j=((i+r-1)%n)+1; s=s" "$j}; print s}'; }
r=1
while [ "$r" -le "$ROUNDS" ]; do
  for th in 1 4; do
    for form in 0 2; do
      for arm in $(rot "on nopack nounpack noreuse off" $r); do
        f="$OUT/raw/r$r-f5-form$form-t$th-$arm.jsonl"
        stamp f5 $arm $form $th $r > "$f"
        env $(armenv $arm) KIC_F5_ECHELON=$form RAYON_NUM_THREADS=$th "$BENCH" 3 24 f5 0 >> "$f" 2>/dev/null
        echo "exit=$?" > "${f%.jsonl}.status"
      done
    done
  done
  for arm in $(rot "on noreuse" $r); do
    f="$OUT/raw/r$r-f4-t1-$arm.jsonl"
    stamp f4 $arm - 1 $r > "$f"
    env $(armenv $arm) RAYON_NUM_THREADS=1 "$BENCH" 3 20 f4 0 >> "$f" 2>/dev/null
    echo "exit=$?" > "${f%.jsonl}.status"
    for cell in "K1_23:--a 1 --n 23 --dim 7 --targets 4" "K0_23:--a 0 --n 23 --dim 8 --targets 4" "K0_31:--a 0 --n 31 --dim 9 --targets 1 --budget 2000"; do
      c=${cell%%:*}; flags=${cell#*:}
      f="$OUT/raw/r$r-inh-$c-$arm.jsonl"
      stamp inh $arm $c 1 $r > "$f"
      env $(armenv $arm) "$PROBE" $flags --engines f4 --reps 3 >> "$f" 2>/dev/null
      echo "exit=$?" > "${f%.jsonl}.status"
    done
  done
  r=$((r + 1))
done
echo done
