#!/bin/sh
# Run 2: see DECISION_RULE.md.  BENCH OUT ROUNDS(15)
BENCH=${BENCH:?}; OUT=${OUT:?}; ROUNDS=${ROUNDS:-15}
mkdir -p "$OUT/raw"
armenv() {
  case $1 in
    on) echo "" ;;
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
      for arm in $(rot "on nounpack noreuse off" $r); do
        f="$OUT/raw/r$r-f5-form$form-t$th-$arm.jsonl"
        stamp f5 $arm $form $th $r > "$f"
        env $(armenv $arm) KIC_F5_ECHELON=$form RAYON_NUM_THREADS=$th "$BENCH" 3 24 f5 0 >> "$f" 2>/dev/null
        echo "exit=$?" > "${f%.jsonl}.status"
      done
    done
  done
  if [ "$r" -le 10 ]; then
    for arm in $(rot "on noreuse" $r); do
      f="$OUT/raw/r$r-f4-t1-$arm.jsonl"
      stamp f4 $arm - 1 $r > "$f"
      env $(armenv $arm) RAYON_NUM_THREADS=1 "$BENCH" 3 18 f4 0 >> "$f" 2>/dev/null
      echo "exit=$?" > "${f%.jsonl}.status"
    done
  fi
  r=$((r + 1))
done
echo done
