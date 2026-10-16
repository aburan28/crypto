#!/bin/bash
# pairs-tsv.sh RUNS A B: one row per (curve, row, round) both arms ran, with
# each arm's cold time (median repetition of set-up plus online, ns), the
# ratio A/B, each run's contended flag, its minor faults and its scalar.
# Curves are written by ICV1 slug; the run tree keeps the suite's names.
set -eu
RUNS=$1; A=$2; B=$3
slug() {
  case $1 in
    k0n53) echo icv1-f2m53-tm56619371-dac20a85 ;;
    k1n59) echo icv1-f2m59-tm943548413-98844ecc ;;
    k0n61) echo icv1-f2m61-t158598901-ab42b6c5 ;;
    *) echo "unknown size $1" >&2; exit 1 ;;
  esac
}
cold() { jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns)] | sort | .[length/2|floor]' "$1"; }
flag() { grep -o '"contended": [a-z]*' "$1" | tail -1 | awk '{print $2}'; }
faults() { grep -o '"minor_faults": [0-9]*' "$1" | tail -1 | awk '{print $2}'; }
printf "curve\trow\tround\t%s_cold_ns\t%s_cold_ns\t%s_over_%s\t%s_contended\t%s_contended\t%s_minor_faults\t%s_minor_faults\t%s_scalar\t%s_scalar\n" $A $B $A $B $A $B $A $B $A $B
for f in "$RUNS/$A"/*/*/r*.price.json; do
  rel=${f#$RUNS/$A/}; g=$RUNS/$B/$rel
  [ -s "$g" ] || continue
  size=${rel%%/*}; rest=${rel#*/}; row=${rest%%/*}; rnd=${rest#*/}; rnd=${rnd%.price.json}
  ca=$(cold "$f"); cb=$(cold "$g")
  printf "%s\t%s\t%s\t%s\t%s\t%.6f\t%s\t%s\t%s\t%s\t%s\t%s\n" "$(slug $size)" $row $rnd $ca $cb \
    "$(echo "$ca $cb" | awk '{print $1/$2}')" \
    "$(flag ${f%.price.json}.isolation.jsonl)" "$(flag ${g%.price.json}.isolation.jsonl)" \
    "$(faults ${f%.price.json}.isolation.jsonl)" "$(faults ${g%.price.json}.isolation.jsonl)" \
    "$(jq -r '.certificates.ic.scalar' "$f")" "$(jq -r '.certificates.ic.scalar' "$g")"
done
