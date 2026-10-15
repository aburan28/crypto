#!/bin/bash
# Lay out cells.tsv and recorded_split.tsv as the README's Markdown table rows
# (field extraction and joining only).
set -eu
cd "$(dirname "$0")"
cell() { # run-prefix m N -> "D / width / steps" (both targets equal), else each target; "(1)" for one target
  awk -F'\t' -v r="$1" -v m="$2" -v n="$3" 'NR>1 && index($1,r)==1 && $4==m && $5==n {
     v=$8" / "$9" / "$10; if($12=="true"||$13=="true") v=v" (stopped)"; c++; s[v]++; o=o (o?"; ":"") v }
     END { if(!c){print "—"; exit} k=0; for(x in s)k++; if(k==1){for(x in s) print x (c==1?" (1)":"")} else print o }' cells.tsv; }
rec() { awk -F'\t' -v m="$1" -v n="$2" 'NR>1 && $4==m && $5==n && ($3=="2013265921"||$3=="786433") {print $6" / "$7" / "$8" ("$1")"; exit}' recorded_split.tsv; }
recp() { awk -F'\t' -v m="$1" -v n="$2" -v p="$3" 'NR>1 && $4==m && $5==n && $3==p {print $6" / "$7" / "$8" ("$9")"; exit}' recorded_split.tsv; }
echo "<!-- m2 -->"
for n in 2 4 6 8 10 12 14 16 18 20; do
  echo "| $n | $(cell R1- 2 $n) | $(cell R2- 2 $n) | $(cell R3- 2 $n) | $(rec 2 $n) | $(cell R4- 2 $n) | $(cell R5- 2 $n) |"
done
echo "<!-- m34 -->"
for mn in "3 6" "3 9" "3 12" "3 15" "4 8" "4 12"; do set -- $mn
  r12=$(cell R12- 3 $2); r13=$(cell R13- 3 $2); t31=$(cell R6- $1 $2); t24=$(cell R7- $1 $2); sp=$(cell R8- $1 $2)
  [ "$2" = 15 ] && { t31=$r12; t24=$r13; sp=$(cell R14- 3 15); }
  echo "| $1 | $2 | $t31 | $t24 | $sp | $(rec $1 $2) | $(cell R9- $1 $2) | $(cell R10- $1 $2) |"
done
echo "<!-- deep -->"
for mn in "2 2" "2 4" "2 6" "2 8" "2 10" "2 12" "2 14" "2 16" "2 18" "3 6" "3 9" "3 12"; do set -- $mn
  echo "| $1 | $2 | $(cell R11- $1 $2) | $(recp $1 $2 2013265921) | $(recp $1 $2 786433) | $(cell R11n- $1 $2) |"
done
