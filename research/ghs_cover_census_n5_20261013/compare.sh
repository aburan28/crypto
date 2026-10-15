#!/bin/sh
# Reproduction checks R1-R3 (PROTOCOL.md): jq projections compared byte-for-byte with cmp.
# Usage: sh compare.sh <run-id> <scratch-dir>   (from the study directory). Prints one PASS/FAIL line per check.
RUN=$1; T=$2; O=results/run-$RUN; M=../ghs_cubic_cover_census_20261012/results/run-20261013T0403Z.jsonl
mkdir -p "$T"
F='{A_idx,B_idx,order,maxprime,min_genus,types_at_min,structures_by_genus,m1_structures,min_genus_2isog,min_genus_class}'
F2='{A_idx,B_idx,order,maxprime,min_genus,types_at_min,structures_by_genus,m1_structures,min_genus_2isog}'
chk() { if cmp -s "$1" "$2"; then echo "PASS $3 ($(wc -l < "$1" | tr -d ' ') rows)"; else echo "FAIL $3"; diff "$1" "$2" | head -5; fi; }
for q in 5 7; do
  jq -c "select(.rec==\"curve\" and .q==$q) | $F" "$M" > "$T/merged-$q.txt"
  jq -c "select(.rec==\"curve\") | $F" "$O/R1-n3-q$q.jsonl" > "$T/new-$q.txt"
  chk "$T/merged-$q.txt" "$T/new-$q.txt" "R1 curve records q=$q"
  jq -c "select(.rec==\"mod4\" and .q==$q) | {scope,order_mod4,curves_by_min_genus}" "$M" | sort > "$T/m4m-$q.txt"
  jq -c "select(.rec==\"group\" and .by==\"order_mod4\") | {scope,order_mod4:.key,curves_by_min_genus}" "$O/R1-n3-q$q.jsonl" | sort > "$T/m4n-$q.txt"
  chk "$T/m4m-$q.txt" "$T/m4n-$q.txt" "R1 mod4 histograms q=$q"
done
jq -c "select(.rec==\"curve\" and .q==11) | $F2" "$M" > "$T/merged-11.txt"
jq -c "select(.rec==\"curve\" and .role==\"sample\") | $F2" "$O/R2-n3-q11-sample120.jsonl" > "$T/new-11.txt"
chk "$T/merged-11.txt" "$T/new-11.txt" "R2 sampled curve records q=11"
jq -c 'select(.rec=="curve") | del(.mismatch)' "$O/D-n5-q3-cosets.jsonl" > "$T/r3a.txt"
jq -c 'select(.rec=="curve") | del(.mismatch)' "$O/R3-n5-q3-cocycles.jsonl" > "$T/r3b.txt"
chk "$T/r3a.txt" "$T/r3b.txt" "R3 cosets vs cocycles n=5 q=3"
jq -c 'select(.rec=="curve") | del(.mismatch)' "$O/pilot-n3-q3-cosets.jsonl" > "$T/p3a.txt"
jq -c 'select(.rec=="curve") | del(.mismatch)' "$O/pilot-n3-q3-cocycles.jsonl" > "$T/p3b.txt"
chk "$T/p3a.txt" "$T/p3b.txt" "pilot cosets vs cocycles n=3 q=3"
