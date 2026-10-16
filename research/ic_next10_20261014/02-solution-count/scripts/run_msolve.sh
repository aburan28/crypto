#!/usr/bin/env bash
# msolve validation arm: F̄_p solution count = degree of the eliminating polynomial of msolve's
# rational parametrization (-P 1). Targets k = 0..4. Full systems: m = 2 (t <= 6), m = 3 (t <= 4).
# Partial systems (one membership equation dropped): m = 2 (t <= 6), m = 3 (t <= 3).
# Each run: single thread, 700 MB resident cap (scripts/rsscap.sh), 300 s timeout.
set -u
V=$1; out=$V/msolve_counts.csv
echo "file,arm,m,t,k,kind,status,degree,seconds,peak_rss_kb" > $out
for f in $V/msolve_in/*.ms; do
  b=$(basename $f .ms)
  arm=${b%%_m*}; rest=${b#${arm}_m}; m=${rest%%_*}; rest=${rest#*_t}; t=${rest%%_*}; rest=${rest#*_k}; k=${rest%%_*}; kind=${rest#*_}
  [ $k -le 4 ] || continue
  if [ $kind = full ]; then { [ $m = 2 ] && [ $t -le 6 ]; } || { [ $m = 3 ] && [ $t -le 4 ]; } || continue; fi
  if [ $kind = partial ]; then { [ $m = 2 ] && [ $t -le 6 ]; } || { [ $m = 3 ] && [ $t -le 3 ]; } || continue; fi
  o=$V/msolve_out/$b.out; mkdir -p $V/msolve_out
  s=$(date +%s.%N)
  peak=$(OMP_NUM_THREADS=1 nice $(dirname "$0")/rsscap.sh 300 msolve -t 1 -P 1 -f $f -o $o 3>&1 >/dev/null 2>&1); rc=$?
  el=$(echo "$(date +%s.%N) - $s" | bc)
  flat=$(tr -d '\n ' < $o 2>/dev/null)
  if [ $rc -ne 0 ]; then st="fail_rc$rc"; deg=NA
  elif [[ $flat == "[-1]"* ]]; then st=ok; deg=0
  elif [[ $flat == "[0,["* ]]; then st=ok; deg=$(echo "$flat" | cut -d, -f4)
  else st="other"; deg=NA; fi
  echo "$b,$arm,$m,$t,$k,$kind,$st,$deg,$el,$peak" >> $out
done
