#!/bin/bash
# pairs2.sh RUNS BASE CAND: per size, the geometric mean of base/cand cold
# time over clean pairs (both runs uncontended), with a 95% t interval;
# and the same over all pairs; and whether every pair's scalars agree.
RUNS=$1; A=$2; B=$3
for f in $RUNS/$A/*/*/r*.price.json; do
  rel=${f#$RUNS/$A/}; g=$RUNS/$B/$rel
  [ -s "$g" ] || { echo "missing $rel" >&2; continue; }
  ca=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns)] | sort | .[length/2|floor]' $f)
  cb=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns)] | sort | .[length/2|floor]' $g)
  sa=$(jq -r '.certificates.ic.scalar' $f); sb=$(jq -r '.certificates.ic.scalar' $g)
  xa=$(grep -o '"contended": [a-z]*' ${f%.price.json}.isolation.jsonl | tail -1 | awk '{print $2}')
  xb=$(grep -o '"contended": [a-z]*' ${g%.price.json}.isolation.jsonl | tail -1 | awk '{print $2}')
  echo "${rel%%/*} $(echo "$ca $cb" | awk '{print $1/$2}') $xa $xb $([ "$sa" = "$sb" ] && echo same || echo DIFF)"
done | awk '
function tq(df) { return df==1?12.706:df==2?4.303:df==3?3.182:df==4?2.776:df==5?2.571:df==6?2.447:df==7?2.365:df==8?2.306:df==9?2.262:df==10?2.228:df==11?2.201:df==12?2.179:df==13?2.160:df==14?2.145:2.131 }
{ all[$1] = all[$1] " " log($2); if ($3 == "false" && $4 == "false") cl[$1] = cl[$1] " " log($2); if ($5 != "same") bad[$1]++ }
END {
  for (s in all) for (pass = 0; pass < 2; pass++) {
    str = pass ? cl[s] : all[s]; n = split(str, x, " "); if (n < 2) continue; m = 0; for (i = 1; i <= n; i++) m += x[i]; m /= n;
    ss = 0; for (i = 1; i <= n; i++) ss += (x[i]-m)^2; sd = sqrt(ss/(n-1)); h = tq(n-1)*sd/sqrt(n);
    printf "%s %-5s n=%2d %.4f [%.4f, %.4f]%s\n", s, pass ? "clean" : "all", n, exp(m), exp(m-h), exp(m+h), bad[s] ? " SCALAR MISMATCH" : ""
  }
}' | sort
