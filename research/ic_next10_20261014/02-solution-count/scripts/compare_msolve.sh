#!/usr/bin/env bash
# Join msolve F̄_p counts with the gp direct grid counts (same base, m, t, target).
# Prints: arm,m,t,k,kind,msolve_degree,grid_count,predicted,verdict
#   full:    predicted = grid count (split L), or "-" for random-nonsplit (no grid analogue)
#   partial: predicted = B^{m-1} * 2^{m-1} (degree of S_{m+1} in x_m), B = 2^t
V=$1
awk -F, 'NR==FNR { if (FNR>1) g[$1","$2","$3","$4]=$6; next }
  FNR>1 { key=$2","$3","$4","$5; B=2^$4;
    if ($6=="full") pred=(($2","$3","$4","$5) in g)? g[key] : "-";
    else pred=B^($3-1)*2^($3-1);
    v=($7!="ok")?$7:(pred=="-"?"n/a":($8==pred?"MATCH":"DIFF"));
    print $2","$3","$4","$5","$6","$8","(key in g?g[key]:"-")","pred","v }' $V/grid_counts.csv $V/msolve_counts.csv
