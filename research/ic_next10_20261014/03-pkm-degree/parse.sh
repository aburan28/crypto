#!/bin/bash
# One CSV line per msolve run: name,m,t,N,x_r,rounds,D_max_round,D_last_productive,
# max_rows,max_cols,verdict,VmHWM_kB,VmPeak_kB,wall_s,exit
cd "$(dirname "$0")"
echo "name,m,t,N,x_r,rounds,D_msolve,D_last_productive,max_rows,max_cols,verdict,VmHWM_kB,VmPeak_kB,wall_s,exit"
for log in runs/*.log; do
  name=$(basename $log .log)
  m=$(echo $name | sed -E 's/.*_m([0-9]+)_.*/\1/'); t=$(echo $name | sed -E 's/.*_t([0-9]+)_.*/\1/'); x=$(echo $name | sed -E 's/.*_x([0-9]+).*/\1/')
  stats=$(awk '/^ +[0-9]+ +[0-9]+ +[0-9]+ +[0-9]+ x [0-9]+/ {n++; d=$1; if(d>D)D=d; if($8>0 && d>P)P=d; if($4>R)R=$4; if($6>C)C=$6} END{printf "%d,%d,%d,%d,%d", n, D, P, R, C}' $log)
  if [ -s runs/$name.out ] && grep -q '^\[-1\]' runs/$name.out; then v=no_solution; elif [ -s runs/$name.out ]; then v=solutions; else v=none; fi
  mem=$(cat runs/$name.mem 2>/dev/null | sed -E 's/[A-Za-z_]+=//g' | tr ' ' ',')
  [ -z "$mem" ] && mem="running,,,"
  echo "$name,$m,$t,$((m*t)),$x,$stats,$v,$mem"
done
