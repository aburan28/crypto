#!/bin/bash
# phases.sh RUNS: per curve and phase, the geometric mean over all pairs of
# v8's median repetition over the gate's, and each arm's median cold time.
set -eu
RUNS=$1
declare -A SLUG=([k0n37]=icv1-f2m37-tm534059-32aad96b [k0n57]=icv1-f2m57-tm747311035-c1f545af [k0n53]=icv1-f2m53-tm56619371-dac20a85)
med() { jq "[.repetitions[] | $2] | sort | .[length/2|floor]" "$1"; }
printf "curve\tphase\tpairs\tv8_over_gate\n"
for s in k0n37 k0n57 k0n53; do
  for ph in collect build select; do
    for f in $RUNS/v8/$s/*/r*.price.json; do
      echo "$(med $f .setup_phases_ns.$ph) $(med $RUNS/gate/${f#$RUNS/v8/} .setup_phases_ns.$ph)"
    done | awk -v c=${SLUG[$s]} -v p=$ph '{l += log($1/$2); n++} END {printf "%s\t%s\t%d\t%.3f\n", c, p, n, exp(l/n)}'
  done
done
