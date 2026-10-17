#!/bin/bash
# phases.sh RUNS ARM: per curve and phase, the geometric mean over all pairs
# of `off`'s median repetition over ARM's.
set -eu
RUNS=$1; B=$2
declare -A SLUG=([k0n53]=icv1-f2m53-tm56619371-dac20a85 [k1n59]=icv1-f2m59-tm943548413-98844ecc [k0n61]=icv1-f2m61-t158598901-ab42b6c5)
med() { jq "[.repetitions[] | $2] | sort | .[length/2|floor]" "$1"; }
printf "curve\tphase\tpairs\toff_over_%s\n" $B
for s in k0n53 k1n59 k0n61; do
  for ph in collect build select; do
    for f in $RUNS/off/$s/*/r*.price.json; do
      echo "$(med $f .setup_phases_ns.$ph) $(med $RUNS/$B/${f#$RUNS/off/} .setup_phases_ns.$ph)"
    done | awk -v c=${SLUG[$s]} -v p=$ph '{l += log($1/$2); n++} END {printf "%s\t%s\t%d\t%.3f\n", c, p, n, exp(l/n)}'
  done
done
