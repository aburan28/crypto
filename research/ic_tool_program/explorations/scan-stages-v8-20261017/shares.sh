#!/bin/bash
# shares.sh RUNS: per curve, over its clean processes, the mean of: the
# scan's probed cycles (over the counter's rate) as a share of the
# collection phase, and the collection, build and selection phases as
# shares of cold time (set-up plus online), each summed over a process's
# repetitions.
set -eu
RUNS=$1
declare -A SLUG=([k0n53]=icv1-f2m53-tm56619371-dac20a85 [k1n59]=icv1-f2m59-tm943548413-98844ecc [k0n61]=icv1-f2m61-t158598901-ab42b6c5)
printf "curve\tprocesses\tscan_of_collection\tcollection_of_cold\tbuild_of_cold\tselection_of_cold\n"
for s in k0n53 k1n59 k0n61; do
  for f in $RUNS/$s/*/r*.price.json; do
    grep -q '"contended": false' ${f%.price.json}.isolation.jsonl || continue
    jq -r '(.scan_probes as $p | ([$p.stages[]] | add) / $p.tsc_hz * 1e9) as $scan
      | ([.repetitions[].setup_phases_ns.collect] | add) as $col
      | ([.repetitions[] | .setup_ns + .ic_online.wall_ns] | add) as $cold
      | "\($scan / $col) \($col / $cold) \(([.repetitions[].setup_phases_ns.build] | add) / $cold) \(([.repetitions[].setup_phases_ns.select] | add) / $cold)"' $f
  done | awk -v c=${SLUG[$s]} '{a+=$1; b+=$2; d+=$3; e+=$4; n++} END {printf "%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n", c, n, a/n, b/n, d/n, e/n}'
done
