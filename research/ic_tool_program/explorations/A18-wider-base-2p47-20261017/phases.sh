#!/bin/bash
# phases.sh RUNS: per arm, the geometric mean over its processes of the
# median repetition's cold time, collection and build, in ms, and the rest
# (cold minus the two).
set -eu
RUNS=$1
printf "arm\tprocesses\tcold_ms\tcollect_ms\tbuild_ms\trest_ms\n"
for a in v8-c320 v8-c400 v8-c470 a17-c320 a17-c400 a17-c470; do
  for f in $RUNS/$a/k0n61/*/r*.price.json; do
    jq -r '[([.repetitions[] | (.setup_ns + .ic_online.wall_ns)] | sort | .[1]), ([.repetitions[].setup_phases_ns.collect] | sort | .[1]), ([.repetitions[].setup_phases_ns.build] | sort | .[1])] | @tsv' $f
  done | awk -v a=$a '{c += log($1); l += log($2); b += log($3); n++} END {C = exp(c/n)/1e6; L = exp(l/n)/1e6; B = exp(b/n)/1e6; printf "%s\t%d\t%.1f\t%.1f\t%.1f\t%.1f\n", a, n, C, L, B, C - L - B}'
done
