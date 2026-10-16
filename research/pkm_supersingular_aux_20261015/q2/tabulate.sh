#!/bin/bash
# Lay out the recorded Q2 rows (runs/*.jsonl) as one TSV, one line per system
# (field extraction only).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
{
  printf 'run\tarm\tcontrol\tp\tm\tt\tN\ttarget_index\tD\twidth_to_solution\tsteps_to_solution\tinconsistent\ttimed_out\toversize\tverify\tmuladds\tms\taux_a\taux_b\taux_order\n'
  for f in runs/*.jsonl; do
    r=$(basename "$f" .jsonl)
    jq -r --arg r "$r" '[$r, .arm, .control, .p, .m, .t, .N, .target_index, .solving_degree_max, .max_cols_to_solution, .steps_to_solution, .inconsistent, .timed_out, .oversize, (.verify // "none"), .muladds, .ms, (.tower.aux_a // ""), (.tower.aux_b // ""), (.tower.aux_order // "")] | @tsv' "$f"
  done
} > cells.tsv
