#!/bin/sh
# Interleaved before/after timing of the k = 4 Macaulay solve (generic
# planted systems, d = 4, 5, 6, two systems each, three rounds, ABAB).
# BEFORE / AFTER are the gaudry_quartic_replay binaries built from commits
# 5a58e9e98 and 029309579; raw/<arm>_r<round>_d<d>.jsonl is what each wrote.
set -e
out=$(dirname "$0")/raw
for r in 1 2 3; do
  for d in 4 5 6; do
    for arm in before after; do
      bin=$BEFORE; [ $arm = after ] && bin=$AFTER
      rm -f $out/${arm}_r${r}_d${d}.jsonl
      $bin --generic $d --systems 2 --json $out/${arm}_r${r}_d${d}.jsonl > /dev/null
    done
  done
done
