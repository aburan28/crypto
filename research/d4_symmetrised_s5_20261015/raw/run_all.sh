#!/bin/sh
set -u
B=/home/user/wt/d4-s5-bin/d4_weighted_s5
R=/home/user/wt/d4-s5/research/d4_symmetrised_s5_20261015/raw
export RAYON_NUM_THREADS=1
run() { nice -n 10 prlimit --as=2684354560 -- "$@"; }
date -u +%Y-%m-%dT%H:%MZ > $R/started_utc.txt
for p in 13 269; do run $B --invariance-control $p --seed 1 --json $R/controls.jsonl >> $R/controls.log 2>&1; done
for p in 269 1033 4093; do run $B --generic $p --seeds 3 --seed 1 --json $R/controls.jsonl >> $R/controls.log 2>&1; done
for p in 13 17 29 37; do run $B --p $p --seed 1 --targets 6 --brute --mitm --json $R/tiny.jsonl >> $R/tiny.log 2>&1; done
for p in 269 1033 4093; do run $B --p $p --seed 1 --targets 8 --mitm --json $R/measured.jsonl >> $R/measured.log 2>&1; done
for p in 269 1033 4093; do run $B --p $p --seed 2 --targets 12 --mitm --no-c8 --json $R/w_only.jsonl >> $R/w_only.log 2>&1; done
date -u +%Y-%m-%dT%H:%MZ > $R/finished_utc.txt
