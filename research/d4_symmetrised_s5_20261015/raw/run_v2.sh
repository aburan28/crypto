#!/bin/sh
set -u
B=/home/user/wt/d4-s5-bin/v2/d4_weighted_s5
R=/home/user/wt/d4-s5/research/d4_symmetrised_s5_20261015/raw
export RAYON_NUM_THREADS=1
run() { nice -n 10 prlimit --as=2684354560 -- "$@"; }
date -u +%Y-%m-%dT%H:%MZ > $R/v2_started_utc.txt
for p in 13 17 29 37; do run $B --p $p --seed 1 --targets 6 --brute --mitm --json $R/v2_tiny.jsonl >> $R/v2_tiny.log 2>&1; done
for p in 269 1033 4093; do run $B --p $p --seed 1 --targets 8 --mitm --json $R/v2_measured.jsonl >> $R/v2_measured.log 2>&1; done
for p in 269 1033 4093; do run $B --p $p --seed 2 --targets 12 --mitm --no-c8 --json $R/v2_w_only.jsonl >> $R/v2_w_only.log 2>&1; done
date -u +%Y-%m-%dT%H:%MZ > $R/v2_finished_utc.txt
