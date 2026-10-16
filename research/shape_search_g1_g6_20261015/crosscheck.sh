#!/bin/bash
# Re-run a subset with the lane-private build and compare with the committed records (secs stripped).
export F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000 KIC_SPARSE_DENSE_FINISH=1 KIC_SPARSE_F5=1 KIC_SPARSE_DENSE_BUDGET_MB=1200
B=/home/user/wt/shape-search-bin/own/factor_base_shape_ladder
R=/home/user/wt/shape-search/research/shape_search_g1_g6_20261015/runs
X=/tmp/claude-0/-home-user/4139a275-8aea-521c-b033-3fb42ac6f932/scratchpad/xcheck; mkdir -p $X
S=20261015
strip() { sed -E 's/,"secs":[0-9.]+//' "$1"; }
run() { name=$1; shift; nice -n 10 prlimit --as=2684354560 -- $B "$@" > $X/$name; }
run g1-5-2.jsonl g1 --cell 5:2:7 --shapes random,random2,gp,fp,coord --instances 4 --targets 2 --seed $S
run g1-7-3.jsonl g1 --cell 7:3:7 --shapes random,random2,gp,fp,coord --instances 4 --targets 2 --seed $S
run g1-9-3-control.jsonl g1 --cell 9:3:7 --shapes subfield,random --instances 4 --targets 2 --seed $S
run g3-5-2.jsonl g3 --cell 5:2:7 --targets 6 --seed $S
run g3-7-3.jsonl g3 --cell 7:3:7 --targets 6 --seed $S
run g2-5-2-7.jsonl g2 --cell 5:2:7 --shapes random,fp --instances 3 --targets 2 --twists 1,2,l --seed $S
run g2-7-3-7.jsonl g2 --cell 7:3:7 --shapes random,fp --instances 3 --targets 2 --twists 1,2,l --seed $S
for c in 7:3:8 9:4:8 11:5:7 13:6:6; do run g5-${c//:/-}.jsonl g5 --cell $c --instances 4 --targets 4 --betas 3 --seed $S; done
for f in $X/*.jsonl; do b=$(basename $f); if diff <(strip $f) <(strip $R/$b) >/dev/null; then echo "$b identical ($(grep -c outcome $f) measured lines)"; else echo "$b DIFFERS"; fi; done
# (11,4): the gp arm alone (its seeds are per shape, so its lines equal the committed gp lines)
run g1-11-4-gp.jsonl g1 --cell 11:4:6 --shapes gp --instances 3 --targets 2 --seed $S
if diff <(strip $X/g1-11-4-gp.jsonl) <(grep '"shape":"gp"' $R/g1-11-4.jsonl | sed -E 's/,"secs":[0-9.]+//') >/dev/null; then echo "g1-11-4 gp arm identical ($(grep -c outcome $X/g1-11-4-gp.jsonl) measured lines)"; else echo "g1-11-4 gp arm DIFFERS"; fi
