#!/bin/sh
# PROTOCOL.md: both arms on three ladders, one binary. From the repository root.
set -e
B=target/release/examples/groebner_stage_bench
R=research/pdp_shape_replay_20261015/raw
for L in chain chain-holdout holdout; do
  $B --ladder $L --label default --out $R/$L-default > $R/$L-default.log 2>&1
  KIC_F4_ROOT_REPLAY=1 $B --ladder $L --label replay --out $R/$L-replay > $R/$L-replay.log 2>&1
done
echo done > $R/DONE
