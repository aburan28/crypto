#!/bin/sh
# Run one f4_lane_probe mode over the four production probe cells.
# BIN, MODE and EXTRA may be set from the environment.
cd /Volumes/SSD990-2/crypto/targets/lit2026-w5/scratch || exit 1
BIN=${BIN:-/Volumes/SSD990-2/crypto/targets/lit2026-w5/release/examples/f4_lane_probe}
MODE=${MODE:-profile}
run() {
  echo "== $1"
  shift
  $BIN $MODE "$@" $EXTRA
}
run K1_17 --a 1 --n 17 --dim 6 --targets 8
run K1_23 --a 1 --n 23 --dim 7 --targets 4
run K0_23 --a 0 --n 23 --dim 8 --targets 4
run K0_31 --a 0 --n 31 --dim 9 --targets 1 --budget 2000
