#!/bin/sh
# One pass of every candidate arm on one cell (default K1/2^17 dim 6), totals only.
cd /Volumes/SSD990-2/crypto/targets/lit2026-w5/scratch || exit 1
B=/Volumes/SSD990-2/crypto/targets/lit2026-w5/bin/f6_ic_probe.cand
FLAGS=${FLAGS:-"--a 1 --n 17 --dim 6 --targets 8"}
LANES=${LANES:-8}
echo prev;  $B $FLAGS --engines f4,f6p | grep total
echo gate;  KIC_F6_GATE_FIRST=1 $B $FLAGS --engines f4,f6p | grep total
echo lanes; KIC_F4_LANES=$LANES $B $FLAGS --engines f4,f6p | grep total
echo rref;  KIC_F4_INHERIT_ROOT=rref $B $FLAGS --engines f4,f6p | grep total
