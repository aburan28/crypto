#!/bin/bash
# Three streams so peak memory stays well under the host's 15 GiB:
#  A: IC L=16384 K=880 (~5 GiB) then IC L=4096 K=660 (~2.5 GiB)
#  B: IC n=61 K=600 (~2.5 GiB)
#  C: rho R3 in every cell, then IC n=37 and n=41 (small)
cd "$(dirname "$0")"
( ./sweep_cell.py callgrind_ic 53 16384 880 n53-strong-sweep-L16384-v1; \
  ./sweep_cell.py callgrind_ic 53 4096 660 n53-strong-sweep-L4096-v1; echo STREAM_A_DONE ) > stream_A.out 2>&1 &
( ./sweep_cell.py callgrind_ic 61 1024 600 n61-strong-sweep-L1024-v1; echo STREAM_B_DONE ) > stream_B.out 2>&1 &
( while read -r N L K CORPUS; do ./sweep_cell.py callgrind_rho "$N" "$L" "$K" "$CORPUS"; done < cells.txt; \
  ./sweep_cell.py callgrind_ic 37 1024 42 n37-strong-sweep-L1024-v1; \
  ./sweep_cell.py callgrind_ic 41 1024 255 n41-strong-sweep-L1024-v1; echo STREAM_C_DONE ) > stream_C.out 2>&1 &
wait
echo ALL_STREAMS_DONE > streams_done.txt
