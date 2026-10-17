#!/bin/sh
# The frozen claim grid (PROTOCOL.md §2-§5). Thin orchestration only: the
# native instrument does the work. Writes runs/m3 and runs/m4 (must not exist).
set -eu
cd "$(dirname "$0")/../.."
cargo build --release --example ic_algebraic_floor
B=target/release/examples/ic_algebraic_floor
S=research/ic_m3_algebraic_vs_mitm_20261014
test ! -e $S/runs/m3 && test ! -e $S/runs/m4
$B run --m 3 --cells 0:13:5,1:17:6,1:19:7,1:23:8,1:29:10,0:31:11 \
  --arms f4-chain,f6-ic,sat,s4-hybrid-f10,s4-hybrid-f12 \
  --pool 2048 --kyes 6 --kno 6 --budget-secs 120 --mem-gib 4 --stop-after 4 \
  --seed 20261014 --out $S/runs/m3 2> $S/runs-m3.log &
$B run --m 4 --cells 1:11:3,0:13:4,1:17:5,1:19:5,1:23:6,1:29:8,0:31:8 \
  --arms f4-chain,f6-ic,sat,s22-hybrid-g8 \
  --pool 1024 --kyes 6 --kno 6 --budget-secs 120 --mem-gib 4 --stop-after 4 \
  --seed 20261014 --out $S/runs/m4 2> $S/runs-m4.log &
wait
