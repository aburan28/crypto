#!/bin/sh
# Amendment 1: the m = 4 grid with the wide (u128) chain solver. Thin
# orchestration only. Writes runs/m4-wide (must not exist).
set -eu
cd "$(dirname "$0")/../.."
cargo build --release --example ic_algebraic_floor
B=target/release/examples/ic_algebraic_floor
S=research/ic_m3_algebraic_vs_mitm_20261014
test ! -e $S/runs/m4-wide
$B run --m 4 --cells 1:11:3,0:13:4,1:17:5,1:19:5,1:23:6,1:29:8,0:31:8 \
  --arms wide-chain \
  --pool 1024 --kyes 6 --kno 6 --budget-secs 120 --mem-gib 4 --stop-after 4 \
  --seed 20261014 --out $S/runs/m4-wide 2> $S/runs-m4-wide.log
