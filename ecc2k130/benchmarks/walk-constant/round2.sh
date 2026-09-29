#!/bin/bash
# WALK-CONSTANT.md section 11: the command lines behind round 2's frozen files.
# Build first: cargo build --release --example ecc2k130_walk_constant (repo
# root) and make -C ecc2k130 walk-constant (the device harness, rule v2).
# About two hours on four cores; each block rewrites its own files.
set -euo pipefail
cd "$(dirname "$0")"
E=../../../target/release/examples/ecc2k130_walk_constant
D=../../build/walk-constant-host-h8

# Merge parting under both rules (emulation) and on the device's own walk.
{
  for n in 23 41 59; do for d in native ecc2k130; do for r in v1 v2; do
    "$E" --n $n --walk table --dist $d --branches 8 --merge 400000 --seed $n --threads 4 --rule $r
  done; done; done
  for d in native ecc2k130; do "$D" --n 23 --walk table --dist $d --merge 40000 --seed 23 --threads 4; done
  "$D" --n 41 --walk table --dist ecc2k130 --merge 4000 --seed 41 --threads 4
} > merge.jsonl 2> merge.log

# The table walk under rule v2: matrix-v2's table rows, same seeds.
{
  for w in "ecc2k130 8" "native 8" "uniform 8" "ecc2k130 16" "native 16"; do
    set -- $w
    "$E" --n 23 --walk table --dist $1 --branches $2 --walks 8 --trials 200000 --seed 23 --threads 4 --rule v2
  done
  for d in ecc2k130 native; do
    "$E" --n 37 --walk table --dist $d --branches 8 --walks 16 --trials 40000 --seed 37 --threads 4 --rule v2
  done
  "$E" --n 59 --walk table --dist ecc2k130 --branches 8 --walks 16 --trials 20000 --seed 59 --threads 4 --rule v2
  for d in ecc2k130 native; do
    "$E" --n 41 --walk table --dist $d --branches 8 --walks 16 --trials 8000 --seed 41 --threads 4 --rule v2
  done
} > matrix-v3.jsonl 2> matrix-v3.log

# The device's own walk: device-v2's table rows, same seeds.
{
  "$D" --n 23 --walk table --dist ecc2k130 --walks 8 --trials 20000 --seed 230 --threads 4
  "$D" --n 23 --walk table --dist native --walks 8 --trials 60000 --seed 230 --threads 4
  "$D" --n 41 --walk table --dist ecc2k130 --walks 16 --trials 300 --seed 410 --threads 4
  # The replicate declared in section 11 item 3, after the first two rows
  # (a new seed; appended by hand to the files the lines above wrote).
  "$D" --n 23 --walk table --dist ecc2k130 --walks 8 --trials 40000 --seed 231 --threads 4
} > device-v3.jsonl 2> device-v3.log

# The residual probe: two uniform branches make rule v2's six-step survivors
# common enough to count.
"$E" --n 23 --walk table --dist uniform --branches 2 --walks 8 --trials 200000 --seed 232 --threads 2 --rule v2 \
  > residual-check.jsonl 2> residual-check.log

python3 fruitless_patterns.py --residual 8
