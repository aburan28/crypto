#!/bin/bash
# A17 exploration: the presence filter copied into a fresh 2 MiB-aligned
# anonymous mapping advised MADV_HUGEPAGE before it is written, on v8.
# Three arms:
#   v8   = ic-r12-v7-23138d20 (R12's candidate, now the baseline)
#   a17  = v8 + the huge-page filter (fb4284d8)
#   heap = the a17 binary with KIC_FILTER_HUGE=0 (the copy, on the heap)
# M1's two rows at the three target sizes, six rounds, the order rotating,
# each process isolated on CPU 2.  A run the busy gate refuses before it
# starts waits and is asked for again.
# Usage: explore.sh <scratch dir holding bin/> <a17 commit>.
set -u
SP=${1:?scratch dir}; AC=${2:?a17 commit}
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
X=$SP/a17
P=$ROOT/research/ic_tool_program/suite/v1/params/S
SUITE=$ROOT/research/ic_tool_program/suite/v1/SUITE.json
declare -A BIN=([v8]=$SP/bin/ic-r12-v7-23138d20 [a17]=$SP/bin/ic-a17-v8-${AC:0:8} [heap]=$SP/bin/ic-a17-v8-${AC:0:8})
declare -A ENVS=([v8]="KIC_ARM=v8" [a17]="KIC_ARM=a17" [heap]="KIC_FILTER_HUGE=0")
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
ORDERS=("v8 a17 heap" "a17 heap v8" "heap v8 a17" "v8 a17 heap" "a17 heap v8" "heap v8 a17")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3 4 5 6; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    for arm in ${ORDERS[$((round - 1))]}; do
      d=$X/runs/$arm/$1/$2; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      for try in 1 2 3 4 5 6 7 8 9 10; do
        env -u KIC_FILTER_HUGE ${ENVS[$arm]} RAYON_NUM_THREADS=1 $SP/bin/isolated_bench run --wait --cpus 2 \
          --out $d/r$round.isolation.jsonl --label "$arm/$1/$2/r$round" -- \
          ${BIN[$arm]} price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
          > /dev/null && break
        echo "REFUSED $arm $1 $2 r$round try $try"; sleep 15
      done
      echo "$(date -u +%T) r$round $arm $1/$2 cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== A17 exploration done $(date -u)"
