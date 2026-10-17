#!/bin/bash
# A19 exploration: A17's huge-page filter plus the folded build's
# 64 partitions (not 256), on v8.
# Two arms, M1's two rows at the three target sizes, six rounds, the
# order alternating, each process isolated on CPU 2.
#   v8   = ic-r12-v7-23138d20 (R12's candidate, now the baseline)
#   a19  = v8 + the huge-page filter + 64 partitions (6c998fb0)
# Usage: explore.sh <scratch dir holding bin/> <a19 commit>.
set -u
SP=${1:?scratch dir}; AC=${2:?a19 commit}
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
X=$SP/a19
P=$ROOT/research/ic_tool_program/suite/v1/params/S
SUITE=$ROOT/research/ic_tool_program/suite/v1/SUITE.json
declare -A BIN=([v8]=$SP/bin/ic-r12-v7-23138d20 [a19]=$SP/bin/ic-a19-v8-${AC:0:8})
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
ORDERS=("v8 a19" "a19 v8" "v8 a19" "a19 v8" "v8 a19" "a19 v8")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3 4 5 6; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    for arm in ${ORDERS[$((round - 1))]}; do
      d=$X/runs/$arm/$1/$2; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      # The gate refuses a run before it starts while the last run's
      # pressure decays (PSI avg10 above 5): wait and ask again.
      for try in 1 2 3 4 5 6 7 8 9 10; do
        RAYON_NUM_THREADS=1 $SP/bin/isolated_bench run --wait --cpus 2 --out $d/r$round.isolation.jsonl \
          --label "$arm/$1/$2/r$round" -- \
          ${BIN[$arm]} price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
          > /dev/null && break
        echo "REFUSED $arm $1 $2 r$round try $try"; sleep 15
      done
      echo "$(date -u +%T) r$round $arm $1/$2 cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== A19 exploration done $(date -u)"
