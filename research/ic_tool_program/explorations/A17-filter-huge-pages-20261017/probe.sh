#!/bin/bash
# Stage probes on v8 and on A17, interleaved: R04's probes on each
# (ic-r04probes-on-v8-786c06f7, ic-r04probes-on-a17-a610cdab), M1's two rows
# at the three target sizes, four rounds, the order alternating, isolated
# on CPU 2.  Usage: probe.sh <scratch dir holding bin/>.
set -u
SP=${1:?scratch dir}
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
X=$SP/a17-probes
P=$ROOT/research/ic_tool_program/suite/v1/params/S
SUITE=$ROOT/research/ic_tool_program/suite/v1/SUITE.json
declare -A BIN=([v8]=$SP/bin/ic-r04probes-on-v8-786c06f7 [a17]=$SP/bin/ic-r04probes-on-a17-a610cdab)
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
ORDERS=("v8 a17" "a17 v8" "v8 a17" "a17 v8")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3 4; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    for arm in ${ORDERS[$((round - 1))]}; do
      d=$X/runs/$arm/$1/$2; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      for try in 1 2 3 4 5 6 7 8 9 10; do
        RAYON_NUM_THREADS=1 $SP/bin/isolated_bench run --wait --cpus 2 --out $d/r$round.isolation.jsonl \
          --label "probes/$arm/$1/$2/r$round" -- \
          ${BIN[$arm]} price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
          > /dev/null && break
        echo "REFUSED $arm $1 $2 r$round try $try"; sleep 15
      done
      echo "$(date -u +%T) r$round $arm $1/$2 scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== A17 probes done $(date -u)"
