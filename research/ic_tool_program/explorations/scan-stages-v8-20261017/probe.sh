#!/bin/bash
# Scan-stage diagnostic on v8: R04's probes ported onto 23138d20 as 786c06f7
# (probes-v8.patch), M1's two rows at the three largest sizes, four rounds,
# isolated on CPU 2.  A run the busy gate refuses before it starts waits and
# is asked for again.  Usage: probe.sh <scratch dir holding bin/>.
set -u
SP=${1:?scratch dir}
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
X=$SP/scanp-v8
P=$ROOT/research/ic_tool_program/suite/v1/params/S
SUITE=$ROOT/research/ic_tool_program/suite/v1/SUITE.json
BIN=$SP/bin/ic-r04probes-on-v8-786c06f7
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3 4; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    d=$X/runs/$1/$2; mkdir -p $d
    [ -s $d/r$round.price.json ] && continue
    for try in 1 2 3 4 5 6 7 8 9 10; do
      RAYON_NUM_THREADS=1 $SP/bin/isolated_bench run --wait --cpus 2 --out $d/r$round.isolation.jsonl \
        --label "scanp-v8/$1/$2/r$round" -- \
        $BIN price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
        > /dev/null && break
      echo "REFUSED $1 $2 r$round try $try"; sleep 15
    done
    echo "$(date -u +%T) r$round $1/$2 scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
  done
done
echo "== scan probes on v8 done $(date -u)"
