#!/bin/bash
# A2 exploration (huge pages advised on the pair tables) on v3, two arms:
# M1's two rows at the three largest sizes, three rounds, order alternating,
# each process isolated on CPU 2 (isorun-r07.sh).
#   base = main 995ea207
#   cand = 995ea207 plus a2-thp (b2e9d484): the builders advise their own tables
set -u
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
X=$SP/a2x
P=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/params/S
SUITE=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/SUITE.json
declare -A BIN=([base]=$SP/bin/ic-main-995ea207 [cand]=$SP/bin/ic-a2-thp-on-995ea207-b2e9d484)
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
ORDERS=("base cand" "cand base" "base cand")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    for arm in ${ORDERS[$((round - 1))]}; do
      d=$X/runs2/$arm/$1/$2; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      RAYON_NUM_THREADS=1 $SP/isorun-r07.sh $d/r$round "$arm/$1/$2/r$round" -- \
        ${BIN[$arm]} price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
        || echo "FAILED $arm $1 $2 r$round"
      echo "$(date -u +%T) r$round $arm $1/$2 cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== a2 exploration (two arms) done $(date -u)"
