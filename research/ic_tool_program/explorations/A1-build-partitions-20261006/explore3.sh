#!/bin/bash
# A1 exploration 3 on v3: three arms, M1's two rows at the three largest
# sizes, three rounds, the order rotating by round, isolated on CPU 2.
#   base = main 995ea207 (v3)
#   v1   = a1-build 41d80f7b (twelve-byte streams)
#   v2b  = a1-build 76a8de3e (tagged four-byte streams, plain per-run filters)
set -u
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
X=$SP/a1x
P=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/params/S
SUITE=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/SUITE.json
declare -A BIN=([base]=$SP/bin/ic-main-995ea207 [v1]=$SP/bin/ic-a1-cand-on-995ea207-41d80f7b [v2b]=$SP/bin/ic-a1-cand-on-995ea207-76a8de3e)
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
ORDERS=("base v1 v2b" "v1 v2b base" "v2b base v1")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    for arm in ${ORDERS[$((round - 1))]}; do
      d=$X/runs3/$arm/$1/$2; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      RAYON_NUM_THREADS=1 $SP/isorun-r07.sh $d/r$round "$arm/$1/$2/r$round" -- \
        ${BIN[$arm]} price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
        || echo "FAILED $arm $1 $2 r$round"
      echo "$(date -u +%T) r$round $arm $1/$2 cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) build_ms=$(jq '[.repetitions[].setup_phases_ns.build/1e6] | sort | .[1]' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== exploration 3 (base, v1, v2b) done $(date -u)"
