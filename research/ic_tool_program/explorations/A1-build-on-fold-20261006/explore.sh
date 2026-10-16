#!/bin/bash
# The folded build's partition streams and per-run filter on the
# folding-kernel candidate, two arms, M1's two rows at the three largest
# sizes, four rounds, order alternating, each process isolated on CPU 2.
#   sub   = r09-sub a4c19ed7 (v3 + key + the folding kernel)
#   build = sub + A1's streams (54ccf235) + the per-run filter (r10-build 4c71a00e)
set -u
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
X=$SP/buildx
P=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/params/S
SUITE=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/SUITE.json
declare -A BIN=([sub]=$SP/bin/ic-r09sub-on-34ad98f9-a4c19ed7 [build]=$SP/bin/ic-r10build-on-a4c19ed7-4c71a00e)
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
ORDERS=("sub build" "build sub" "sub build" "build sub")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3 4; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    for arm in ${ORDERS[$((round - 1))]}; do
      d=$X/runs/$arm/$1/$2; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      RAYON_NUM_THREADS=1 $SP/isorun-r07.sh $d/r$round "buildx/$arm/$1/$2/r$round" -- \
        ${BIN[$arm]} price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
        || echo "FAILED $arm $1 $2 r$round"
      echo "$(date -u +%T) r$round $arm $1/$2 cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) build_ms=$(jq '[.repetitions[].setup_phases_ns.build/1e6] | sort | .[1]' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== build exploration done $(date -u)"
