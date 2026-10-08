#!/bin/bash
# Pipelined-scan exploration on the Cascade Lake host (no GFNI, VBMI2 or
# VPCLMULQDQ): v3 against v3 plus the cross-trial pipelined scan.  Two arms,
# M1's two rows at the three largest sizes, four rounds, the order
# alternating, each process isolated on CPU 2.
#   base = main's ic at 995ea207 (v3)
#   pipe = 995ea207 + the pipelined scan
set -u
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
X=$SP/pipex
P=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/params/S
SUITE=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/SUITE.json
PIPE=${PIPE_BIN:?set PIPE_BIN to the candidate binary}
declare -A BIN=([base]=$SP/bin/ic-main-995ea207 [pipe]=$PIPE)
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
ORDERS=("base pipe" "pipe base" "base pipe" "pipe base")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3 4; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    for arm in ${ORDERS[$((round - 1))]}; do
      d=$X/runs/$arm/$1/$2; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      RAYON_NUM_THREADS=1 $SP/isorun-r07.sh $d/r$round "pipe/$arm/$1/$2/r$round" -- \
        ${BIN[$arm]} price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
        || echo "FAILED $arm $1 $2 r$round"
      echo "$(date -u +%T) r$round $arm $1/$2 cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) collect_ms=$(jq '[.repetitions[].setup_phases_ns.collect/1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== pipe exploration done $(date -u)"
