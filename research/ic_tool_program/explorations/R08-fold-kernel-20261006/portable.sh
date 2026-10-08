#!/bin/bash
# The portable path: R06's candidate (the folding kernel's base) against the
# folding-kernel candidate, both with KIC_SCAN_SIMD=0, so both run the
# scalar scan a CPU without AVX-512F and VPCLMULQDQ runs.  M1's two rows at
# the three largest sizes, three rounds, the order alternating, each process
# isolated on CPU 2.
set -u
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
X=$SP/portx
P=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/params/S
SUITE=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/SUITE.json
declare -A BIN=([base]=$SP/bin/ic-r06key-on-995ea207-34ad98f9 [sub]=$(cat $SP/subx/SUB_BIN))
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
ORDERS=("base sub" "sub base" "base sub")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    for arm in ${ORDERS[$((round - 1))]}; do
      d=$X/runs/$arm/$1/$2; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      KIC_SCAN_SIMD=0 RAYON_NUM_THREADS=1 $SP/isorun-r07.sh $d/r$round "port/$arm/$1/$2/r$round" -- \
        ${BIN[$arm]} price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
        || echo "FAILED $arm $1 $2 r$round"
      echo "$(date -u +%T) r$round $arm $1/$2 cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) collect_ms=$(jq '[.repetitions[].setup_phases_ns.collect/1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== portable exploration done $(date -u)"
