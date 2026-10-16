#!/bin/bash
# The scan's and the build's additions by carry-less folds, on R06's
# candidate (r09-sub), four arms, M1's two rows at the three largest sizes,
# four rounds in a Latin square of orders, each process isolated on CPU 2.
#   base  = r06-key 34ad98f9 (v3 + the key's GFNI and funnel shifts)
#   sub   = base + the folding kernel, two chains, coordinate slices, x₃ to
#           keys, the build's rows past their own orbit (r09-sub a4c19ed7)
#   fold  = the sub binary with KIC_SCAN_SOA=0: the folding kernel on the
#           point list everywhere
#   shift = the sub binary with KIC_SCAN_FOLD=shift: the shift kernel, which
#           takes no slices, so the base's paths (a check that the
#           candidate's fallback costs nothing)
set -u
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
X=$SP/subx
P=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/params/S
SUITE=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/SUITE.json
SUB=${SUB_BIN:?set SUB_BIN to the candidate binary}
declare -A BIN=([base]=$SP/bin/ic-r06key-on-995ea207-34ad98f9 [sub]=$SUB [fold]=$SUB [shift]=$SUB)
declare -A ENVS=([base]="" [sub]="" [fold]="KIC_SCAN_SOA=0" [shift]="KIC_SCAN_FOLD=shift")
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
ORDERS=("base sub fold shift" "sub fold shift base" "fold shift base sub" "shift base sub fold")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3 4; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    for arm in ${ORDERS[$((round - 1))]}; do
      d=$X/runs/$arm/$1/$2; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      env ${ENVS[$arm]} RAYON_NUM_THREADS=1 $SP/isorun-r07.sh $d/r$round "$arm/$1/$2/r$round" -- \
        ${BIN[$arm]} price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
        || echo "FAILED $arm $1 $2 r$round"
      echo "$(date -u +%T) r$round $arm $1/$2 cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) collect_ms=$(jq '[.repetitions[].setup_phases_ns.collect/1e6] | sort | .[1]' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== sub exploration done $(date -u)"
