#!/bin/bash
# Scan-stage diagnostic on v3: R04's probes (cherry-picked onto 995ea207 as
# 65cc4f6e), M1's two rows at the three largest sizes, two rounds, isolated.
set -u
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
X=$SP/scanp
P=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/params/S
SUITE=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/SUITE.json
BIN=$SP/bin/ic-r04probes-on-995ea207-65cc4f6e
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    d=$X/runs/$1/$2; mkdir -p $d
    [ -s $d/r$round.price.json ] && continue
    RAYON_NUM_THREADS=1 $SP/isorun-r07.sh $d/r$round "scanp/$1/$2/r$round" -- \
      $BIN price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
      || echo "FAILED $1 $2 r$round"
    echo "$(date -u +%T) r$round $1/$2 scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
  done
done
echo "== scan probes done $(date -u)"
