#!/bin/bash
# A1 exploration: the folded build's stage timings on v3 (995ea207 + probes),
# M1's two rows at the three largest sizes, two rounds, isolated on CPU 2.
set -u
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
X=$SP/a1x
P=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/params/S
SUITE=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/SUITE.json
BIN=$SP/bin/ic-a1-probes-on-995ea207
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    d=$X/runs/probe/$1/$2; mkdir -p $d
    [ -s $d/r$round.price.json ] && continue
    KIC_BUILD_PROBES=1 RAYON_NUM_THREADS=1 $SP/isorun-r07.sh $d/r$round "probe/$1/$2/r$round" -- \
      $BIN price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
      || echo "FAILED $1 $2 r$round"
    echo "$(date -u +%T) r$round $1/$2 build_ms=$(jq '.repetitions[0].setup_phases_ns.build/1e6' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
  done
done
echo "== probes done $(date -u)"
