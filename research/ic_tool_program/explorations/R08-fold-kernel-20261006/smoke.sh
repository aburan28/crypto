#!/bin/bash
# Smoke diagnostic: one isolated process per arm, M1-T01 at n = 53 and 61.
set -u
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
X=$SP/subx/smoke
P=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/params/S
SUITE=$SP/r02b-run-wt/research/ic_tool_program/suite/v1/SUITE.json
BASE=$SP/bin/ic-r06key-on-995ea207-34ad98f9
SUB=$(cat $SP/subx/SUB_BIN)
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for row in "k0n53 M1-T01" "k0n61 M1-T01" "k1n59 M1-T01"; do
  set -- $row; seed=$(seed_of $1 $2)
  for arm in base sub fold shift; do
    case $arm in base) B=$BASE; E=();; sub) B=$SUB; E=();; fold) B=$SUB; E=(KIC_SCAN_SOA=0);; shift) B=$SUB; E=(KIC_SCAN_FOLD=shift);; esac
    d=$X/$arm/$1/$2; mkdir -p $d
    env "${E[@]}" RAYON_NUM_THREADS=1 $SP/isorun-r07.sh $d/r1 "smoke/$arm/$1/$2" -- $B price --params $P/$1/$2.json --json --out $d/r1.price.json --single-target --rho-seed $seed || echo "FAILED $arm $1"
    echo "$arm $1 cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r1.price.json) collect_ms=$(jq '[.repetitions[].setup_phases_ns.collect/1e6] | sort | .[1]' $d/r1.price.json) build_ms=$(jq '[.repetitions[].setup_phases_ns.build/1e6] | sort | .[1]' $d/r1.price.json) scalar=$(jq -r '.certificates.ic.scalar' $d/r1.price.json) contended=$(tail -1 $d/r1.isolation.jsonl | jq .run.contended)"
  done
done
