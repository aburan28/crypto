#!/bin/bash
# A18 exploration: the factor base widened at 2^47.2, alone and with A17's
# huge-page filter.  `M1`'s two rows at icv1-f2m61-t158598901-ab42b6c5, the
# suite's recipe (320 signed-Frobenius orbits, 39,040 points) against 400
# and 470 orbits, every other parameter the suite's (params/).  Two
# binaries:
#   v8  = ic-r12-v7-23138d20 (R12's candidate, now the baseline)
#   a17 = ic-a17-v8-fb4284d8 (v8 + the huge-page filter)
# Six arms, four rounds, the order rotating, each process isolated on CPU 2.
# A run the busy gate refuses before it starts waits and is asked for again.
# Usage: explore.sh <scratch dir holding bin/>.
set -u
SP=${1:?scratch dir}
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)
X=$SP/a18
SUITE=$ROOT/research/ic_tool_program/suite/v1/SUITE.json
declare -A BIN=([v8]=$SP/bin/ic-r12-v7-23138d20 [a17]=$SP/bin/ic-a17-v8-fb4284d8)
ARMS=(v8-c320 v8-c400 v8-c470 a17-c320 a17-c400 a17-c470)
seed_of() { jq -r --arg id "k0n61-$1" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3 4; do
  for t in M1-T01 M1-T02; do
    seed=$(seed_of $t)
    k=$(( (round - 1) % 6 ))
    for j in 0 1 2 3 4 5; do
      arm=${ARMS[$(( (k + j) % 6 ))]}
      b=${arm%%-*}; c=${arm#*-}
      d=$X/runs/$arm/k0n61/$t; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      for try in 1 2 3 4 5 6 7 8 9 10; do
        RAYON_NUM_THREADS=1 $SP/bin/isolated_bench run --wait --cpus 2 --out $d/r$round.isolation.jsonl \
          --label "$arm/$t/r$round" -- \
          ${BIN[$b]} price --params $HERE/params/icv1-f2m61-t158598901-ab42b6c5-$t-$c.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
          > /dev/null && break
        echo "REFUSED $arm $t r$round try $try"; sleep 15
      done
      echo "$(date -u +%T) r$round $arm $t cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== A18 exploration done $(date -u)"
