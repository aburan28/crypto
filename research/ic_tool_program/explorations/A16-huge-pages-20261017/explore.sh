#!/bin/bash
# A16 exploration: the folded pair table's scan arrays collapsed onto 2 MiB
# pages after the build (MADV_HUGEPAGE + MADV_COLLAPSE), on v8.  One binary,
# three arms chosen by KIC_HUGE_PAGES:
#   off    = KIC_HUGE_PAGES=0       (v8's pages)
#   filter = KIC_HUGE_PAGES=filter  (the presence filter alone)
#   all    = unset                  (filter, bucket starts, stored keys)
# M1's two rows at the three target sizes, six rounds, the order rotating,
# each process isolated on CPU 2.  A run the busy gate refuses before it
# starts waits and is asked for again.
# Usage: explore.sh <scratch dir holding bin/> <a16 commit>.
set -u
SP=${1:?scratch dir}; AC=${2:?a16 commit}
ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
X=$SP/a16
P=$ROOT/research/ic_tool_program/suite/v1/params/S
SUITE=$ROOT/research/ic_tool_program/suite/v1/SUITE.json
BIN=$SP/bin/ic-a16-v8-${AC:0:8}
declare -A ENVS=([off]="KIC_HUGE_PAGES=0" [filter]="KIC_HUGE_PAGES=filter" [all]="KIC_HUGE_PAGES_UNSET=1")
ROWS=("k0n53 M1-T01" "k0n53 M1-T02" "k1n59 M1-T01" "k1n59 M1-T02" "k0n61 M1-T01" "k0n61 M1-T02")
ORDERS=("off filter all" "filter all off" "all off filter" "off filter all" "filter all off" "all off filter")
seed_of() { jq -r --arg id "$1-$2" '.rows[] | select(.id == $id) | .rho_seed' $SUITE; }
for round in 1 2 3 4 5 6; do
  for row in "${ROWS[@]}"; do
    set -- $row
    seed=$(seed_of $1 $2)
    for arm in ${ORDERS[$((round - 1))]}; do
      d=$X/runs/$arm/$1/$2; mkdir -p $d
      [ -s $d/r$round.price.json ] && continue
      for try in 1 2 3 4 5 6 7 8 9 10; do
        env -u KIC_HUGE_PAGES ${ENVS[$arm]} RAYON_NUM_THREADS=1 $SP/bin/isolated_bench run --wait --cpus 2 \
          --out $d/r$round.isolation.jsonl --label "$arm/$1/$2/r$round" -- \
          $BIN price --params $P/$1/$2.json --json --out $d/r$round.price.json --single-target --rho-seed $seed \
          > /dev/null && break
        echo "REFUSED $arm $1 $2 r$round try $try"; sleep 15
      done
      echo "$(date -u +%T) r$round $arm $1/$2 cold_ms=$(jq '[.repetitions[] | (.setup_ns + .ic_online.wall_ns) / 1e6] | sort | .[length/2|floor]' $d/r$round.price.json 2>/dev/null) scalar=$(jq -r '.certificates.ic.scalar' $d/r$round.price.json 2>/dev/null)"
    done
  done
done
echo "== A16 exploration done $(date -u)"
