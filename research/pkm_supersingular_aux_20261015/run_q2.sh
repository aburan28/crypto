#!/bin/bash
# Q2 of PROTOCOL.md: F4 solving degree, one run at a time, single-threaded,
# niced, 4 GB address space. No run starts after 180 minutes.
#   BIN=path/to/examples/pkm_ss_tower research/pkm_supersingular_aux_20261015/run_q2.sh
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:?set BIN to the pkm_ss_tower binary}
cd "$here/q2/runs"
export RAYON_NUM_THREADS=1
P0=786433; T24=14680063; T31=2046820351; S25=18874367
common=(--targets 2 --budget 3600 --max-nnz 2000000000)
start=$(date +%s); deadline=$((start + 180 * 60))
run() {
  name=$1; shift
  if [ "$(date +%s)" -gt "$deadline" ]; then
    echo "$(date -u +%FT%TZ) not run (deadline) $name: $*" >> ../../progress.txt; return
  fi
  echo "$(date -u +%FT%TZ) start $name: $*" >> ../../progress.txt
  nice -n 10 prlimit --as=4294967296 -- "$bin" f4 "${common[@]}" "$@" --out "$name.jsonl" 2> "$name.log"
  echo "$(date -u +%FT%TZ) end $name exit $?" >> ../../progress.txt
}
M2=(--m 2 --t-min 1 --t-max 9)
run P0-ordinary --p $P0 --arm ordinary --m 2 --t-min 1 --t-max 5
run P0-split --p $P0 --arm split --m 2 --t-min 1 --t-max 5
for pr in T24 T31 S25; do
  p=${!pr}
  arms="ss-j1728 ss-j1728-twist ordinary torus"
  [ "$pr" = S25 ] && arms="ss-j1728 ss-j1728-twist ss-j0 ss-j0-twist ordinary torus"
  for arm in $arms; do run "$pr-$arm-m2" --p "$p" --arm "$arm" "${M2[@]}"; done
done
CH=(--m 3,4 --t-min 2 --t-max-m3 4 --t-max-m4 3)
for arm in ss-j1728 ordinary torus; do run "T31-$arm-m34" --p $T31 --arm $arm "${CH[@]}"; done
run T31-ss-j1728-null-m2 --p $T31 --arm ss-j1728 --control null "${M2[@]}"
run T31-ss-j1728-null-m34 --p $T31 --arm ss-j1728 --control null "${CH[@]}"
for arm in ss-j1728 ordinary torus; do run "T31-$arm-m2-N20" --p $T31 --arm $arm --m 2 --t-min 10 --t-max 10; done
echo "$(date -u +%FT%TZ) q2 all done" >> ../../progress.txt
