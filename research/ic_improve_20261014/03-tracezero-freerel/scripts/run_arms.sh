#!/bin/bash
# Usage: run_arms.sh <n> <K> <seed_first> <seed_last>
# Runs the five arms (inc, A, F, R, FR) per seed, interleaved, one process each.
set -u
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
D=research/ic_improve_20261014/03-tracezero-freerel
DATA=$D/data
BIN=${IC03_BIN:?set IC03_BIN to a directory holding frozen copies of both binaries (sha256 in data/binaries.sha256)}
CAP=3221225472
n=$1; K=$2
head -8 research/ic_next10_20261014/09-equal-advice/data/targets_n$n.txt > $DATA/targets_n$n.txt
for seed in $(seq $3 $4); do
  for arm in inc A F R FR; do
    out=$DATA/n${n}_K${K}_s${seed}_${arm}
    case $arm in
      inc) env_args=""; exe=$BIN/koblitz_orbit_dlp_fast_online ;;
      A)   env_args="IC03_TRACE_FILTER=0 IC03_FREE_RELATIONS=0"; exe=$BIN/ic_improve_03_tracezero_freerel ;;
      F)   env_args="IC03_TRACE_FILTER=1 IC03_FREE_RELATIONS=0"; exe=$BIN/ic_improve_03_tracezero_freerel ;;
      R)   env_args="IC03_TRACE_FILTER=0 IC03_FREE_RELATIONS=1"; exe=$BIN/ic_improve_03_tracezero_freerel ;;
      FR)  env_args="IC03_TRACE_FILTER=1 IC03_FREE_RELATIONS=1"; exe=$BIN/ic_improve_03_tracezero_freerel ;;
    esac
    echo "$(date -u +%FT%TZ) start $out" >> $DATA/jobs.log
    env $env_args timeout 3000 nice prlimit --as=$CAP -- $exe construct:$n:0:$K $DATA/targets_n$n.txt $seed $out.jsonl > $out.summary.json 2> $out.stderr
    echo "$(date -u +%FT%TZ) end $out exit=$?" >> $DATA/jobs.log
  done
done
