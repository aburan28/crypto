#!/bin/bash
# Usage: run_job.sh <arm> <n> <K> <seed>
#   arm = inc   frozen incumbent binary, default environment
#         free  candidate binary, KIC_FREE_RELATIONS=1
#         aa    candidate binary, KIC_FREE_RELATIONS unset, with KIC_DUMP_RANK (A/A)
#         incd  frozen incumbent binary, with KIC_DUMP_RANK (A/A partner of aa)
# Binaries: $IC2_02_BIN/{incumbent,candidate} (hashes in data/binaries.sha256).
# Rank traces (A/A only) go to $IC2_02_TRACES (scratch; too large for the repo).
set -u
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
D=research/ic_improve_20261015/02-freerel-production
DATA=$D/data/runs
BIN=${IC2_02_BIN:?set IC2_02_BIN}
TR=${IC2_02_TRACES:-/tmp}
CAP=3221225472
arm=$1; n=$2; K=$3; seed=$4
mkdir -p $DATA
T=$D/data/targets_n$n.txt
[ -s $T ] || head -8 research/ic_next10_20261014/09-equal-advice/data/targets_n$n.txt > $T
out=$DATA/n${n}_K${K}_s${seed}_${arm}
unset KIC_FREE_RELATIONS KIC_DUMP_RANK KIC_PAIR_SYMMETRY KIC_ROOT_BLOOM KIC_RANK_PROBE_CAP
case $arm in
  inc)  exe=$BIN/incumbent; envs=() ;;
  free) exe=$BIN/candidate; envs=(KIC_FREE_RELATIONS=1) ;;
  aa)   exe=$BIN/candidate; envs=(KIC_DUMP_RANK=$TR/n${n}_K${K}_s${seed}_aa.rank.jsonl) ;;
  incd) exe=$BIN/incumbent; envs=(KIC_DUMP_RANK=$TR/n${n}_K${K}_s${seed}_incd.rank.jsonl) ;;
  *) echo "bad arm $arm"; exit 2 ;;
esac
echo "$(date -u +%FT%TZ) start $out" >> $D/data/jobs.log
env "${envs[@]}" timeout 3000 nice prlimit --as=$CAP -- $exe construct:$n:0:$K $T $seed $out.jsonl > $out.summary.json 2> $out.stderr
rc=$?
echo "$(date -u +%FT%TZ) end $out exit=$rc" >> $D/data/jobs.log
[ -s $out.stderr ] || rm -f $out.stderr
