#!/bin/bash
# Usage: run_stream.sh <stream-name> ; runs one preregistered stream of jobs sequentially.
set -u
cd /home/user/crypto
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
DIR=research/ic_next10_20261014/09-equal-advice/data
EX=target/release/examples
CAP=1073741824
run() { local name=$1; shift; echo "$(date -u +%FT%TZ) start $name" >> $DIR/jobs.log
  timeout 2700 nice prlimit --as=$CAP -- "$@" > $DIR/$name.stdout 2> $DIR/$name.stderr
  echo "$(date -u +%FT%TZ) end $name exit=$?" >> $DIR/jobs.log; }
bytes() { python3 -I research/ic_next10_20261014/09-equal-advice/scripts/ic_bytes.py $1 $2; }
rho() { local n=$1 count=$2; local i=0
  while read -r line && [ $i -lt $count ]; do
    xy=$(echo "$line" | tr -d '[] '); 
    KIC_RHO_TARGET_POINT=$xy KIC_RHO_RUNG=3 KIC_RHO_LANES=32 KIC_RHO_DP_BITS=4 \
      timeout 600 nice prlimit --as=$CAP -- $EX/koblitz_rho_batch_ks_strong_online $n 0 signed_frobenius 1 $((610000+i)) \
      > $DIR/rho_n${n}_t$(printf %02d $i).json 2>> $DIR/rho_n${n}.stderr
    i=$((i+1)); done < $DIR/targets_n${n}.txt; }
case $1 in
 s1)
  # (deviation 1) ic_n41_K100 and bl_n41 ran in the first launch; wait for bl_n41.
  # (deviation 3) null control at 1/16 bytes to save precompute
  BL_ONLINE_MAP_SEED=mismatch BL_PREFIXES=1 run bl_n41_mismatch $EX/ic_next10_09_equal_advice bl 41 0 $(( $(bytes 100 41) / 16 )) $DIR/targets_n41.txt $DIR/bl_n41_mismatch.jsonl
  rho 41 64
  run ic_n53_K200 $EX/koblitz_orbit_dlp_fast_online construct:53:0:200 $DIR/targets_n53.txt 7 $DIR/ic_n53_K200.jsonl
  run bl_n53 $EX/ic_next10_09_equal_advice bl 53 0 $(( $(bytes 200 53) / 16 )) $DIR/targets_n53.txt $DIR/bl_n53.jsonl
  rho 53 32
  ;;
 s2)
  # (deviation 2) ic_n61_K256 ran in the first launch; BL at 1/64 of its bytes.
  run bl_n61_64th $EX/ic_next10_09_equal_advice bl 61 0 $(( $(bytes 256 61) / 64 )) $DIR/targets_n61.txt $DIR/bl_n61_64th.jsonl
  run ic_n61_K360 $EX/koblitz_orbit_dlp_fast_online construct:61:0:360 $DIR/targets_n61.txt 7 $DIR/ic_n61_K360.jsonl
  rho 61 64
  ;;
esac
echo "$(date -u +%FT%TZ) stream $1 done" >> $DIR/jobs.log
