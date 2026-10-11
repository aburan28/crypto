#!/bin/sh
S=/tmp/claude-0/-home-user/5db8738c-4b29-5210-9888-79d70024753c/scratchpad
REPO=/home/user/crypto
RID=run-$(date -u +%Y%m%dT%H%MZ)-verify
W=$S/$RID; mkdir -p $W/l1probe $W/toyprobe $W/failrate
echo $RID > $S/verify_rid.txt
git -C $REPO rev-parse HEAD > $W/commit.txt
git -C $REPO status --short -- eigen_seta > $W/worktree_status.txt
cd $REPO/eigen_seta
CARGO_TARGET_DIR=$S/tgt_verify CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release --example ct_probe --example dudect_decaps --example run > $W/build.log 2>&1
B=$S/tgt_verify/release/examples
cp $B/ct_probe $B/dudect_decaps $B/run $W/
date -u +%FT%TZ > $W/start.txt
probe() { # dir set key ct mode
  ( cd $W/$1 && valgrind --tool=callgrind --toggle-collect='*probe_core*' --callgrind-out-file=cg.$3.$4.$5 ../ct_probe $3 $4 $5 $2 > out.$3.$4.$5 2> err.$3.$4.$5 ) ; }
# phase 1: L1 probes (4 parallel) + native failure-rate stream
probe l1probe L1 1 1 valid &
probe l1probe L1 1 2 valid &
probe l1probe L1 2 1 valid &
probe l1probe L1 1 1 foreign &
( cd $W/failrate && ../run Toy16 3000 > toy16.txt 2>&1; ../run Toy32 600 > toy32.txt 2>&1; ../run L1 20 > l1.txt 2>&1 ) &
wait
date -u +%FT%TZ > $W/phase1_end.txt
# phase 2: Toy16 grid, 3 streams
( probe toyprobe Toy16 1 1 valid; probe toyprobe Toy16 1 2 valid; probe toyprobe Toy16 1 1 foreign; probe toyprobe Toy16 4 1 valid ) &
( probe toyprobe Toy16 2 1 valid; probe toyprobe Toy16 2 2 valid; probe toyprobe Toy16 2 1 foreign; probe toyprobe Toy16 4 2 valid ) &
( probe toyprobe Toy16 3 1 valid; probe toyprobe Toy16 3 2 valid; probe toyprobe Toy16 3 1 foreign; probe toyprobe Toy16 4 1 foreign ) &
wait
$S/ctsummary.sh $W/l1probe > $W/l1probe/ctprobe_summary.txt
$S/ctsummary.sh $W/toyprobe > $W/toyprobe/ctprobe_summary.txt
# phase 3: timing test alone
( cd $W && ./dudect_decaps 1000 > dudect_output.txt 2>&1 )
date -u +%FT%TZ > $W/end.txt
echo done > $W/DONE
