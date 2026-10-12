#!/bin/bash
# Track B's chain arms on main: build each arm's ic with its commit, run
# its tests (B4's tree is 6a787e07's, already tested) and the conformance
# suite at its level, natively.  Everything heavy waits for the lock.
SP=/tmp/claude-0/-home-user-crypto/56047a67-af13-5c18-babc-0bd1175ee451/scratchpad
taskset -cp 0-3 $$ >/dev/null
LOG=$SP/tbarm-check.log
B=$SP/bin/isolated_bench-r07
STEPS=""
cd $SP/tbarm-wt
for step in B0 B1 B3 B2 B2b B7a B3b B4; do
  STEPS=${STEPS:+$STEPS,}$step
  H=$(git rev-parse tbarm-$step)
  BIN=$SP/bin/ic-tbarm-$step-${H:0:9}
  {
    echo "== $step $H start $(date -u)"
    git checkout -q --detach $H && cp $SP/main-wt/Cargo.lock . || { echo "checkout failed"; continue; }
    if [ ! -x $BIN ]; then
      IC_BUILD_COMMIT=$H CARGO_TARGET_DIR=$SP/tbarm-target $B busy --wait -- cargo build --release --bin ic 2>&1 | grep -E '^(error|warning)' -A6 | head -40
      cp $SP/tbarm-target/release/ic $BIN && echo "built $(sha256sum $BIN | cut -c1-16)"
    fi
    if [ $step != B4 ]; then
      CARGO_TARGET_DIR=$SP/tbarm-target $B busy --wait -- bash -c 'cargo test --release --lib -- koblitz_two_word koblitz_multi koblitz_index_calculus koblitz_fast rho_bignum curve_id ic_boundary 2>&1 | grep -E "test result|FAILED|panicked|^error"; cargo test --release --bin ic 2>&1 | grep -E "test result|FAILED|panicked|^error"'
    fi
    $B busy --wait -- $SP/bin/icprog-n4dev conformance --ic $BIN --steps $STEPS --root $SP/r07-run-wt --build-commit $H --out $SP/tbarm-conf-$step.json > /dev/null
    echo "conformance $STEPS exit $? $(jq -c '{passed, cases: (.results|length), failed: [.results[]|select(.pass|not)|.id]}' $SP/tbarm-conf-$step.json)"
    echo "== $step done $(date -u)"
  } >> $LOG 2>&1
done
echo "== all done $(date -u)" >> $LOG
