#!/bin/bash
# Arm B5a on main: build its ic with its commit, run B5a's tests and the
# conformance suite through B5a (the arm's root holds v2-b5a), natively.
SP=$SP
taskset -cp 0-3 $$ >/dev/null
LOG=$SP/tbarm-b5a-check.log
B=$SP/bin/isolated_bench-r07
cd $SP/tbarm-wt
H=$(git rev-parse tbarm-B5a)
BIN=$SP/bin/ic-tbarm-B5a-${H:0:9}
{
  echo "== B5a $H start $(date -u)"
  git checkout -q --detach $H && cp $SP/main-wt/Cargo.lock . || { echo "checkout failed"; exit 1; }
  IC_BUILD_COMMIT=$H CARGO_TARGET_DIR=$SP/tbarm-target $B busy --wait -- cargo build --locked --release --bin ic 2>&1 | grep -E '^(error|warning)' -A8 | head -60
  cp $SP/tbarm-target/release/ic $BIN && echo "built $(sha256sum $BIN | cut -c1-16)"
  CARGO_TARGET_DIR=$SP/tbarm-target $B busy --wait -- bash -c 'cargo test --locked --release --lib -- fpk_curve curve_id gaudry_cubic rho_bignum koblitz_two_word koblitz_multi koblitz_index_calculus koblitz_fast ic_boundary 2>&1 | grep -E "test result|FAILED|panicked|^error"; cargo test --locked --release --bin ic 2>&1 | grep -E "test result|FAILED|panicked|^error"; cargo test --locked --release --test curve_id 2>&1 | grep -E "test result|FAILED|panicked|^error"'
  $B busy --wait -- $SP/bin/icprog-b5a-dev conformance --ic $BIN --steps B0,B1,B3,B2,B2b,B7a,B3b,B4,B5a --root $SP/b5an-wt --build-commit $H --out $SP/tbarm-conf-B5a.json > /dev/null
  echo "conformance exit $? $(jq -c '{passed, cases: (.results|length), failed: [.results[]|select(.pass|not)|.id]}' $SP/tbarm-conf-B5a.json)"
  echo "== B5a done $(date -u)"
} >> $LOG 2>&1
