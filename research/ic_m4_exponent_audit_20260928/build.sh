#!/usr/bin/env bash
# Build the registered engine for the m = 4 exponent audit (PREREGISTRATION.md §2).
#
# The engine is the source tree of commit 2809b498 -- the commit the frozen
# research/chain_split_order_20260924 totals were built from, and the only tree on
# which they reproduce (§5.1) -- plus one new file, examples/m4_exponent_audit.rs
# from this branch.  Nothing under src/ is changed.  The lock file is the pinned
# copy next to this script (Cargo.lock is git-ignored in the repository).
#
#   WORK=<scratch dir> research/ic_m4_exponent_audit_20260928/build.sh
#
# Produces $WORK/bin/m4_exponent_audit-2809b498 and $WORK/bin/groebner_stage_bench-2809b498.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(git -C "$here" rev-parse --show-toplevel)"
work="${WORK:?set WORK to a scratch directory}"
commit=2809b498f3e45bc0379c0b832c4672c63cae5a73
src="$work/src2809"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$work/target}"

mkdir -p "$src" "$work/bin"
if [ ! -e "$src/Cargo.toml" ]; then
  # src/ include_str!s docs/ic/calibration.json; nothing else outside src/examples/benches.
  git -C "$root" archive "$commit" Cargo.toml src examples benches docs/ic/calibration.json \
    | tar -x -C "$src"
fi
cp "$here/Cargo.lock.pinned" "$src/Cargo.lock"
cp "$root/examples/m4_exponent_audit.rs" "$src/examples/m4_exponent_audit.rs"

(cd "$src" && cargo build --release --locked --example m4_exponent_audit --example groebner_stage_bench)
cp "$CARGO_TARGET_DIR/release/examples/m4_exponent_audit" "$work/bin/m4_exponent_audit-2809b498"
cp "$CARGO_TARGET_DIR/release/examples/groebner_stage_bench" "$work/bin/groebner_stage_bench-2809b498"
sha256sum "$src/examples/m4_exponent_audit.rs" "$src/Cargo.lock" \
  "$work/bin/m4_exponent_audit-2809b498" "$work/bin/groebner_stage_bench-2809b498"
