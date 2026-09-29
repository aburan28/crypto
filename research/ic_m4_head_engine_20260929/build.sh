#!/usr/bin/env bash
# Build the head engine for the m = 4 head-engine rerun (PREREGISTRATION.md §2).
#
# The engine is the source tree of commit 4ff512f2 (origin/main when this rerun was
# registered), which already contains examples/m4_exponent_audit.rs.  Nothing is
# changed or copied in.  The lock file is the first audit's pinned copy (Cargo.lock
# is git-ignored in the repository); it built this tree with --locked unchanged.
# The build runs under the benchmark lock (tools/isolated_bench.py busy), so it can
# never overlap a timed run.
#
#   WORK=<scratch dir> research/ic_m4_head_engine_20260929/build.sh
#
# Produces $WORK/bin/m4_exponent_audit-4ff512f2 and $WORK/bin/groebner_stage_bench-4ff512f2.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(git -C "$here" rev-parse --show-toplevel)"
work="${WORK:?set WORK to a scratch directory}"
commit=4ff512f25813cb66c896860890eb405704f8bf00
src="$work/src4ff512f2"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$work/target}"

mkdir -p "$src" "$work/bin"
if [ ! -e "$src/Cargo.toml" ]; then
  # src/ include_str!s docs/ic/calibration.json; nothing else outside src/examples/benches.
  git -C "$root" archive "$commit" Cargo.toml src examples benches docs/ic/calibration.json \
    | tar -x -C "$src"
fi
cp "$here/../ic_m4_exponent_audit_20260928/Cargo.lock.pinned" "$src/Cargo.lock"
# examples/m4_exponent_audit.rs is part of this commit (merged in #925); nothing is copied in.

(cd "$src" && python3 "$root/tools/isolated_bench.py" busy -- cargo build --release --locked --example m4_exponent_audit --example groebner_stage_bench)
cp "$CARGO_TARGET_DIR/release/examples/m4_exponent_audit" "$work/bin/m4_exponent_audit-4ff512f2"
cp "$CARGO_TARGET_DIR/release/examples/groebner_stage_bench" "$work/bin/groebner_stage_bench-4ff512f2"
sha256sum "$src/examples/m4_exponent_audit.rs" "$src/Cargo.lock" \
  "$work/bin/m4_exponent_audit-4ff512f2" "$work/bin/groebner_stage_bench-4ff512f2"
