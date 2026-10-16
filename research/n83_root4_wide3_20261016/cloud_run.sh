#!/usr/bin/env bash
# Execute the frozen eight-arm screen on one dedicated Linux allocation.
set -euo pipefail

: "${N83_SOURCE_COMMIT:?}"
: "${N83_OUTPUT_DIR:?}"
: "${N83_INSTANCE_ID:?}"
: "${N83_AVAILABILITY_ZONE:?}"
: "${N83_ALLOCATION_STARTED_UTC:?}"

test "$(git rev-parse HEAD)" = "$N83_SOURCE_COMMIT"
test "$(sha256sum examples/koblitz_orbit_dlp_fast.rs | cut -d' ' -f1)" = \
  f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc
mkdir -p "$N83_OUTPUT_DIR"
date -u +%FT%TZ > "$N83_OUTPUT_DIR/started_utc.txt"
uname -a > "$N83_OUTPUT_DIR/uname.txt"
rustc --version > "$N83_OUTPUT_DIR/rustc-version.txt"
cargo --version > "$N83_OUTPUT_DIR/cargo-version.txt"
git rev-parse HEAD > "$N83_OUTPUT_DIR/source-commit.txt"

export PYTHONDONTWRITEBYTECODE=1
for prior in \
  n83_k3250_sparse_fresh_20261015 \
  n83_k3350_sparse_fresh_20261015 \
  n83_k_abba_20261015 \
  n83_root4_highk_20261015 \
  n83_best_e2e_20261015 \
  n83_legacy_challenger_20261015 \
  n83_root4_bloom_20261015 \
  n83_root11_highk_20261015 \
  n83_root11_edgek_20261015 \
  n83_final_tournament_20261015 \
  n83_edge_confirmation_20261015 \
  n83_final_tier1_20261015 \
  n83_unresolved_tier1_20261015
do
  python3 "research/$prior/generate.py" > "$N83_OUTPUT_DIR/fixture-$prior.json"
done
python3 research/n83_root4_wide3_20261016/generate.py > \
  "$N83_OUTPUT_DIR/fixture-check.json"
python3 research/n83_root4_wide3_20261016/run.py --preflight > \
  "$N83_OUTPUT_DIR/source-preflight.json"

export CARGO_TARGET_DIR="${N83_CARGO_TARGET_DIR:-/mnt/n83-wide3-target}"
export CARGO_BUILD_JOBS=8
cargo build --release --locked --example koblitz_orbit_dlp_fast > \
  "$N83_OUTPUT_DIR/build.stdout.txt" 2> "$N83_OUTPUT_DIR/build.stderr.txt"
sha256sum "$CARGO_TARGET_DIR/release/examples/koblitz_orbit_dlp_fast" > \
  "$N83_OUTPUT_DIR/binary.sha256"

set +e
python3 research/n83_root4_wide3_20261016/run.py \
  --binary "$CARGO_TARGET_DIR/release/examples/koblitz_orbit_dlp_fast" \
  --output "$N83_OUTPUT_DIR/results" > "$N83_OUTPUT_DIR/run.stdout.txt" \
  2> "$N83_OUTPUT_DIR/run.stderr.txt"
producer_status=$?
if test -f "$N83_OUTPUT_DIR/results/source.json"; then
  python3 research/n83_root4_wide3_20261016/verify.py \
    "$N83_OUTPUT_DIR/results" > "$N83_OUTPUT_DIR/verify.stdout.txt" \
    2> "$N83_OUTPUT_DIR/verify.stderr.txt"
  verifier_status=$?
else
  verifier_status=99
fi
set -e
printf '{"producer_exit_status":%s,"verifier_exit_status":%s}\n' \
  "$producer_status" "$verifier_status" > "$N83_OUTPUT_DIR/exit-status.json"
date -u +%FT%TZ > "$N83_OUTPUT_DIR/finished_utc.txt"
test "$producer_status" -eq 0 && test "$verifier_status" -eq 0
