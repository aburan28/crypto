#!/bin/bash
# Build the widened worker: the round-0020 winner plus round21-wide-pair-table.patch.
#
#     bash campaign_20260916/round21_build.sh [OUTDIR]
#
# The promoted worker refuses degree 37, so the crossover measurement needs
# this one. Requires the restored round-0020 evidence, because candidate
# sources are archived rather than committed.
set -eu
cd "$(dirname "$0")/.."
BASE=runs/round-0020/source_candidates/both/source
OUT="${1:-${TMPDIR:-/tmp}/round21-wide}"
[ -d "$BASE" ] || { echo "restore round-0020 first: evidence/restore.py --archive round-0020"; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"
for f in Cargo.toml Cargo.lock build.rs .cargo/config.toml; do
  [ -f "$BASE/$f" ] && mkdir -p "$OUT/$(dirname "$f")" && cp "$BASE/$f" "$OUT/$f"
done
for d in src examples gpu; do [ -d "$BASE/$d" ] && cp -r "$BASE/$d" "$OUT/$d"; done
patch -p1 -s -d "$OUT" -i "$PWD/campaign_20260916/round21-wide-pair-table.patch"

grep -q 'MAX_DEGREE: u32 = 61' "$OUT/src/cryptanalysis/koblitz_tiny_ic.rs" || { echo "module ceiling not lifted"; exit 3; }
grep -q '(5..=61)' "$OUT/examples/ic_tournament_worker.rs" || { echo "worker dispatch bound not lifted"; exit 3; }
grep -q 'sum_x: u64' "$OUT/src/cryptanalysis/koblitz_tiny_ic.rs" || { echo "pair table not widened"; exit 3; }
echo "patched: pair table widened, both degree bounds lifted"

( cd "$OUT" && cargo build --release --example ic_tournament_worker --target x86_64-unknown-linux-musl )
W="$OUT/target/x86_64-unknown-linux-musl/release/examples/ic_tournament_worker"
echo "worker: $W"
echo "sha256: $(sha256sum < "$W" | cut -d' ' -f1)"
