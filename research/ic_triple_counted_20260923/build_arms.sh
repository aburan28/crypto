#!/usr/bin/env bash
# Rebuild both IC arms of the counted-sizing check from committed inputs.
#
#   triple  = round-0020 `both` + round21-wide-pair-table.patch + round23-scaled-base.patch
#             + ../ic_triple_table_20260923/triple-table.patch
#   counted = triple + counted-sizing.patch (this directory)
#
# Rho is the `triple` build's rho mode; neither patch touches it.
#
#   research/ic_triple_counted_20260923/build_arms.sh /absolute/work/dir
#
# Needs zstd, valgrind 3.22.0 for check.py, and the x86_64-unknown-linux-musl
# target, which the source tree's .cargo/config.toml selects.
set -euo pipefail

WORK=${1:?usage: build_arms.sh /absolute/work/dir}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
CAMPAIGN="$REPO/research/ic_candidate_tournament_20260915/campaign_20260916"
TRIPLE="$REPO/research/ic_triple_table_20260923/triple-table.patch"

mkdir -p "$WORK"
python3 "$REPO/research/ic_candidate_tournament_20260915/evidence/restore.py" \
    --archive round-0020 --out "$WORK/evidence"
BOTH="$WORK/evidence/runs/round-0020/source_candidates/both/source"

build() {
    local name=$1; shift
    rm -rf "${WORK:?}/$name"
    mkdir -p "$WORK/$name"
    (cd "$BOTH" && tar --exclude=./target -cf - .) | (cd "$WORK/$name" && tar -xf -)
    for p in "$@"; do
        (cd "$WORK/$name" && patch -p1 -s < "$p")
    done
    (cd "$WORK/$name" && cargo build --release --example ic_tournament_worker)
    echo "$name: $WORK/$name/target/x86_64-unknown-linux-musl/release/examples/ic_tournament_worker"
}

build triple  "$CAMPAIGN/round21-wide-pair-table.patch" "$CAMPAIGN/round23-scaled-base.patch" "$TRIPLE"
build counted "$CAMPAIGN/round21-wide-pair-table.patch" "$CAMPAIGN/round23-scaled-base.patch" "$TRIPLE" \
              "$HERE/counted-sizing.patch"
