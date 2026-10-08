#!/bin/sh
# Reproduces research/ic_rho_koblitz_20260923: the Koblitz references of
# RESEARCH_IC_BOUNDARY_LEDGER.md §19.  Run from the repository root after
#   cargo build --release --bin ic --example koblitz_reference_prices
# `ic` never overwrites a report, so move the frozen files aside to rerun.
set -e
IC=./target/release/ic
PRICES=./target/release/examples/koblitz_reference_prices
D=research/ic_rho_koblitz_20260923
mkdir -p $D/prices $D/batch

# Provenance: the commit, a clean tree, the binaries.
{
    echo "commit $(git rev-parse HEAD)"
    echo "tree_clean $(test -z "$(git status --porcelain -- src examples Cargo.toml Cargo.lock)" && echo yes || echo NO)"
    echo "ic_sha256 $(sha256sum $IC | cut -d' ' -f1)"
    echo "prices_sha256 $(sha256sum $PRICES | cut -d' ' -f1)"
    echo "host $(nproc) cpus, $(uname -m)"
    date -u +"started %Y-%m-%dT%H:%M:%SZ"
} > $D/provenance.txt

# 1. Prices, §19.1 (b) and (c): timing, so alone on the host and one after
#    another.  Seven interleaved rounds each, K_0 throughout.
#    Declared: (b) and (c) on the thread's n = 41 base, |F| = 15,744.
$PRICES 41 15744 0 7 > $D/prices/n41-F15744.json
#    (b) on the other two curves the re-read covers.
$PRICES 53 0 0 7 > $D/prices/n53.json
$PRICES 61 0 0 7 > $D/prices/n61.json
#    Not declared: (c) on the bases of the other quoted figures, so no
#    build price is carried from one base to another.
$PRICES 41 5248 0 7 > $D/prices/n41-F5248.json
$PRICES 53 15264 0 7 > $D/prices/n53-F15264.json
for w in 6832 9760 12688 18544; do
    $PRICES 61 $w 0 7 > $D/prices/n61-F$w.json
done

# 2. Batch rho, §19.1 (a): counts, so the three curves run side by side.
#    Seed 0xBA7C4 = 764868; fresh targets every batch.
$IC rho --batch-koblitz 0/41 --batch-sizes 1,4,16,32 --batches 16 --seed 764868 \
    --out $D/batch/k0n41.json 2> $D/batch/k0n41.stderr.log &
$IC rho --batch-koblitz 0/53 --batch-sizes 1,4,16,32 --batches 16 --seed 764868 \
    --out $D/batch/k0n53.json 2> $D/batch/k0n53.stderr.log &
$IC rho --batch-koblitz 0/61 --batch-sizes 1,32 --batches 8 --seed 764868 \
    --out $D/batch/k0n61.json 2> $D/batch/k0n61.stderr.log &
wait
date -u +"finished %Y-%m-%dT%H:%M:%SZ" >> $D/provenance.txt

# 3. Added after steps 1-2 ran, before any analysis: (c) on the 16,400-point
#    base of two rows of the n = 41 ladder, so none of its rows carries a
#    build price from another base.  Timing again, alone on the host.
$PRICES 41 16400 0 7 > $D/prices/n41-F16400.json
date -u +"addendum finished %Y-%m-%dT%H:%M:%SZ" >> $D/provenance.txt
