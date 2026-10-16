#!/bin/sh
# Every command behind RESULTS.md, from the repository root (single-threaded).
set -e
export CARGO_TARGET_DIR=${CARGO_TARGET_DIR:-target}
cargo build --release --example r4_16_6_genus2_census
B=$CARGO_TARGET_DIR/release/examples/r4_16_6_genus2_census
D=research/isogeny_fibers_20261010/round4/16_6_genus2_law/raw
# curve search: p = 1 mod 10, all ten twists y^2 = x^5 + b, cofactor h <= 16
for w in "1000 1200" "2000 2300" "4000 4300"; do set -- $w; $B search --lo $1 --hi $2 --maxh 16 > $D/search_p$1_$2.jsonl; done
# C2 (all five 2-torsion points rational, h = 16)
for pb in "1051 1042" "2131 32" "4111 2172"; do set -- $pb; $B census --p $1 --b $2 --size 80 --reps 40 --seed 1 --arms plain,A,M2,M4iso,M4non,G --json $D/census_c2_p$1_s80.json 2> $D/census_c2_p$1_s80.log; done
$B census --p 4111 --b 2172 --size 160 --reps 40 --seed 3 --arms plain,A,M2,M4iso,M4non --json $D/census_c2_p4111_s160.json 2> $D/census_c2_p4111_s160.log
# C1 (rational 5-torsion T5 = (0, sqrt b) - inf, h = 5)
for pb in "1051 49" "2131 4" "4111 3984"; do set -- $pb; $B census --p $1 --b $2 --size 100 --reps 40 --seed 2 --arms plain,A,M5,AM5,G --json $D/census_c1_p$1_s100.json 2> $D/census_c1_p$1_s100.log; done
# A-arm size scan on C1 (finite-size correction check)
for s in 50 200; do $B census --p 4111 --b 3984 --size $s --reps 40 --seed 4 --arms plain,A --json $D/census_c1_p4111_s${s}_Aonly.json 2> $D/census_c1_p4111_s${s}_Aonly.log; done
# plots and PDF (python only reads the recorded JSON)
python research/isogeny_fibers_20261010/round4/16_6_genus2_law/plot_law.py
