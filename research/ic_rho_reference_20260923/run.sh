#!/bin/sh
# Reproduces research/ic_rho_reference_20260923: the matched rho reference
# of RESEARCH_IC_BOUNDARY_LEDGER.md §18.  Run from the repository root
# after `cargo build --release --bin ic`.  `ic` never overwrites a report,
# so each step writes a new file; move the frozen ones aside to rerun.
set -e
IC=./target/release/ic
D=research/ic_rho_reference_20260923

# 1. Calibration (§18.2): the jump count, on curves no evaluation uses
#    (seed 0xCA11B = 827675), J in {4, 8, 16} for both tuned walks.
$IC rho --prime-bits 12,14,16,18,20,22,24,26 --runs 128 --seed 827675 \
    --jumps 4,8,16 --out $D/calibration/prime.json
$IC rho --char2-degrees 13,15,17,19,21,23,25,27 --runs 128 --seed 827675 \
    --jumps 4,8,16 --out $D/calibration/char2.json

# 2. Evaluation ladder: fresh curves (seed 0xE7A1 = 59297), every prime
#    curve generated, the three walks paired, J by the calibrated rule.
$IC rho --prime-bits 12,14,16,18,20,22,24,26 --generated-primes \
    --char2-degrees 13,15,17,19,21,23,25,27 --runs 128 --seed 59297 \
    --out $D/evaluation/ladder.json

# 3. Re-pricing: every frozen report whose vs-rho the scoreboard or the
#    ledger quotes on a prime or random binary curve.  Each replays the
#    frozen walk on the recorded seeds first and fails if any run differs.
for f in ic-boundary-ledger-round5-2026-09-22 \
         ic-boundary-ledger-round5-headline-2026-09-22 \
         ic-boundary-ledger-round5-holdout-2026-09-22; do
    $IC rho --reprice docs/ic/runs/$f.json --out $D/reprice/$f.json
done
for w in W13-s20260922 W13-s20260923 W13-s20260924 \
         W15-s20260922 W15-s20260923 W15-s20260924 W17-s20260922; do
    $IC rho --reprice research/ic_framework_engines_20260922/results/baseline_v1/$w.json \
        --out $D/reprice/$w.json
done
