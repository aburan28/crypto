#!/bin/sh
# Thin orchestration: copy the built example, record host load, run all parts.
S=/private/tmp/claude-501/-Volumes-SSD990-crypto-elliptic-curve-crypto-research-88c7cc/ab9482d3-a14b-4bea-8cfa-130f6184453a/scratchpad
cp /Volumes/SSD990-2/crypto/targets/lit2026-w10/release/examples/lit2026_w10_coppersmith_bound "$S/w10bin2"
cd "$S" || exit 1
date -u +%Y%m%dT%H%MZ > run_start.txt
uptime > run_uptime_start.txt
/usr/bin/time -l ./w10bin2 --part all --max-ell 3 --json run.json > run.md 2> run.err
echo "exit $?" >> run.err
uptime > run_uptime_end.txt
date -u +%Y%m%dT%H%MZ > run_end.txt
