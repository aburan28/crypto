#!/bin/sh
D=/home/user/wt/lane-g/research/pkm_isogeny_tower_p192_20261015/data
# C3 then C4, as PROTOCOL.md (one heavy run at a time).
set -e
date -u +%Y%m%dT%H%MZ > $D/c3/started_utc.txt
nice -n 0 prlimit --as=4294967296 -- taskset -c 2 /home/user/cargo-target-lane-g/release/pkm-aux-search --prime nist-p192 --t 64 --family z2xz2 --arm fast --seed-lo 0 --seed-hi 33554432 --stop none --record-min 16 --summary $D/c3/summary.json --checkpoint $D/c3/checkpoint.json > $D/c3/stdout.json
date -u +%Y%m%dT%H%MZ > $D/c4/started_utc.txt
nice -n 0 prlimit --as=4294967296 -- taskset -c 2 /home/user/cargo-target-lane-g/release/pkm-aux-search --prime nist-p192 --t 64 --family z2xz8 --arm fast --seed-lo 0 --seed-hi 268435456 --stop none --record-min 16 --summary $D/c4/summary.json --checkpoint $D/c4/checkpoint.json > $D/c4/stdout.json
date -u +%Y%m%dT%H%MZ > $D/c4/finished_utc.txt
