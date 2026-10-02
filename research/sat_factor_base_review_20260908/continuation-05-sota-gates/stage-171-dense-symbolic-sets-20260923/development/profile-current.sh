#!/bin/zsh
set -eu

binary=/Volumes/SSD990/kic-f4-stage171-frozen-target/release/examples/koblitz_pdp_backend
instance=/Volumes/SSD990/crypto-f4-single-target-stage171/research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-170-parallel-fixed-x1-construction-20260923/selected-run/tasks/000000-b-421e22a9c1c3b9d56396c8bbd0e46185bbebc32de0306f2bee6e9585703c2be4/instance/manifest.json
output=/Volumes/SSD990/kic-f4-stage171-dev/profile-current

"$binary" native-f4 "$instance" 300 > "$output/stdout.json" 2> "$output/backend.stderr" &
backend_pid=$!
/usr/bin/sample "$backend_pid" 10 -file "$output/sample.txt" > "$output/sample.stdout" 2> "$output/sample.stderr" || true
wait "$backend_pid"
