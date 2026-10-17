#!/bin/sh
. "$(dirname "$0")/env.sh"
cd "$STUDY/calibration"
for B in 18 20 22 24 26; do
  $ECB run --wait --cpus auto --spec b$B/spec.json --out sessions/b$B > sessions/b$B.log 2>&1
  $ECB verify --dir sessions/b$B --replay 12 --out sessions/b$B.audit.json > sessions/b$B.verify.log 2>&1
done
echo ALLDONE
