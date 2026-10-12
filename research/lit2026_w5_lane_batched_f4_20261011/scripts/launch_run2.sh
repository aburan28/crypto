#!/bin/sh
# Second W5 session: orig, prev, gate, rref and gate+rref, 5 rounds, detached.
here=$(cd "$(dirname "$0")" && pwd)
OUT=/Volumes/SSD990-2/crypto/targets/lit2026-w5/run-20261012T0008Z \
ROUNDS=5 ARMS="orig prev gate rref gaterref" \
  nohup sh "$here/run_arms.sh" > /Volumes/SSD990-2/crypto/targets/lit2026-w5/arms2.log 2>&1 &
echo launched
