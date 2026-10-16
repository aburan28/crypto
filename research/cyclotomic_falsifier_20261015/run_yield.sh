#!/usr/bin/env bash
# Yield leg (PROTOCOL.md §2). Thin orchestration only.
set -euo pipefail
D="$(cd "$(dirname "$0")" && pwd)"
BIN="${BIN:-/home/user/wt/cyclotomic-falsifier-bin/cyclotomic_falsifier_yield}"
for an in "1 23" "0 23" "1 11"; do
  set -- $an
  nice -n 10 "$BIN" --a $1 --n $2 --nulls 1999 --seed 20261015 --alpha 0.05 --family 15 \
    --out "$D/yield/yield-a$1-n$2.json" 2> "$D/yield/log-a$1-n$2.txt"
done
echo yield-done
