#!/bin/sh
# lit2026 W11 (a): counted field operations per visited point, per kernel.
#   run_opcount.sh <bin> <curves.jsonl> <out.jsonl> <bits...>
set -eu
bin=$1; curves=$2; out=$3; shift 3
: > "$out"
for b in "$@"; do
  line=$(grep "\"bits\":$b," "$curves")
  p=$(printf '%s' "$line" | sed -E 's/.*"p":([0-9]+).*/\1/')
  n=$(printf '%s' "$line" | sed -E 's/.*"n":([0-9]+).*/\1/')
  gx=$(printf '%s' "$line" | sed -E 's/.*"Qx":([0-9]+).*/\1/')
  gy=$(printf '%s' "$line" | sed -E 's/.*"Qy":([0-9]+).*/\1/')
  "$bin" "$p" 0 7 "$gx" "$gy" "$n" "secp${b}k1" >> "$out"
done
