#!/bin/sh
# D2 substitute for prime_fast (descriptive): `step_pf64` is never called, so
# no binary links it; disassemble the production rho's per-thread loop,
# prime_fast::rho_worker (and its closure), with the same counts as disasm.sh.
#   disasm_rho_worker.sh BINARY LABEL OUT_DIR   (appends to OUT_DIR/rho_worker_counts.csv)
set -eu
bin=$1; label=$2; out=$3
mkdir -p "$out"
[ -f "$out/rho_worker_counts.csv" ] || echo "label,function,instructions,jcc,cmov,vector,mulx" > "$out/rho_worker_counts.csv"
for sym in 'rho_worker' 'rho_worker::{closure#0}'; do
  line=$(nm -C -S "$bin" | grep -F " crypto_lib::cryptanalysis::prime_fast::$sym" | grep -vF "$sym::" | head -1)
  [ -n "$line" ] || { echo "$label: $sym not in binary" >&2; continue; }
  addr=$(echo "$line" | awk '{print $1}'); size=$(echo "$line" | awk '{print $2}')
  name=$(echo "$sym" | tr -c 'a-z0-9_\n' '_')
  objdump -d --no-show-raw-insn -C "$bin" --start-address=0x"$addr" \
    --stop-address="$(printf '0x%x' $((0x$addr + 0x$size)))" | grep -E '^\s+[0-9a-f]+:' > "$out/$label-$name.s"
  f="$out/$label-$name.s"
  ins=$(grep -cE '^\s+[0-9a-f]+:' "$f" || true)
  jcc=$(grep -cE '^\s+[0-9a-f]+:\s+j(a|ae|b|be|e|ne|g|ge|l|le|s|ns|o|no|p|np)\s' "$f" || true)
  cmov=$(grep -cE '^\s+[0-9a-f]+:\s+cmov' "$f" || true)
  vec=$(grep -cE '%[xyz]mm' "$f" || true)
  mulx=$(grep -cE '^\s+[0-9a-f]+:\s+mulx' "$f" || true)
  echo "$label,prime_fast::$sym,$ins,$jcc,$cmov,$vec,$mulx" >> "$out/rho_worker_counts.csv"
done
