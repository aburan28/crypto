#!/bin/sh
# D2 (descriptive): disassemble the timed kernel functions of an hwcost binary
# and count instruction classes. Text processing only; no result is computed.
#   disasm.sh BINARY LABEL OUT_DIR
# Writes OUT_DIR/LABEL-<fn>.s for each function found and appends rows to
# OUT_DIR/counts.csv (label,function,instructions,jcc,cmov,vector,mulx).
set -eu
bin=$1; label=$2; out=$3
mkdir -p "$out"
[ -f "$out/counts.csv" ] || echo "label,function,instructions,jcc,cmov,vector,mulx" > "$out/counts.csv"
for fn in step_fp64 step_fp64cm step_pf64; do
  addr=$(nm -C "$bin" | grep -E " [Tt] crypto_lib::cryptanalysis::hwcost::kernel::$fn\$" | awk '{print $1}' || true)
  [ -n "$addr" ] || { echo "$label: $fn not in binary" >&2; continue; }
  objdump -d --no-show-raw-insn -C "$bin" --start-address=0x"$addr" \
    --stop-address="$(printf '0x%x' $((0x$addr + 4096)))" \
    | awk '/^[0-9a-f]+ <.*>:$/ {n++} n==1' > "$out/$label-$fn.s"
  f="$out/$label-$fn.s"
  ins=$(grep -cE '^\s+[0-9a-f]+:' "$f" || true)
  jcc=$(grep -cE '^\s+[0-9a-f]+:\s+j(a|ae|b|be|e|ne|g|ge|l|le|s|ns|o|no|p|np)\s' "$f" || true)
  cmov=$(grep -cE '^\s+[0-9a-f]+:\s+cmov' "$f" || true)
  vec=$(grep -cE '%[xyz]mm' "$f" || true)
  mulx=$(grep -cE '^\s+[0-9a-f]+:\s+mulx' "$f" || true)
  echo "$label,$fn,$ins,$jcc,$cmov,$vec,$mulx" >> "$out/counts.csv"
done
