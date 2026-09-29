#!/bin/sh
# Ledger §21.4: is the unit's code the same in both binaries?
#
# Disassembles FastCurve::add_many (the unit) and its one out-of-line
# callee, Gf2::sqr, from each `ic` binary, with addresses normalised so
# that only a change in the instructions shows in a diff.
#
#     IC_BASELINE=… IC_CANDIDATE=… sh unit_shift/disasm.sh
set -eu
here=$(dirname "$0")
for arm in baseline candidate; do
  eval bin=\$IC_$(echo $arm | tr a-z A-Z)
  for sym in 'koblitz_fast.*FastCurve.*add_many' 'semaev_decomp3Gf23sqr'; do
    line=$(nm -S --defined-only "$bin" | grep "$sym" | grep -v lazy | head -1)
    addr=$(echo "$line" | cut -d' ' -f1)
    size=$(echo "$line" | cut -d' ' -f2)
    name=$(echo "$sym" | sed 's/.*add_many/add_many/; s/.*sqr/sqr/')
    echo "$arm $name address 0x$addr size 0x$size mod64 $((0x$addr % 64))"
    objdump -d --no-show-raw-insn --start-address=0x$addr --stop-address=$((0x$addr + 0x$size)) "$bin" \
      | sed -n '8,$p' \
      | sed -E 's/^ *[0-9a-f]+:\t//; s/[0-9a-f]{6,} <[^>]*>/ADDR/g; s/0x[0-9a-f]{5,}/ADDR/g' \
      > "$here/$name-$arm.s"
  done
done
for name in add_many sqr; do
  echo "$name: $(wc -l < "$here/$name-baseline.s") instructions;" \
       "$(diff "$here/$name-baseline.s" "$here/$name-candidate.s" | grep -c '^[<>]' || true) lines differ"
done
