#!/bin/sh
# Restart-safe continuation of PROTOCOL_EDWARDS.md's DLP arms (DEVIATION.md
# items 3-4): an arm whose record already ends in a final `ic` row with
# k determined or the trial cap reached is skipped; every other arm runs
# with --resume (from its checkpoint when one exists, fresh otherwise).
# Each arm is capped at 10 hours of wall time summed over its segments
# (--max-wall-hours, read from the checkpoint); a capped arm ends with an
# `ic` row marked wall_capped and is skipped thereafter.
set -eu
D=research/joux_vitse_n1_20261015/data/edwards
BIN=./target/release/examples/joux_vitse_edwards
for ps in 31:2 41:2 71:1 41:1; do
  p=${ps%:*}; s=${ps#*:}
  out="$D/dlp_p${p}_s${s}.jsonl"
  if [ -f "$out" ] && grep -q '"kind":"ic"' "$out"; then continue; fi
  $BIN dlp --p $p --seed $s --resume --max-wall-hours 10 --json "$out" > /dev/null || echo "exit $? p=$p s=$s" >> "$D/exits.txt"
done
