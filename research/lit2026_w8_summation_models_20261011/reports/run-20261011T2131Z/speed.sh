#!/bin/sh
# Thin orchestration: interleave the original binary (harness built against the
# program-baseline library 4b6e09b78, `sym` arm only) with the candidate binary
# (`sym` = previous, `eg` = candidate, interleaved in-process), ABAB with the
# order flipped every round.  Each process runs one round over every target.
ORIG=$1
CAND=$2
OUT=$3
CELLS=$4
ROUNDS=${5:-5}
mkdir -p "$OUT"
export KIC_F4_MAX_DEGREE_ONLY=0
i=1
while [ $i -le $ROUNDS ]; do
  if [ $((i % 2)) -eq 1 ]; then order="orig cand"; else order="cand orig"; fi
  for who in $order; do
    echo "round $i $who $(date -u +%FT%TZ) $(uptime)"
    if [ $who = orig ]; then
      "$ORIG" eg --rounds 1 --targets 8 --only "$CELLS" --no-degrees --arms sym --out "$OUT/r$i-orig.json" 2> "$OUT/r$i-orig.stderr"
    else
      "$CAND" eg --rounds 1 --targets 8 --only "$CELLS" --no-degrees --arms sym,eg --out "$OUT/r$i-cand.json" 2> "$OUT/r$i-cand.stderr"
    fi
    echo "exit $?"
  done
  i=$((i + 1))
done
