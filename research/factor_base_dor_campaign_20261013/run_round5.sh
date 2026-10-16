#!/bin/bash
# Round 5: (a) l=1 cells - does yield persist at the smallest base (D*=2)?
# (b) fine l sweep at fixed n=29/37; (c) standard-base obstruction
# diagnostic handled by analyze_standard.sage (no engine change).
set -u
cd "$(dirname "$0")"
export TMPDIR=/Volumes/SSD990-2/tmp
mkdir -p "$TMPDIR"
SAGE=../../sage/factor_base_dor.sage
RAW=raw
run() {
  local name="$1"; shift
  if [ -s "$RAW/$name.json" ]; then echo "skip $name"; return; fi
  rm -f "$RAW/$name.json"
  echo "=== $name : $*"
  sage "$SAGE" "$@" --out "$RAW/$name.json" > "$RAW/$name.log" 2>&1 \
    || echo "FAILED $name (see $RAW/$name.log)"
}

# l=1: smallest possible base per summand
run "sym4_n23_l1" --field 23 --dim 1 --symmetric-s4 --structures geometric,standard,random --targets 512 --dmax 4 --dstar-sample 4
run "sym4_n29_l1" --field 29 --dim 1 --symmetric-s4 --structures geometric,standard,random --targets 512 --dmax 4 --dstar-sample 4
run "sym4_n37_l1" --field 37 --dim 1 --symmetric-s4 --structures geometric,standard,random --targets 512 --dmax 4 --dstar-sample 4

# fine l sweep at fixed n (fill l=4 for n=29 already done; l=4 at n=37)
run "sym4_n37_l4" --field 37 --dim 4 --symmetric-s4 --structures geometric,random --targets 64 --dmax 6 --dstar-sample 4

echo "round5 done"
