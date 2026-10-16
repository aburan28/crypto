#!/bin/bash
# Round 3: pursue the arity+size (heavy-delta symmetric-S4) lever toward the
# challenge degree. n=29 has ord_29(2)=28=n-1 (131's full-length order).
# Yield floor 2^(3l-n) decays with delta, so cells also probe where the
# trade lands. TMPDIR redirected off the full sage temp volume.
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

# symmetric-S4 delta ladder at fixed n=21 (subfield exists for l=3)
run "sym4_n21_l2" --field 21 --dim 2 --symmetric-s4 --structures geometric,standard,random --targets 16 --dmax 5 --dstar-sample 6
run "sym4_n21_l4" --field 21 --dim 4 --symmetric-s4 --structures geometric,standard,random --targets 16 --dmax 6 --dstar-sample 6

# heavy-delta cells at larger fields
run "sym4_n19_l2" --field 19 --dim 2 --symmetric-s4 --structures geometric,standard,random --targets 16 --dmax 5 --dstar-sample 6
run "sym4_n23_l3" --field 23 --dim 3 --symmetric-s4 --structures geometric,standard,random --targets 16 --dmax 5 --dstar-sample 6
run "sym4_n23_l2" --field 23 --dim 2 --symmetric-s4 --structures geometric,standard,random --targets 32 --dmax 5 --dstar-sample 8
run "sym4_n29_l3" --field 29 --dim 3 --symmetric-s4 --structures geometric,standard,random --targets 32 --dmax 5 --dstar-sample 8
run "sym4_n29_l4" --field 29 --dim 4 --symmetric-s4 --structures geometric,standard,random --targets 16 --dmax 6 --dstar-sample 6

# complete the n=23 m=2 grid (top dims)
for seed in 1 2; do
  run "m2_n23_d${d:-9}_s${seed}" --field 23 --dim 9 --structures standard,random,nb_window,nb_break --targets 16 --seed $seed --dmax 6
  run "m2_n23_d10_s${seed}" --field 23 --dim 10 --structures standard,random,nb_window,nb_break --targets 16 --seed $seed --dmax 6
done

echo "round3 done"
