#!/bin/bash
# Round 4: yield-amplification precision. High target counts on the
# symmetric-S4 heavy-delta cells to fit A(n,l) = pfrac * 2^n / 2^(3l);
# D* = l+1 law extension to n=31/37/41 (ord_37(2)=36=n-1 like 131);
# m=2 grid extension to n=29.
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

# yield precision: l=2 cells with 256 targets (cheap systems)
run "sym4_n23_l2b" --field 23 --dim 2 --symmetric-s4 --structures geometric,standard,random --targets 256 --dmax 5 --dstar-sample 4
run "sym4_n29_l2"  --field 29 --dim 2 --symmetric-s4 --structures geometric,standard,random --targets 256 --dmax 5 --dstar-sample 4
run "sym4_n31_l2"  --field 31 --dim 2 --symmetric-s4 --structures geometric,standard,random --targets 256 --dmax 5 --dstar-sample 4
run "sym4_n37_l2"  --field 37 --dim 2 --symmetric-s4 --structures geometric,standard,random --targets 256 --dmax 5 --dstar-sample 4
run "sym4_n41_l2"  --field 41 --dim 2 --symmetric-s4 --structures geometric,standard,random --targets 256 --dmax 5 --dstar-sample 4

# yield precision: l=3 cells with 128 targets
run "sym4_n29_l3b" --field 29 --dim 3 --symmetric-s4 --structures geometric,standard,random --targets 128 --dmax 5 --dstar-sample 4
run "sym4_n37_l3"  --field 37 --dim 3 --symmetric-s4 --structures geometric,standard,random --targets 128 --dmax 5 --dstar-sample 4

# m=2 grid extension to n=29 (standard/random only, capped dims)
for seed in 1 2; do
  for d in 4 5 6 7 8 9; do
    run "m2_n29_d${d}_s${seed}" --field 29 --dim $d \
      --structures standard,random --targets 16 --seed $seed --dmax 6
  done
done

echo "round4 done"
