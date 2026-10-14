#!/bin/bash
# Thin orchestration of sage/factor_base_dor.sage (the engine does all math).
# Usage: bash run_campaign.sh
set -u
cd "$(dirname "$0")"
SAGE=../../sage/factor_base_dor.sage
RAW=raw
mkdir -p "$RAW"
run() { # run <name> <args...>
  local name="$1"; shift
  if [ -s "$RAW/$name.json" ]; then echo "skip $name"; return; fi
  echo "=== $name : $*"
  sage "$SAGE" "$@" --out "$RAW/$name.json" > "$RAW/$name.log" 2>&1 \
    || echo "FAILED $name (see $RAW/$name.log)"
}

# --- m=2 structure x dimension grid (ECC2K-130 shape: n=11,13 prime x^n-1) ---
for seed in 1 2 3; do
  for n in 11 13; do
    for d in 2 3 4 5 6; do
      if [ "$d" -le $(( (n-1)/2 )) ]; then
        dmax=6
        if [ "$n$d" = "135" ] || [ "$n$d" = "136" ]; then dmax=5; fi
        run "m2_n${n}_d${d}_s${seed}" --field $n --dim $d \
          --structures standard,random,nb_window,nb_break \
          --targets 24 --seed $seed --dmax $dmax
      fi
    done
  done
done

# --- subfield closure contrast (composite n; subfield exists) ---
for seed in 1 2 3; do
  run "m2_n9_d3_s${seed}"   --field 9  --dim 3 --structures standard,random,subfield --targets 24 --seed $seed --dmax 6
  run "m2_n15_d3_s${seed}"  --field 15 --dim 3 --structures standard,random,subfield --targets 24 --seed $seed --dmax 6
  run "m2_n15_d5_s${seed}"  --field 15 --dim 5 --structures standard,random,subfield --targets 24 --seed $seed --dmax 6
done

# --- Shantz-Teske delta sweeps (standard base, m=2) ---
run "delta_n11_m2" --field 11 --structures standard --delta-sweep=-1,0,1,2,3,4 --dmax 6 --targets 24
run "delta_n13_m2" --field 13 --structures standard --delta-sweep=-1,0,1,2,3,4 --dmax 6 --targets 24
run "delta_n11_m3" --field 11 --structures standard --summands 3 --delta-sweep=-1,0,1,2,3 --dmax 6 --targets 24 --dstar-sample 6

# --- m=3 chained Semaev (2 summands listed per chain pair) ---
for seed in 1 2 3; do
  run "m3_n11_d3_s${seed}" --field 11 --dim 3 --summands 3 --dstar-sample 6 --targets 24 --seed $seed --dmax 6
  run "m3_n13_d4_s${seed}" --field 13 --dim 4 --summands 3 --dstar-sample 6 --targets 24 --seed $seed --dmax 6
done

# --- symmetric S4 (m=3 direct encoding, Galbraith-Gebregiyorgis structures) ---
run "sym4_n9_l2"  --field 9  --dim 2 --symmetric-s4 --structures geometric,standard,subfield,random --targets 16 --dmax 5 --dstar-sample 6
run "sym4_n15_l3" --field 15 --dim 3 --symmetric-s4 --structures geometric,standard,subfield,random --targets 16 --dmax 5 --dstar-sample 6

# --- n=131 applicability checks ---
run "frobwin_n131" --frobenius-windows --field 131 --dim 4 --summands 3
run "stack_n11"    --stacking 1,2 --field 11 --dim 4
run "stack_n13"    --stacking 1,2 --field 13 --dim 4

echo "campaign done"
