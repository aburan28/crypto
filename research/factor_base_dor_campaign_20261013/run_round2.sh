#!/bin/bash
# Round 2: scaling the m=2 grid toward n=131-fidelity primes (n=17,19,23:
# x^n-1 = (x+1)*irreducible, like 131; ord_19(2)=18=n-1 like ord_131(2)=130),
# delta sweeps at larger n, and symmetric-S4 dim/subfield cells (arity+size).
set -u
cd "$(dirname "$0")"
SAGE=../../sage/factor_base_dor.sage
RAW=raw
mkdir -p "$RAW"
run() {
  local name="$1"; shift
  if [ -s "$RAW/$name.json" ]; then echo "skip $name"; return; fi
  echo "=== $name : $*"
  sage "$SAGE" "$@" --out "$RAW/$name.json" > "$RAW/$name.log" 2>&1 \
    || echo "FAILED $name (see $RAW/$name.log)"
}

# m=2 grid at n=17, 19, 23 (2 seeds, 16 targets)
for seed in 1 2; do
  for d in 2 3 4 5 6 7 8; do
    run "m2_n17_d${d}_s${seed}" --field 17 --dim $d \
      --structures standard,random,nb_window,nb_break --targets 16 --seed $seed --dmax 6
  done
  for d in 2 3 4 5 6 7 8 9; do
    run "m2_n19_d${d}_s${seed}" --field 19 --dim $d \
      --structures standard,random,nb_window,nb_break --targets 16 --seed $seed --dmax 6
  done
  for d in 2 3 4 5 6 7 8; do
    run "m2_n23_d${d}_s${seed}" --field 23 --dim $d \
      --structures standard,random,nb_window,nb_break --targets 16 --seed $seed --dmax 6
  done
done

# delta sweeps at larger n (standard, m=2)
run "delta_n15_m2" --field 15 --structures standard --delta-sweep=-1,1,3,5 --dmax 6 --targets 16
run "delta_n17_m2" --field 17 --structures standard --delta-sweep=-1,1,3,5,7 --dmax 6 --targets 16
run "delta_n19_m2" --field 19 --structures standard --delta-sweep=-1,1,3,5,7,9 --dmax 6 --targets 16

# symmetric S4: arity+size cells
run "sym4_n21_l3" --field 21 --dim 3 --symmetric-s4 --structures geometric,standard,subfield,random --targets 16 --dmax 5 --dstar-sample 6
run "sym4_n19_l3" --field 19 --dim 3 --symmetric-s4 --structures geometric,standard,random --targets 16 --dmax 5 --dstar-sample 6
run "sym4_n15_l2" --field 15 --dim 2 --symmetric-s4 --structures geometric,standard,random --targets 16 --dmax 5 --dstar-sample 6

echo "round2 done"
