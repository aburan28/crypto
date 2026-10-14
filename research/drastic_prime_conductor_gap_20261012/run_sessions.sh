#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
study=research/drastic_prime_conductor_gap_20261012
runner=${ECBENCH_BINARY:-./target/release/ecbench}
export RAYON_NUM_THREADS=1
for cell in d19-ell508433 d43-ell508433 d19-ell1016839 d43-ell1016839 d19-ell2033609 d43-ell2033609; do
  "$runner" run --wait --cpus auto --allow-busy \
    --spec "$study/specs/$cell.json" --out "$study/sessions/$cell" \
    > "$study/sessions/$cell.run.log" 2>&1
done
for cell in d19-ell508433 d43-ell508433 d19-ell1016839 d43-ell1016839 d19-ell2033609 d43-ell2033609; do
  "$runner" verify --dir "$study/sessions/$cell" --replay-all \
    --out "$study/sessions/$cell.audit.json" \
    > "$study/sessions/$cell.verify.log" 2>&1
  "$runner" compare --dir "$study/sessions/$cell" --a rho-neg --b rho-neg-aa --save \
    > "$study/sessions/$cell.aa-compare.log" 2>&1
done
./target/release/examples/drastic_gap_analyze \
  --sessions "$study/sessions" --out "$study/analysis"
