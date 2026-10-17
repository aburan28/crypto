#!/usr/bin/env bash
# PS-1: algebra spec, bit-degree, root-test formulations (protocols/PS-1.md).
# None of these is a timing measurement; they can run on a busy machine.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
OUT=research/pdp_oracle_speedup_20261011/ps1_algebra/raw
B=${PS1_BIN:-target/release/examples/pdp_speedup}
C=research/pdp_oracle_speedup_20261011/baseline_v1/corpus
mkdir -p "$OUT"; : > "$OUT"/degrees.jsonl; : > "$OUT"/torsion.jsonl; : > "$OUT"/agree.jsonl
: > "$OUT"/sample.jsonl; : > "$OUT"/bscheck.jsonl; : > "$OUT"/corpus_check.jsonl; : > "$OUT"/gates.jsonl

# 1. bit-degree of the quartic's coefficients (joint in x1,x2 and in x2 alone with x1 fixed)
for l in 3 4 5 6 7 8; do "$B" degrees --l "$l" --xr 6 | tee -a "$OUT/degrees.jsonl"; done

# 2. translation by the 2-torsion point is x -> 1/x (used by the exactness theorem)
for n in 5 6 7 8 9 10 11 12 24 27 30 33 36; do "$B" torsion --n "$n" | tee -a "$OUT/torsion.jsonl"; done

# 3. exhaustive agreement of the group-law formulation with the summation-polynomial
#    one: every field of degree 5..9, every l <= n, every rational target
for n in 5 6 7 8 9; do for l in $(seq 1 "$n"); do "$B" agree --n "$n" --l "$l" | tee -a "$OUT/agree.jsonl"; done; done

# 4. the bit-sliced 64-lane pipeline against the exact reference, per pair, same range
for n in 5 6 7 8 9; do for l in $(seq 1 "$n"); do "$B" bs-check --n "$n" --l "$l" | tee -a "$OUT/bscheck.jsonl"; done; done

# 5. sampled agreement at n = 11, 13 (10^6 random pairs, random rational targets)
for nl in "11 3" "11 4" "11 5" "13 4" "13 5" "13 6"; do set -- $nl; "$B" sample --n "$1" --l "$2" --pairs 1000000 | tee -a "$OUT/sample.jsonl"; done

# 6. the exact reference against the PS-0 corpus root sets, all rungs
for l in 8 9 10 11 12; do "$B" corpus-check --l "$l" --corpus "$C/corpus_l$l.txt" | tee -a "$OUT/corpus_check.jsonl"; done

# 7. gate counts per rung (needs PS-0's opcount/counts.jsonl for the baseline formulation)
for l in 8 9 10 11 12; do "$B" gates --l "$l" | tee -a "$OUT/gates.jsonl"; done
