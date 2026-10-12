#!/bin/bash
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
TMP=$(mktemp -d "${TMPDIR:-/tmp}/sigma-fused-star.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

CXX=${CXX:-g++}
"$CXX" -O2 -std=c++17 -Wall -Wextra -Werror "$HERE/star_log_check.cpp" -o "$TMP/log-check"
"$CXX" -O2 -std=c++17 -Wall -Wextra -Werror "$HERE/star_summarize.cpp" -o "$TMP/summarize"
ARMS=()
while IFS=$'\t' read -r arm _; do ARMS+=("$arm"); done < <("$TMP/log-check" list)
[ "${#ARMS[@]}" = 11 ]
[ "${ARMS[0]}" = baseline ]
CANDIDATES=("${ARMS[@]:1}")
DIGEST=0000000000000000000000000000000000000000000000000000000000000000
UPDATES=201863462912
printf '%s\n' 'PASS: 300/300 replay and sorted v1 corpus identity across 11 arms' > "$TMP/preflight.txt"

write_row() {
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$UPDATES" "$DIGEST" '2400, 500, 55'
}

write_panel() {
  local mode=$1 path=$2 pair arm ordinal candidateRate aRate bRate
  printf 'phase\tcomparison\tpair\torder\tvariant\tbinary\trateMps\tupdates\tlogSha256\tgpuState\n' > "$path"
  for arm in "${ARMS[@]}"; do write_row warmup "$arm" 0 1 warmup "$arm" 15000 >> "$path"; done
  for pair in 1 2 3 4 5; do
    aRate=15000; bRate=15003
    [ "$mode" = noisy ] && bRate=15300
    if [ $((pair % 2)) = 1 ]; then
      write_row aa baseline-aa "$pair" 1 a baseline "$aRate" >> "$path"
      write_row aa baseline-aa "$pair" 2 b baseline "$bRate" >> "$path"
    else
      write_row aa baseline-aa "$pair" 1 b baseline "$bRate" >> "$path"
      write_row aa baseline-aa "$pair" 2 a baseline "$aRate" >> "$path"
    fi
  done
  ordinal=0
  for arm in "${CANDIDATES[@]}"; do
    candidateRate=15150
    [ "$mode" = winner ] && [ "$arm" = pair-clmul ] && candidateRate=15450
    for pair in 1 2 3; do
      if [ $(((ordinal + pair) % 2)) = 0 ]; then
        write_row screen "$arm" "$pair" 1 baseline baseline 15000 >> "$path"
        write_row screen "$arm" "$pair" 2 candidate "$arm" "$candidateRate" >> "$path"
      else
        write_row screen "$arm" "$pair" 1 candidate "$arm" "$candidateRate" >> "$path"
        write_row screen "$arm" "$pair" 2 baseline baseline 15000 >> "$path"
      fi
    done
    ordinal=$((ordinal + 1))
  done
}

write_panel none "$TMP/none.tsv"
"$TMP/summarize" "$TMP/none.tsv" "$TMP/preflight.txt" "$TMP/none.json" > "$TMP/none.txt"
grep -q '"selected_arm": null' "$TMP/none.json"
grep -q 'SELECT NONE' "$TMP/none.txt"

write_panel winner "$TMP/winner.tsv"
"$TMP/summarize" "$TMP/winner.tsv" "$TMP/preflight.txt" "$TMP/winner.json" > "$TMP/winner.txt"
grep -q '"selected_arm": "pair-clmul"' "$TMP/winner.json"
grep -q 'SELECT pair-clmul' "$TMP/winner.txt"

sed '$d' "$TMP/winner.tsv" > "$TMP/missing.tsv"
if "$TMP/summarize" "$TMP/missing.tsv" "$TMP/preflight.txt" "$TMP/missing.json" >/dev/null 2>&1; then
  echo "summarizer admitted a missing timing row" >&2
  exit 1
fi

write_panel noisy "$TMP/noisy.tsv"
if "$TMP/summarize" "$TMP/noisy.tsv" "$TMP/preflight.txt" "$TMP/noisy.json" >/dev/null 2>&1; then
  echo "summarizer admitted excessive A/A drift" >&2
  exit 1
fi

cat > "$TMP/baseline-bench.log" <<'EOF'
packed kernel: 126 registers/thread, 0 local bytes/thread, 1792 shared bytes/block, single-product multiplier
packed launch bounds: 256 threads, 2 min blocks
packed denominator cache: 1
packed multiply by value: 1
packed Frobenius network: 3
packed polynomial chain: 1
packed polynomial state: 1
packed unrolled inversion: 1
packed paired products: 1
packed pair ilp: 0
packed pair clmul: 0
packed clmul flat: 0
packed top hoist: 0
packed onb inv: 0
packed from reduced: 0
packed inline polynomial: 3
packed slot unroll: 1
packed chains: 1
packed slot prefetch: 0
packed slot pipeline: 0
packed sigma fused: 1
packed sigma fused late y: 0
packed witness: 0
packed L2 persist: 0
packed direct reduction: 1
packed generated product: 1
packed native carryless multiply: 1
packed native carryless square: 0
packed three-limb Karatsuba: 0
packed weighted prefix: 2
packed compact state: 1
packed shared sigma: 1
packed top clmad: 0
packed state tile: 256
packed add combine: 0
packed alu square: 0
packed alu onb square: 0
packed square table: 0
packed polynomial inversion: 0
packed profile ranges: 0
packed table walk: 0 (0 branches, 1792 shared bytes)
backend cuda-packed131: 385024 threads x 16 slots x 1 lanes = 6160384 walks, dp weight 0, 1024 steps per launch
finished: 15000.000 M it/s (0 verified against the reference, 0 dropped)
EOF
[ "$("$TMP/log-check" bench baseline "$TMP/baseline-bench.log")" = 15000 ]
sed -e 's/packed L2 persist: 0/packed L2 persist: 1/' \
    -e '/packed L2 persist: 1/a\
packed L2 persist window: 83886080 of 100000000 field bytes, cap 83886080' \
    "$TMP/baseline-bench.log" > "$TMP/l2-bench.log"
[ "$("$TMP/log-check" bench l2-persist "$TMP/l2-bench.log")" = 15000 ]
grep -v 'packed L2 persist window:' "$TMP/l2-bench.log" > "$TMP/l2-no-window.log"
if "$TMP/log-check" bench l2-persist "$TMP/l2-no-window.log" >/dev/null 2>&1; then
  echo "log checker admitted L2 persist without an installed window" >&2
  exit 1
fi
sed 's/packed launch bounds: 256 threads/packed launch bounds: 128 threads/' \
  "$TMP/baseline-bench.log" > "$TMP/bad-bounds.log"
if "$TMP/log-check" bench baseline "$TMP/bad-bounds.log" >/dev/null 2>&1; then
  echo "log checker admitted wrong launch bounds" >&2
  exit 1
fi

# These are the star's arithmetic-changing arms.  L2 policy, loop unrolling
# and late-Y placement are schedule-only and are covered by the device replay.
ARITH_ARMS=(baseline pair-ilp pair-clmul from-reduced inv-poly1 inv-poly2 clmul-flat alu-square)
COMMON_DEFS=(
  -DECC_PACKED_SINGLE_PRODUCT=1 -DECC_PACKED_BY_VALUE=1
  -DECC_PACKED_PERM_SIGMA=3 -DECC_PACKED_DIRECT_REDUCE=1
  -DECC_PACKED_GENERATED_PRODUCT=1 -DECC_PACKED_INLINE_POLY=3
  -DECC_PACKED_CLMAD=1
)
for arm in "${ARITH_ARMS[@]}"; do
  read -r -a varied <<< "$("$TMP/log-check" flags "$arm")"
  defs=()
  for setting in "${varied[@]}"; do defs+=("-DECC_${setting}"); done
  "$CXX" -O2 -std=c++17 "${COMMON_DEFS[@]}" "${defs[@]}" \
    -Wno-unknown-pragmas "$ROOT/src/testpacked.cpp" -o "$TMP/testpacked-$arm"
  "$TMP/testpacked-$arm" > "$TMP/testpacked-$arm.txt"
  grep -q '^PASS: packed multiplication' "$TMP/testpacked-$arm.txt"
done

bash -n "$HERE/gpujob-star.sh"
echo "PASS: fused sigma star native contract, summarizer, marker and arithmetic tests"
