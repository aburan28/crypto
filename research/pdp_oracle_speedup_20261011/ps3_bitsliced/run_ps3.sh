#!/usr/bin/env bash
# PS-3: bit-sliced CPU kernel (protocols/PS-3.md).  Correctness, counts and programs run on a busy
# machine; the interleaved timing stage needs a QUIET one (tools/isolated_bench.py).
#   PS3_STAGES = all | correctness | timing | threads
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
R=research/pdp_oracle_speedup_20261011
OUT=$R/ps3_bitsliced/raw
C=$R/baseline_v1/corpus
SCR=${PS3_SCRATCH:-/tmp/ps3_bin}
CPU=${PS3_CPU:-3}
ROUNDS=${PS3_ROUNDS:-5}
STAGES=${PS3_STAGES:-all}
mkdir -p "$OUT" "$SCR"
declare -A TARGETS=( [8]=16 [9]=16 [10]=16 [11]=16 [12]=6 [13]=3 )
declare -A PASSES=(  [8]=40 [9]=10 [10]=3  [11]=1  [12]=1  [13]=1 )   # the baseline's frozen passes (PS-0); l = 13 one
ARMS="baseline baseline2 group-scalar bs64 bs256 bs512"            # baseline2 = the baseline again: the A/A of every round
S="$SCR/pdp_speedup"
BUSY="python3 tools/isolated_bench.py busy --wait --"   # untimed heavy work takes the benchmark lock, so it never overlaps a timed run

if [ "$STAGES" = all ] || [ "$STAGES" = correctness ]; then
  python3 tools/isolated_bench.py busy --wait -- cargo build --release --example pdp_speedup
  cp target/release/examples/pdp_speedup "$S"
  sha256sum "$S" > "$OUT/binary_hashes.txt"
  git rev-parse HEAD > "$OUT/revision.txt"; git status --short >> "$OUT/revision.txt"

  # 1. identical verdict digest and root sets on the whole PS-0 corpus, every width, 1 and 4 threads
  : > "$OUT/verify_arms.jsonl"
  for l in 8 9 10 11 12; do
    for arm in bs64 bs256 bs512 bs64-t4 bs256-t4 bs512-t4; do
      $BUSY "$S" verify-arm --arm "$arm" --l "$l" --corpus "$C/corpus_l$l.txt" | tee -a "$OUT/verify_arms.jsonl"
    done
  done

  # 2. the l = 13 sample (timing only): verdicts decided against the baseline arm on the same targets
  if [ ! -s "$OUT/corpus_l13.txt" ]; then
    $BUSY "$S" corpus --l 13 --out "$OUT/corpus_l13.txt" --uniform 3 --planted 0 --threads 3 \
      | tee "$OUT/corpus_l13.log"
  fi
  : > "$OUT/verify_l13.jsonl"
  for arm in baseline group-scalar bs64 bs256 bs512 bs512-t4; do
    $BUSY "$S" verify-arm --arm "$arm" --l 13 --corpus "$OUT/corpus_l13.txt" | tee -a "$OUT/verify_l13.jsonl"
  done

  # 3. the recorded program is still PS-2's (same hash on every rung); the vector lanes run the same function
  : > "$OUT/programs.jsonl"
  for l in 8 9 10 11 12; do $BUSY "$S" export-program --l "$l" --tasks 8 | tee -a "$OUT/programs.jsonl"; done
  $BUSY "$S" export-program --l 11 --tasks 64 | tee -a "$OUT/programs.jsonl"
  grep -o '"sha256":"[0-9a-f]*"' "$R/ps2_scalar/raw/programs.jsonl" > "$OUT/programs_ps2.sha"
  grep -o '"sha256":"[0-9a-f]*"' "$OUT/programs.jsonl" > "$OUT/programs_ps3.sha"
  if diff -q "$OUT/programs_ps2.sha" "$OUT/programs_ps3.sha" >/dev/null; then echo "program hashes equal PS-2: true" | tee "$OUT/programs_equal_ps2.txt"
  else echo "program hashes equal PS-2: FALSE" | tee "$OUT/programs_equal_ps2.txt"; exit 1; fi

  # 4. slots, valid pairs and lane utilisation per width
  : > "$OUT/utilisation.jsonl"
  for l in 8 9 10 11 12 13; do $BUSY "$S" bs-util --l "$l" | tee -a "$OUT/utilisation.jsonl"; done

  # 5. the product on its own and the batch pipeline, per width (diagnostic, not the headline)
  : > "$OUT/micro.jsonl"
  for l in 9 11; do $BUSY "$S" bs-micro --l "$l" --iters 3000 | tee -a "$OUT/micro.jsonl"; done
fi

[ "$STAGES" = correctness ] && exit 0
# 6. interleaved timing, one core: every round runs every arm once, in the order above.
[ -f "$OUT/iso.jsonl" ] || : > "$OUT/iso.jsonl"; [ -f "$OUT/times.jsonl" ] || : > "$OUT/times.jsonl"
timed_run() {   # tag arm l round cpus corpus
  local tag=$1 arm=$2 l=$3 r=$4 cpus=$5 corp=$6 label="$1-l$3-r$4" extra=()
  grep -q "\"label\":\"$label\"" "$OUT/times.jsonl" && return 0
  if [ "$arm" = baseline ]; then extra=(--passes "${PASSES[$l]}"); else extra=(--passes 1 --min-secs 3); fi
  for try in $(seq 1 200); do
    if python3 tools/isolated_bench.py run --cpus "$cpus" --wait --out "$OUT/iso.jsonl" --label "$label" -- \
         "$S" time --arm "$arm" --l "$l" --corpus "$corp" --targets "${TARGETS[$l]}" \
         "${extra[@]}" --label "$label" --out "$OUT/times.jsonl"; then return 0; fi
    echo "$label attempt $try refused to start (machine busy)" >> "$OUT/refused_starts.log"; sleep 3
  done
  echo "$label never started" >&2; return 1
}
if [ "$STAGES" = all ] || [ "$STAGES" = timing ]; then
  for l in 8 9 10 11 12 13; do
    corp=$C/corpus_l$l.txt; [ "$l" = 13 ] && corp=$OUT/corpus_l13.txt
    for r in $(seq 1 "$ROUNDS"); do
      for tag in $ARMS; do
        arm=$tag; [ "$tag" = baseline2 ] && arm=baseline
        timed_run "$tag" "$arm" "$l" "$r" "$CPU" "$corp"
      done
    done
  done
  {
    echo "## A/A: baseline vs baseline (ns per target)"; "$S" summarize --times "$OUT/times.jsonl" --iso "$OUT/iso.jsonl" --a baseline --b baseline2 --metric ns_per_target
    for b in group-scalar bs64 bs256 bs512; do
      echo; echo "## baseline vs $b (ns per target, then ns per own pair)"
      "$S" summarize --times "$OUT/times.jsonl" --iso "$OUT/iso.jsonl" --a baseline --b "$b" --metric ns_per_target
      "$S" summarize --times "$OUT/times.jsonl" --iso "$OUT/iso.jsonl" --a baseline --b "$b" --metric ns_per_pair
    done
    for b in bs64 bs256 bs512; do
      echo; echo "## group-scalar vs $b (ns per pair)"
      "$S" summarize --times "$OUT/times.jsonl" --iso "$OUT/iso.jsonl" --a group-scalar --b "$b" --metric ns_per_pair
    done
  } > "$OUT/../ab_summary.md"
  {
    echo "## contended runs excluded"
    for b in group-scalar bs64 bs256 bs512; do
      echo; echo "## baseline vs $b, contended rounds excluded (ns per target)"
      "$S" summarize --times "$OUT/times.jsonl" --iso "$OUT/iso.jsonl" --a baseline --b "$b" --metric ns_per_target --exclude-contended
    done
    echo; echo "## A/A, contended rounds excluded"
    "$S" summarize --times "$OUT/times.jsonl" --iso "$OUT/iso.jsonl" --a baseline --b baseline2 --metric ns_per_target --exclude-contended
  } > "$OUT/../ab_summary_excl_contended.md"
fi

# 7. threads: the same arm on 1, 2 and 4 pinned cores, interleaved (efficiency = per-core throughput ratio).
if [ "$STAGES" = all ] || [ "$STAGES" = threads ]; then
  for l in 11 12; do
    for r in $(seq 1 "$ROUNDS"); do
      timed_run "bs512-1c" bs512 "$l" "$r" 3 "$C/corpus_l$l.txt"
      timed_run "bs512-2c" bs512-t2 "$l" "$r" 2-3 "$C/corpus_l$l.txt"
      timed_run "bs512-4c" bs512-t4 "$l" "$r" 0-3 "$C/corpus_l$l.txt"
    done
  done
fi
