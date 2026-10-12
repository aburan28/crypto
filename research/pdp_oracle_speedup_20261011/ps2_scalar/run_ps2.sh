#!/usr/bin/env bash
# PS-2: scalar kernel by the group law (protocols/PS-2.md).  Correctness and counts can run on a
# busy machine; the interleaved timing stage needs a QUIET one (tools/isolated_bench.py).
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
OUT=research/pdp_oracle_speedup_20261011/ps2_scalar/raw
C=research/pdp_oracle_speedup_20261011/baseline_v1/corpus
SCR=${PS2_SCRATCH:-/tmp/ps2_bin}
CPU=${PS2_CPU:-3}
ROUNDS=${PS2_ROUNDS:-5}
STAGES=${PS2_STAGES:-all}        # all | correctness (builds, digests, counts, programs) | timing (needs a quiet machine)
mkdir -p "$OUT" "$SCR"
declare -A TARGETS=( [8]=16 [9]=16 [10]=16 [11]=16 [12]=6 )
declare -A PASSES=(  [8]=40 [9]=10 [10]=3  [11]=1  [12]=1 )   # the baseline's frozen passes (PS-0)
ARMS="baseline baseline2 so-full so-liftable group-scalar"      # baseline2 = the baseline again: the A/A of every round

if [ "$STAGES" != timing ]; then
  python3 tools/isolated_bench.py busy -- cargo build --release --example pdp_speedup
  CARGO_TARGET_DIR=target/opcount python3 tools/isolated_bench.py busy -- \
    cargo build --release --example pdp_speedup --features pdp-opcount
  cp target/release/examples/pdp_speedup "$SCR/pdp_speedup"
  cp target/opcount/release/examples/pdp_speedup "$SCR/pdp_speedup_opcount"
  sha256sum "$SCR/pdp_speedup" "$SCR/pdp_speedup_opcount" > "$OUT/binary_hashes.txt"

  # 1. identical verdict digest (and, for the group kernel, identical root sets) on the whole PS-0 corpus
  : > "$OUT/verify_arms.jsonl"
  for l in 8 9 10 11 12; do
    for arm in group-scalar so-liftable so-full baseline; do
      "$SCR/pdp_speedup" verify-arm --arm "$arm" --l "$l" --corpus "$C/corpus_l$l.txt" | tee -a "$OUT/verify_arms.jsonl"
    done
  done

  # 2. counted field operations per pair (instrumented build) and per quartic formation
  : > "$OUT/counts.jsonl"
  for l in 8 9 10 11 12; do
    for arm in group-scalar so-liftable; do
      "$SCR/pdp_speedup_opcount" count --l "$l" --corpus "$C/corpus_l$l.txt" --arm "$arm" --targets 2 | tee -a "$OUT/counts.jsonl"
    done
    "$SCR/pdp_speedup_opcount" count --l "$l" --corpus "$C/corpus_l$l.txt" --what formation | tee -a "$OUT/counts.jsonl"
  done

  # 3. the straight-line program of a batch of 8 tasks (hash, size, interpreter check), every rung; 64 tasks at l = 11
  : > "$OUT/programs.jsonl"
  for l in 8 9 10 11 12; do "$SCR/pdp_speedup" export-program --l "$l" --tasks 8 | tee -a "$OUT/programs.jsonl"; done
  "$SCR/pdp_speedup" export-program --l 11 --tasks 64 | tee -a "$OUT/programs.jsonl"
fi

[ "$STAGES" = correctness ] && exit 0
# 4. interleaved timing: every round runs every arm once, in the order above; two copies of the baseline give the A/A.
[ -f "$OUT/iso.jsonl" ] || : > "$OUT/iso.jsonl"; [ -f "$OUT/times.jsonl" ] || : > "$OUT/times.jsonl"
timed_run() {   # tag arm l round
  local tag=$1 arm=$2 l=$3 r=$4 label="$1-l$3-r$4" extra=()
  grep -q "\"label\":\"$label\"" "$OUT/times.jsonl" && return 0
  if [ "$arm" = baseline ]; then extra=(--passes "${PASSES[$l]}"); else extra=(--passes 1 --min-secs 3); fi
  for try in $(seq 1 100); do
    if python3 tools/isolated_bench.py run --cpus "$CPU" --wait --out "$OUT/iso.jsonl" --label "$label" -- \
         "$SCR/pdp_speedup" time --arm "$arm" --l "$l" --corpus "$C/corpus_l$l.txt" --targets "${TARGETS[$l]}" \
         "${extra[@]}" --label "$label" --out "$OUT/times.jsonl"; then return 0; fi
    echo "$label attempt $try refused to start (machine busy)" >> "$OUT/refused_starts.log"; sleep 3
  done
  echo "$label never started" >&2; return 1
}
for l in 8 9 10 11 12; do
  for r in $(seq 1 "$ROUNDS"); do
    for tag in $ARMS; do
      arm=$tag; [ "$tag" = baseline2 ] && arm=baseline
      timed_run "$tag" "$arm" "$l" "$r"
    done
  done
done
S="$SCR/pdp_speedup"
{
  echo "## A/A: baseline vs baseline (ns per target)";            "$S" summarize --times "$OUT/times.jsonl" --iso "$OUT/iso.jsonl" --a baseline --b baseline2 --metric ns_per_target
  for b in so-full so-liftable group-scalar; do
    echo; echo "## baseline vs $b (ns per target, then ns per own pair)"
    "$S" summarize --times "$OUT/times.jsonl" --iso "$OUT/iso.jsonl" --a baseline --b "$b" --metric ns_per_target
    "$S" summarize --times "$OUT/times.jsonl" --iso "$OUT/iso.jsonl" --a baseline --b "$b" --metric ns_per_pair
  done
} > "$OUT/../ab_summary.md"
cat "$OUT/../ab_summary.md"
