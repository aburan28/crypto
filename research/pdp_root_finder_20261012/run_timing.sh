#!/bin/bash
# The isolated timing stage of research/pdp_root_finder_20261012 (PROTOCOL.md):
# every timed command runs alone on CPU 3 through tools/isolated_bench.py,
# which records the conditions of each run in raw/isolated.jsonl.  Five rounds,
# arms interleaved (ABAB...), an A/A pair of the baseline in every round.
# From the repository root:  research/pdp_root_finder_20261012/run_timing.sh [rounds]
set -u
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
bin=$root/target/release/examples/pdp_rootfind
ps0=$root/research/pdp_oracle_speedup_20261011/baseline_v1/corpus
iso() { python3 "$root/tools/isolated_bench.py" run --cpus 3 --out "$here/raw/isolated.jsonl" --label "$1" -- "${@:2}"; }
rounds=${1:-5}
export RAYON_NUM_THREADS=1
mkdir -p "$here/raw"
corpus() { # ℓ → corpus path
  case $1 in
    8|9|10|11|12) echo "$ps0/corpus_l$1.txt" ;;
    *) echo "$here/corpus/corpus_l$1.txt" ;;
  esac
}
for r in $(seq 1 "$rounds"); do
  echo "$(date -u +%FT%TZ) round $r" >> "$here/raw/timing_progress.txt"
  # A/A: the baseline against itself at ℓ = 10.
  iso "round$r-aa1-baseline-l10" $bin bin-time --l 10 --corpus "$(corpus 10)" --arm baseline --targets 16 --passes 1 --label "round$r-aa1-baseline-l10" --out "$here/raw/bin_time.jsonl"
  iso "round$r-aa2-baseline-l10" $bin bin-time --l 10 --corpus "$(corpus 10)" --arm baseline --targets 16 --passes 1 --label "round$r-aa2-baseline-l10" --out "$here/raw/bin_time.jsonl"
  for l in 8 9 10 11 12 13 14; do
    [ -f "$(corpus $l)" ] || continue
    k=16; p=3
    if [ "$l" -ge 13 ]; then k=4; p=1; fi
    if [ "$l" -ge 12 ]; then p=1; fi
    for arm in baseline group-scalar rf-pairs rf-mitm; do
      iso "round$r-$arm-l$l" $bin bin-time --l "$l" --corpus "$(corpus $l)" --arm "$arm" --targets "$k" --passes "$p" --label "round$r-$arm-l$l" --out "$here/raw/bin_time.jsonl"
    done
  done
  # The (23, 8) cell the Gröbner references ran on.
  for arm in baseline group-scalar rf-pairs rf-mitm; do
    iso "round$r-$arm-n23l8" $bin bin-time --l 8 --n 23 --corpus "$here/corpus/corpus_n23_l8.txt" --arm "$arm" --targets 16 --passes 3 --label "round$r-$arm-n23l8" --out "$here/raw/bin_time.jsonl"
  done
  # Prime ladders, full-scan timing on 8 uniform targets per rung.
  iso "round$r-m3" $bin prime-ladder --p 2013265921 --m 3 --t-min 4 --t-max 13 --arms pairs3,mitm3 --check-max-t 0 --alg-max-t 0 --targets 8 --passes 1 --min-secs 0.3 --label "round$r-m3" --out "$here/raw/prime_time.jsonl"
  iso "round$r-m2" $bin prime-ladder --p 2013265921 --m 2 --t-min 4 --t-max 12 --arms subtract2 --check-max-t 0 --alg-max-t 0 --targets 8 --passes 1 --min-secs 0.3 --label "round$r-m2" --out "$here/raw/prime_time.jsonl"
  iso "round$r-m4" $bin prime-ladder --p 2013265921 --m 4 --t-min 4 --t-max 11 --arms mitm4 --check-max-t 0 --alg-max-t 0 --targets 8 --passes 1 --min-secs 0.3 --label "round$r-m4" --out "$here/raw/prime_time.jsonl"
done
echo "$(date -u +%FT%TZ) done" >> "$here/raw/timing_progress.txt"
