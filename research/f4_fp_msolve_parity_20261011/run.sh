#!/usr/bin/env bash
# f4_fp before and after this study's changes, against msolve and Singular,
# on the f4_fp_bench systems and the larger fp-xl ones.
#
#   run.sh <export_main> <export_new> <outdir> [rounds] [singular_timeout_s]
#
# <export_main>, <export_new>: release builds of examples/gb_cas_export.rs
# with src/cryptanalysis/f4_fp.rs as on main (4b6e09b78) and as on this
# branch.  Each round runs, interleaved: both builds at 1 and 4 rayon
# threads (median of 3 repetitions inside; 1 for fp-xl), msolve -g 2 at -t 1 and -t 4,
# and Singular std (which also checks the new build's basis).
set -uo pipefail
MAIN=$1
NEW=$2
OUT=$3
ROUNDS=${4:-3}
SING_TIMEOUT=${5:-600}
mkdir -p "$OUT/raw" "$OUT/inputs"
now_ms() { echo $(( $(date +%s%N) / 1000000 )); }
{
  echo "msolve $(dpkg-query -W -f='${Version}' msolve 2>/dev/null)"
  echo "singular $(dpkg-query -W -f='${Version}' singular 2>/dev/null)"
  echo "host $(lscpu | sed -n 's/^Model name: *//p') nproc $(nproc)"
  echo "export_main sha256 $(sha256sum "$MAIN" | cut -c1-64)"
  echo "export_new sha256 $(sha256sum "$NEW" | cut -c1-64)"
} > "$OUT/raw/versions.txt"

for r in $(seq 1 "$ROUNDS"); do
  for tier in fp fp-xl; do
    for arm in main new; do
      bin=$MAIN; [ "$arm" = new ] && bin=$NEW
      for th in 1 4; do
        reps=3; [ "$tier" = fp-xl ] && reps=1
        RAYON_NUM_THREADS=$th "$bin" "$tier" "$OUT/inputs" $reps > "$OUT/raw/ours_${arm}_t${th}_${tier}.round$r.jsonl"
      done
    done
    # the new build wrote the inputs last: its basis is the one Singular checks
    cases=$(sed -n 's/.*"case":"\([^"]*\)".*/\1/p' "$OUT/raw/ours_new_t1_${tier}.round$r.jsonl")
    for c in $cases; do
      for th in 1 4; do
        a=$(now_ms)
        msolve -g 2 -t "$th" -f "$OUT/inputs/$c.ms" -o "$OUT/raw/$c.msolve_t$th.out" > /dev/null 2>&1
        b=$(now_ms)
        printf '%s\tmsolve_t%s\twall_ms\t%s\n' "$c" "$th" $((b - a)) >> "$OUT/raw/external.round$r.tsv"
      done
      if out=$(timeout "$SING_TIMEOUT" Singular -q "$OUT/inputs/$c.sing" | tr '\n' ' '); then st=ok; else st="exit$?"; fi
      printf '%s\tsingular_std\t%s\t%s\n' "$c" "$st" "$out" >> "$OUT/raw/external.round$r.tsv"
    done
  done
done
gzip -9 -n -f "$OUT"/raw/*.out
