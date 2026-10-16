#!/usr/bin/env bash
# Full census: every curve in data/curves.csv, m = 2 and 3, t = 2..10.
set -euo pipefail
cd "$(dirname "$0")/../../../.."
D=research/ic_next10_20261014/02-solution-count
BIN=target/release/examples/ic_next10_02_solution_count
export RAYON_NUM_THREADS=1 OMP_NUM_THREADS=1
out=$D/data/census.csv; : > $out; first=1
tail -n +2 $D/data/curves.csv | while IFS=, read bits sign p a b n; do
  prlimit --as=700000000 -- nice $BIN --p $p --a $a --b $b --n $n --tmin 2 --tmax 10 \
    --targets2 2000 --targets3 500 --seed $p > $D/data/.tmp.csv
  if [ $first = 1 ]; then head -1 $D/data/.tmp.csv | sed 's/^/bits,sign,/' >> $out; first=0; fi
  tail -n +2 $D/data/.tmp.csv | sed "s/^/$bits,$sign,/" >> $out
  echo "done $bits $sign $p" >&2
done
rm -f $D/data/.tmp.csv
