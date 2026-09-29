#!/usr/bin/env bash
# T3 of RESEARCH_CHAIN_SPLIT_ORDER.md §3: the oracle-pricing ladder of
# RESEARCH_INHERITED_F4.md §3.4, both arms on one `ic` binary.  Every oracle
# decides the same targets; a disagreement between any two is a failed run.
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here}/oracle_ladder"
ic="${IC:-target/release/ic}"
mkdir -p "$out"
for arm in reference candidate; do
  if [ "$arm" = reference ]; then env=(KIC_CHAIN_ORDER=layout KIC_LINEAR_ELIM=0 KIC_F4_DROP=complete)
  else env=(KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete); fi
  env "${env[@]}" "$ic" boundary --regime koblitz --koblitz-degrees 11 --repeats 1 \
    --seed 123212651130 --no-fold-max-degree 31 --s4-max-degree 31 --oracles --oracle-targets 8 \
    --json --out "$out/$arm.json" > /dev/null 2> "$out/$arm.stderr"
  echo "$arm: exit $?"
done
