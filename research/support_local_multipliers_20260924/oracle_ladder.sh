#!/usr/bin/env bash
# T3 of RESEARCH_SUPPORT_LOCAL_MULTIPLIERS.md §3: the oracle-pricing ladder,
# both arms on one `ic` binary; a disagreement between any two oracles fails.
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here}/oracle_ladder"
ic="${IC:-target/release/ic}"
mkdir -p "$out"
for spec in "reference occurring" "candidate support"; do
  read -r arm multipliers <<<"$spec"
  KIC_F4_MULTIPLIERS=$multipliers KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete \
    "$ic" boundary --regime koblitz --koblitz-degrees 11 --repeats 1 \
    --seed 123212651130 --no-fold-max-degree 31 --s4-max-degree 31 --oracles --oracle-targets 8 \
    --json --out "$out/$arm.json" > /dev/null 2> "$out/$arm.stderr"
  echo "$arm: exit $?"
done
