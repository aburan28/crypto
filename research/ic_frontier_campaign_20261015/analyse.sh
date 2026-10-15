#!/bin/sh
# Phase 1 analysis (PROTOCOL.md "Analysis"), thin orchestration only: every
# figure is computed by native `ecbench`.
#   1. replay-all audit of every Phase 1 session (counts must reproduce exactly)
#   2. the AGENTS.md §5a band per family, with pairwise verdicts
#   3. exponent fits per arm, tier and band weight (`bound fit --weight`)
set -u
cd "$(dirname "$0")"
E=${ECBENCH:-../../target/release/ecbench}
S=phase1/sessions
A=analysis
mkdir -p phase1/audits "$A/bounds"

for d in "$S"/*; do
  s=$(basename "$d")
  [ -f "phase1/audits/$s.json" ] || "$E" verify --dir "$d" --replay-all --out "phase1/audits/$s.json" > /dev/null 2>&1 \
    || echo "audit FAILED: $s (receipt kept)"
done

"$E" band --dir "$S"/K* --ranks --json "$A/band-koblitz.json" > "$A/band-koblitz.md"
"$E" band --dir "$S"/P20 "$S"/P24 "$S"/P28 "$S"/P32 "$S"/P36 "$S"/P40 --ranks --json "$A/band-prime.json" > "$A/band-prime.md"
"$E" band --dir "$S"/PKM* --ranks --json "$A/band-pkm.json" > "$A/band-pkm.md"

RULE="per-size parameters chosen on the train panel by PROTOCOL.md A4 (research/ic_frontier_campaign_20261015/params.json)"
fit() { # FAMILYGLOB ARM TIER WEIGHT OUT
  dirs=$(for d in $1; do grep -q "\"name\": *\"$2\"" "$d/spec.json" && echo "$d"; done)
  case "$2" in ic-compact|ic-negfold|ic-pair|ic-pkm) rule="--param-rule"; ;; *) rule=""; ;; esac
  # shellcheck disable=SC2086
  "$E" bound fit --dir $dirs --arm "$2" --tier "$3" --weight "$4" --root ../.. ${rule:+$rule "$RULE"} \
    --label "ic-frontier $2 $3 w=$4" --out "$A/bounds/$5.json" > /dev/null 2>> "$A/bounds/fit.log" \
    || echo "fit refused: $5 (see fit.log)" >> "$A/bounds/fit.log"
}
: > "$A/bounds/fit.log"
for w in 0 0.1 1; do
  for arm in rho-strong rho-frob bsgs-neg ic-compact ic-negfold ic-f4; do
    for t in toy medium; do fit "$S/K*" "$arm" "$t" "$w" "K-$arm-$t-w$w"; done
  done
  for arm in rho-neg bsgs-neg kangaroo ic-pair; do
    for t in toy medium; do fit "$S/P[0-9]*" "$arm" "$t" "$w" "P-$arm-$t-w$w"; done
  done
  for arm in rho-neg ic-pkm; do fit "$S/PKM*" "$arm" toy "$w" "PKM-$arm-toy-w$w"; done
done
echo done
