#!/usr/bin/env bash
# Writes the frozen claim specs (one session per curve). Thin orchestration only.
# Usage: bash research/ecbench_frobenius_orbit_vs_rho_20261014/make_specs.sh
set -euo pipefail
ROOT="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
OUT="$ROOT/docs/ecbench/specs"
# a n columns negfold-dimension(or -) target_seed
PANEL="0 13 1 4 202610142001
1 17 1 6 202610142003
0 19 2 6 202610142005
0 23 3 8 202610142007
1 29 1 - 202610142009
0 31 2 - 202610142011
0 37 8 - 202610142013
0 41 90 - 202610142015"
ic() { # name role factor_base oracle targets
  jq -n --arg name "$1" --arg role "$2" --arg fb "$3" --arg or "$4" --arg t "$5" \
    '{name:$name, role:$role, method:{id:"ic.pipeline_full", params:{factor_base:$fb, oracle:$or, max_trials:"1000000000", targets:$t}}}'
}
while read -r a n c d seed; do
  orbit="compact-orbit-scan:columns=$c,raw_x_cap=100000000"
  arms=$(jq -n '[
    {name:"rho-frob", role:"reference", method:{id:"rho.signed_frobenius_full"}},
    {name:"rho-strong", role:"baseline", method:{id:"rho.signed_frobenius_strong_full"}}]')
  arms=$(jq --argjson x "$(ic "ic-orbit-c$c-renew" candidate "$orbit" "mitm-frobenius-counted:m=2" stride-renew)" '. + [$x]' <<<"$arms")
  arms=$(jq --argjson x "$(ic "ic-orbit-c$c-walk" baseline "$orbit" "mitm-frobenius-counted:m=2" walk)" '. + [$x]' <<<"$arms")
  if [ "$d" != "-" ]; then
    arms=$(jq --argjson x "$(ic "ic-negfold-d$d-renew" baseline "koblitz-standard-subspace:dimension=$d" "pdp3-koblitz:m=2,engine=rootfind-negfold" stride-renew)" '. + [$x]' <<<"$arms")
  fi
  arms=$(jq '. + [{name:"rho-frob-aa", role:"control", method:{id:"rho.signed_frobenius_full"}}]' <<<"$arms")
  jq -n --argjson arms "$arms" --argjson a "$a" --argjson n "$n" --argjson seed "$seed" '{
    schema:"ecbench.spec/v1",
    label:"Frobenius-orbit IC, fully priced, vs rho: Koblitz a=\($a) n=\($n)",
    description:"Frozen in research/ecbench_frobenius_orbit_vs_rho_20261014/PROTOCOL.md before any of its data. Every arm is a *_full method: all native counters priced from docs/ecbench/full_prices.json. Engineering; no security claim.",
    workloads:{curves:[{kind:"koblitz", a:$a, n:$n}], targets_per_curve:16, target_seed:$seed, target_kind:"public"},
    arms:$arms,
    measurement:{rounds:3, warmup:1, order:"alternate", seed:($seed + 1000), isolation_required:"L0", timeout_seconds:3600}}' \
    > "$OUT/claim-koblitz-frob-orbit-full-a$a-n$n.json"
  echo "$OUT/claim-koblitz-frob-orbit-full-a$a-n$n.json"
done <<<"$PANEL"
