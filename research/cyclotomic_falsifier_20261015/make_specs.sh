#!/usr/bin/env bash
# Writes the frozen cost-leg specs of PROTOCOL.md §3. Thin orchestration only.
# Methods: ic.pipeline and rho.signed_frobenius (operations only; amendment A2: no full_prices.json row for a=1 n=23).
# K_full from `ecbench fb`: a=1 n=23 cyclotomic:w=2 -> 30 columns; a=0 n=23 -> 17.
set -euo pipefail
OUT="$(cd "$(dirname "$0")" && pwd)/specs"
ic() { # name role factor_base
  jq -n --arg name "$1" --arg role "$2" --arg fb "$3" \
    '{name:$name, role:$role, method:{id:"ic.pipeline", params:{factor_base:$fb, oracle:"mitm-frobenius-counted:m=2", max_trials:"1000000000", targets:"stride-renew"}}}'
}
spec() { # a n seed arms...
  local a=$1 n=$2 seed=$3; shift 3
  local arms='[{"name":"rho-frob","role":"reference","method":{"id":"rho.signed_frobenius"}}]'
  while [ $# -gt 0 ]; do arms=$(jq --argjson x "$(ic "$1" "$2" "$3")" '. + [$x]' <<<"$arms"); shift 3; done
  jq -n --argjson arms "$arms" --argjson a "$a" --argjson n "$n" --argjson seed "$seed" '{
    schema:"ecbench.spec/v1",
    label:"D4 falsifier cost leg: cyclotomic vs compact-orbit-scan at equal columns, Koblitz a=\($a) n=\($n)",
    description:"Frozen in research/cyclotomic_falsifier_20261015/PROTOCOL.md section 3 before any of its data. Operation counts only; L0; not a speed claim.",
    workloads:{curves:[{kind:"koblitz", a:$a, n:$n}], targets_per_curve:16, target_seed:$seed, target_kind:"public"},
    arms:$arms,
    measurement:{rounds:2, warmup:1, order:"alternate", seed:($seed + 1000), isolation_required:"L0", timeout_seconds:3600}}' > "$OUT/cost-a$a-n$n.json"
  echo "$OUT/cost-a$a-n$n.json"
}
C=raw_x_cap=100000000
spec 1 23 202610151001 \
  ic-cyc-w1 candidate "cyclotomic:w=1" \
  ic-cos-c1 baseline "compact-orbit-scan:columns=1,$C" \
  ic-cyc-w2 candidate "cyclotomic:w=2" \
  ic-cos-c30 baseline "compact-orbit-scan:columns=30,$C"
spec 0 23 202610151003 \
  ic-cyc-w2 candidate "cyclotomic:w=2" \
  ic-cos-c17 baseline "compact-orbit-scan:columns=17,$C" \
  ic-cyc-w3-c3 candidate "cyclotomic:w=3,columns=3" \
  ic-cos-c3 baseline "compact-orbit-scan:columns=3,$C"
