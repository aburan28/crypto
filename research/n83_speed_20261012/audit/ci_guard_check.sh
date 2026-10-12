#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
command -v rg >/dev/null
if rg -n '\b(sorry|admit|sorryAx|axiom|native_decide)\b|^#eval' ScalarArithmetic.lean AxiomAudit.lean; then
  echo 'Incomplete proof or unverified evaluation is not permitted.'
  exit 1
fi
lake env lean AxiomAudit.lean | tee /tmp/n83-scalar-axioms.log
test "$(rg -c 'depends on axioms:' /tmp/n83-scalar-axioms.log)" -eq 4
unexpected=$(sed -n 's/.*depends on axioms: \[\(.*\)\]/\1/p' /tmp/n83-scalar-axioms.log |
  tr ',' '\n' | sed 's/ //g' |
  rg -v '^$|^(propext|Classical.choice|Quot.sound)$' || true)
if test -n "$unexpected"; then
  echo 'Unexpected theorem axioms:'
  printf '%s\n' "$unexpected"
  exit 1
fi
echo 'PASS: local CI source and transitive-axiom guards'
