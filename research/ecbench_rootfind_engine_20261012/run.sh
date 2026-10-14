#!/usr/bin/env bash
# The sessions of this directory, in the order they ran.  Every session
# writes under sessions/<name>; nothing here deletes or edits one.
#
#   ECBENCH=target/release/ecbench LOCK=/tmp/crypto-bench.lock bash research/ecbench_rootfind_engine_20261012/run.sh
#
# On this macOS host there is no affinity control (every run is L0) and the
# host was shared with other sessions' builds and sweeps, so the runs were
# made with --cpus none --allow-busy and the wall figures are descriptive.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
ECBENCH=${ECBENCH:-target/release/ecbench}
LOCK=${LOCK:-/tmp/crypto-bench.lock}
D=research/ecbench_rootfind_engine_20261012
mkdir -p "$D/sessions" "$D/audits" "$D/comparisons" "$D/tables" "$D/logs"
"$ECBENCH" host > "$D/host.json"

run() {
  local name=$1 spec=$2 reference=$3
  shift 3
  if [ -d "$D/sessions/$name" ]; then
    echo "== $name: exists, skipping run"
  else
    echo "== $name: run"
    nice -n 10 "$ECBENCH" run --spec "$spec" --out "$D/sessions/$name" \
      --cpus none --allow-busy --lock "$LOCK" 2>&1 | tee "$D/logs/run-$name.log"
  fi
  echo "== $name: verify --replay-all"
  nice -n 10 "$ECBENCH" verify --dir "$D/sessions/$name" --replay-all --exit-code \
    --out "$D/audits/audit-$name.json" 2>&1 | tee "$D/logs/audit-$name.log"
  echo "== $name: table"
  "$ECBENCH" table --dir "$D/sessions/$name" > "$D/tables/table-$name.md"
  for arm in "$@"; do
    echo "== $name: compare $reference vs $arm"
    "$ECBENCH" compare --dir "$D/sessions/$name" --a "$reference" --b "$arm" --save \
      2>&1 | tee -a "$D/logs/compare-$name.log"
  done
  echo "== $name: resources"
  "$ECBENCH" resources --dir "$D/sessions/$name" --out "$D/comparisons/$name.resources.json" \
    2>&1 | tee "$D/logs/resources-$name.log"
}

run smoke-koblitz docs/ecbench/specs/smoke-koblitz-rootfind.json rho-strong \
  ic-f4-m3 ic-rootfind-m3 ic-f4-m2 ic-rootfind-m2 ic-mitm-m2
run smoke-prime docs/ecbench/specs/smoke-prime-rootfind.json rho-neg \
  ic-subtract ic-mitm ic-rootfind-m2 ic-rootfind-m3
run claim-n13 docs/ecbench/specs/claim-koblitz-rootfind-n13.json rho-strong \
  rho-frob ic-f4-m3 ic-rootfind-m3 ic-f4-m2 ic-rootfind-m2 rho-strong-aa
run panel-prime docs/ecbench/specs/panel-prime-rootfind.json rho-neg \
  ic-subtract ic-mitm ic-rootfind-m2 ic-rootfind-m3 rho-neg-aa
run claim-n19 docs/ecbench/specs/claim-koblitz-rootfind-n19.json rho-strong \
  rho-frob ic-f4-m3 ic-rootfind-m3 ic-f4-m2 ic-rootfind-m2 rho-strong-aa
echo "== done"
