#!/usr/bin/env bash
# Cost leg (PROTOCOL.md §3). Thin orchestration only.
set -euo pipefail
D="$(cd "$(dirname "$0")" && pwd)"
ECB="${ECB:-/home/user/wt/cyclotomic-falsifier-bin/ecbench}"
mkdir -p "$D/sessions" "$D/audits" "$D/tables"
for c in a1-n23 a0-n23; do
  out="$D/sessions/$c"
  [ -e "$out" ] || nice -n 10 "$ECB" run --wait --cpus none --allow-busy --spec "$D/specs/cost-$c.json" --out "$out"
  "$ECB" verify --dir "$out" --replay-all --out "$D/audits/audit-$c.json" | tail -3
  "$ECB" table --dir "$out" > "$D/tables/table-$c.md"
done
cmp() { "$ECB" compare --dir "$D/sessions/$1" --a "$2" --b "$3" --save > "$D/tables/compare-$1-$3-vs-$2.txt"; }
cmp a1-n23 ic-cos-c1 ic-cyc-w1; cmp a1-n23 ic-cos-c30 ic-cyc-w2
cmp a0-n23 ic-cos-c17 ic-cyc-w2; cmp a0-n23 ic-cos-c3 ic-cyc-w3-c3
for a in ic-cyc-w1 ic-cos-c1 ic-cyc-w2 ic-cos-c30; do cmp a1-n23 rho-frob $a; done
for a in ic-cyc-w2 ic-cos-c17 ic-cyc-w3-c3 ic-cos-c3; do cmp a0-n23 rho-frob $a; done
echo cost-leg-done
