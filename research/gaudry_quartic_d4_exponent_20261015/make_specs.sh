#!/usr/bin/env bash
# Writes specs/SPEC-<job>.json and jobs.txt (the queue, in run order) for
# PROTOCOL.md. Thin orchestration only: every number is ecbench's.
set -euo pipefail
cd "$(dirname "$0")"
ORACLE='gaudry-d4:macaulay=26,weighted=1,chains=1,sparsest=1,replay=1,charge_unsolved=1'
MSEED=202610150
curve() { printf '{"kind":"fp4_torsion","p":%d,"seed":1,"max_cofactor":16}' "$1"; }
ic_arm() { # name oracle-suffix
  printf '{"name":"%s","role":"candidate","method":{"id":"ic.pipeline","params":{"factor_base":"fghr-y-line:zero=1","oracle":"%s%s","linalg":"wiedemann","max_trials":"400000"}}}' "$1" "$ORACLE" "$2"
}
rho_arm='{"name":"rho-neg","role":"reference","method":{"id":"rho.negation"}}'
spec() { # job p targets target_seed arms-json label
  cat > "specs/SPEC-$1.json" <<J
{"schema":"ecbench.spec/v1","label":"$6","description":"research/gaudry_quartic_d4_exponent_20261015/PROTOCOL.md, job $1","workloads":{"curves":[$(curve "$2")],"targets_per_curve":$3,"target_seed":$4},"arms":[$5],"measurement":{"rounds":1,"warmup":0,"order":"alternate","seed":$MSEED,"isolation_required":"L0","timeout_seconds":172800}}
J
  echo "$1" >> jobs.txt
}
tseed() { echo $((2026101500000 + $1 * 10 + $2)); }
: > jobs.txt
main() { spec "main-p$1-t$2" "$1" 1 "$(tseed "$1" "$2")" "$rho_arm,$(ic_arm ic-d4-w '')" "k=4 IC (D4 oracle, Wiedemann) and rho on E(F_$1^4), target $2"; }
lp2()  { spec "lp2-p$1-t$2"  "$1" 1 "$(tseed "$1" "$2")" "$rho_arm,$(ic_arm ic-d4-lp2 ',large_primes=2')" "k=4 IC double large prime and rho on E(F_$1^4), target $2"; }
lp1()  { spec "lp1-p$1-t$2"  "$1" 1 "$(tseed "$1" "$2")" "$rho_arm,$(ic_arm ic-d4-lp1 ',large_primes=1')" "k=4 IC single large prime and rho on E(F_$1^4), target $2"; }
rho()  { spec "rho-p$1" "$1" 16 $((2026101590000 + $1)) "$rho_arm" "rho on E(F_$1^4), 16 targets"; }
# Run order (longest first within each tier; see PROTOCOL.md "Budget").
for p in 1009 701 509 401; do main $p 0; main $p 1; done
rho 1009
for p in 293 197; do lp2 $p 0; lp2 $p 1; done
for p in 101 61; do lp1 $p 0; lp1 $p 1; done
for p in 293 197 149 101 61 53 37 29; do main $p 0; main $p 1; done
for p in 149 101 61 53 37 29; do lp2 $p 0; lp2 $p 1; done
for p in 53 37 29; do lp1 $p 0; lp1 $p 1; done
for p in 293 197 149 101 61 53 37 29; do main $p 2; done
for p in 701 509 401 293 197 149 101 61 53 37 29; do rho $p; done
lp1 149 0; lp1 149 1
