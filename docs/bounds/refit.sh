#!/usr/bin/env sh
# Regenerate every seed bound record and the frontier page from the committed
# sessions.  Thin orchestration only (AGENTS.md): every figure is computed by
# `ecbench`.  Run from the repository root after `cargo build --release --bin
# ecbench`.  Records are content-addressed, so an unchanged session yields an
# unchanged file; a changed one yields a new id, which `git status` shows.
set -eu
cd "$(dirname "$0")/../.."
E=${ECBENCH:-./target/release/ecbench}
R=docs/bounds/records
A=research/ecbench_all_candidates_20261003
mkdir -p "$R"

fit() { # session arm receipt
  "$E" bound fit --dir "$A/sessions/$1" --arm "$2" --audit "$A/$3" --root . \
    --out "$R/$1-$2.json" 2>&1 | grep -v '^wrote' || true
}

for arm in rho-neg rho-plain rho-frozen bsgs bsgs-il bsgs-neg kangaroo ic-mitm ic-subtract; do
  fit prime "$arm" audit-prime.json
done
for arm in rho-strong rho-frob rho-neg rho-plain rho-frozen bsgs bsgs-il bsgs-neg kangaroo ic-frob-m2 ic-frob-m3 ic-subtract; do
  fit koblitz "$arm" audit-koblitz.json
done

"$E" frontier build --bounds "$R" --out docs/bounds/frontier.json --markdown docs/bounds/FRONTIER.md
"$E" bound check --root . --record "$R"/*.json
if ls docs/bounds/challenges/*.json >/dev/null 2>&1; then
  "$E" challenge check --file docs/bounds/challenges/*.json
fi
