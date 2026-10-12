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
mkdir -p "$R"

# fit TOPIC SESSION ARM RECEIPT OUT [TIER]: one arm of research/TOPIC/sessions/SESSION
# into $R/OUT.json, citing the receipt research/TOPIC/RECEIPT.  TIER keeps one
# tier of a session whose sizes span two; a mixed fit is refused.
fit() {
  "$E" bound fit --dir "research/$1/sessions/$2" --arm "$3" --audit "research/$1/$4" \
    ${6:+--tier "$6"} --root . --out "$R/$5.json" 2>&1 | grep -v '^wrote' || true
}

# Every candidate (research/ecbench_all_candidates_20261003): four sizes per
# family, all toy.
A=ecbench_all_candidates_20261003
for arm in rho-neg rho-plain rho-frozen bsgs bsgs-il bsgs-neg kangaroo ic-mitm ic-subtract; do
  fit "$A" prime "$arm" audit-prime.json "prime-$arm"
done
for arm in rho-strong rho-frob rho-neg rho-plain rho-frozen bsgs bsgs-il bsgs-neg kangaroo ic-frob-m2 ic-frob-m3 ic-subtract; do
  fit "$A" koblitz "$arm" audit-koblitz.json "koblitz-$arm"
done

# Calibration (research/ecbench_calibration_20261002): six sizes per family.
# The Koblitz session spans tiers (degrees 17, 19, 23 are toy; 37, 41, 43 are
# medium), so it is fitted once per tier, three sizes each: constant bounds.
# The prime session (16 to 26 bits) is toy throughout.
C=ecbench_calibration_20261002
for arm in rho-frob rho-neg rho-plain bsgs-il bsgs-neg kangaroo; do
  fit "$C" koblitz "$arm" audit-koblitz.json "calibration-koblitz-medium-$arm" medium
  fit "$C" koblitz "$arm" audit-koblitz.json "calibration-koblitz-toy-$arm" toy
done
for arm in rho-neg rho-plain rho-frozen bsgs bsgs-il bsgs-neg kangaroo; do
  fit "$C" prime "$arm" audit-prime.json "calibration-prime-$arm"
done

# Pair claw (research/ecbench_pair_claw_20261003): six Koblitz degrees 41 to
# 61, all medium, the first exponent-level medium bounds.  The reference, the
# two claw shapes and the generic-table comparator; the A/A control is not a
# bound.
P=ecbench_pair_claw_20261003
for arm in rho-strong claw-c1 claw-pr bsgs-neg; do
  fit "$P" koblitz "$arm" audit-koblitz.json "pair-claw-koblitz-$arm"
done

# Challenge verdicts (research/ecbench_bounds_challenges_20261006): the two
# epoch-1 sessions judged again, every run replayed.  A verdict is
# re-derivable, so the verdict file and the bound it wrote come back byte for
# byte (the command reports on stderr); the committed audit receipts carry a
# timestamp and are left alone.
# An inadmissible verdict stops the script (`--exit-code` under `set -e`).
B=research/ecbench_bounds_challenges_20261006
verdict() { # challenge session epoch bound
  "$E" challenge verdict --challenge "docs/bounds/challenges/$1.json" --dir "$B/sessions/$2" \
    --epoch "$3" --replay-all --bounds "$R" --root . --out "$B/verdict-$2.json" \
    --bound-out "$R/$4.json" --exit-code
}
verdict prime-toy-reference prime-toy-reference-e1 1 challenge-prime-toy-reference-e1-bsgs-negation
verdict koblitz-toy-reference koblitz-toy-reference-e1 1 challenge-koblitz-toy-reference-e1-rho-signed-frobenius

"$E" frontier build --bounds "$R" --out docs/bounds/frontier.json --markdown docs/bounds/FRONTIER.md
"$E" bound check --root . --record "$R"/*.json
if ls docs/bounds/challenges/*.json >/dev/null 2>&1; then
  "$E" challenge check --file docs/bounds/challenges/*.json
fi
