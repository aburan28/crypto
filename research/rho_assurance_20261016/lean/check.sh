#!/usr/bin/env bash
# Build the RhoFormal Lean proofs and audit their axioms.
# Installs elan if missing. Fails on any build error, on any sorry/admit/
# native_decide in the sources, and on any axiom outside the standard three.
set -euo pipefail
cd "$(dirname "$0")"

if ! command -v lake >/dev/null 2>&1; then
  if [ ! -x "$HOME/.elan/bin/lake" ]; then
    curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh \
      | sh -s -- -y --default-toolchain none
  fi
  export PATH="$HOME/.elan/bin:$PATH"
fi

if grep -RnwE 'sorry|admit|native_decide' RhoFormal RhoFormal.lean; then
  echo "check.sh: forbidden tactic in sources" >&2; exit 1
fi
if grep -RnE '^\s*axiom\s' RhoFormal RhoFormal.lean; then
  echo "check.sh: custom axiom declared" >&2; exit 1
fi

lake build 2>&1 | tee build.log
build_status=${PIPESTATUS[0]}
[ "$build_status" -eq 0 ] || { echo "check.sh: lake build failed" >&2; exit 1; }
if grep -q "declaration uses 'sorry'" build.log; then
  echo "check.sh: sorry in build" >&2; exit 1
fi

lake env lean RhoFormal/Axioms.lean 2>&1 | tee axioms.log
[ "${PIPESTATUS[0]}" -eq 0 ] || { echo "check.sh: axioms check failed" >&2; exit 1; }
if grep -q "sorryAx" axioms.log; then
  echo "check.sh: sorryAx in axioms" >&2; exit 1
fi
# Every axiom listed must be one of propext, Classical.choice, Quot.sound.
bad=$(grep -o "depends on axioms: \[[^]]*\]" axioms.log \
  | sed 's/.*\[\(.*\)\]/\1/' | tr ',' '\n' | sed 's/^ *//; s/ *$//' \
  | grep -vxE 'propext|Classical\.choice|Quot\.sound' || true)
if [ -n "$bad" ]; then
  echo "check.sh: non-standard axioms: $bad" >&2; exit 1
fi
expected=$(grep -c '^#print axioms' RhoFormal/Axioms.lean)
got=$(grep -cE "depends on axioms|does not depend on any axioms" axioms.log)
[ "$expected" -eq "$got" ] || { echo "check.sh: expected $expected axiom reports, got $got" >&2; exit 1; }
echo "check.sh: OK ($got theorems audited)"
