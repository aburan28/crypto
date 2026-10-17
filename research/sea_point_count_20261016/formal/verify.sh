#!/usr/bin/env bash
# Thin tool orchestration. Arithmetic, measurements and verification are native.
set -euo pipefail
study=$(cd "$(dirname "$0")/.." && pwd)
jar=${TLA_JAR:-${RUNNER_TEMP:-/tmp}/sea-tla2tools-1.8.0.jar}
logs=${FORMAL_LOG_DIR:-${RUNNER_TEMP:-/tmp}/sea-formal-replay}
mkdir -p "$logs"
if [ ! -f "$jar" ]; then
  curl --proto '=https' --tlsv1.2 -fsSL \
    https://github.com/tlaplus/tlaplus/releases/download/v1.8.0/tla2tools.jar -o "$jar"
fi
printf '15d103b7e5617780bcc904c0af12f1ad242e2b345c0301136c5be4588f17edac  %s\n' "$jar" | sha256sum -c -
cd "$study/formal"
check() {
  local model=$1 config=$2 expected=$3 invariant=${4:-}
  local code=0
  timeout 300 java -XX:+UseParallelGC -Xmx2g -cp "$jar" tlc2.TLC \
    -workers 1 -seed 1 -fp 0 -metadir "$logs/states-$config" \
    -config "$config.cfg" "$model.tla" > "$logs/$config.txt" 2>&1 || code=$?
  printf '%s\n' "$code" > "$logs/$config.exit"
  if [ "$code" -ne "$expected" ]; then cat "$logs/$config.txt"; exit 1; fi
  if [ "$expected" -eq 0 ]; then
    grep -q 'Model checking completed. No error has been found.' "$logs/$config.txt"
  else
    grep -q "Invariant $invariant is violated" "$logs/$config.txt"
  fi
  printf '%s: expected exit %s verified\n' "$config" "$expected"
}
check Recombination Fixed 0
check Recombination Wide 0
check Recombination Legacy 12 NontrivialTrueSurvives
check Certification Certification 0
check Certification SampleOnly 12 NontrivialSampleSound
cd "$study/lean"
lake build > "$logs/lean-build.txt" 2>&1
lake env lean SeaProofs.lean > "$logs/lean-axioms.txt" 2>&1
if grep -Eq 'sorryAx|^error:' "$logs/lean-axioms.txt"; then cat "$logs/lean-axioms.txt"; exit 1; fi
printf 'Lean build and axiom inspection passed\n'
