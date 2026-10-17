#!/usr/bin/env bash
# Per-family external comparison: isogeny_algos (examples/ext_families.rs) against PARI/GP
# (external/pari_families.gp) on identical inputs with identical outputs: kernel -> codomain,
# kernel -> rational maps, l^n chains. Every timed command runs through tools/isolated_bench.py
# on CPU $CPU; two rounds, the tool order reversed in the second.
#
#   OUT=results/p18-families external/families.sh   (from isogeny_algos/; needs gp)
set -euo pipefail
cd "$(dirname "$0")/.."
CPU=${CPU:-3}
MAX_OTHER=${MAX_OTHER:-0.5}
OUT=${OUT:-results/p18-families}
ISO="python3 ../tools/isolated_bench.py"
mkdir -p "$OUT"
cargo build --release -q --example ext_families
OURS=target/release/examples/ext_families
$OURS gen > "$OUT/cases.gp"
{ echo "host: $(uname -srm); $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs); $(nproc) logical CPUs"
  echo "rustc: $(rustc --version)"; echo "pari: $(echo 'print(version())' | gp -q)"
  echo "head: $(git rev-parse HEAD)"
  echo "isolation: tools/isolated_bench.py run --cpus $CPU --max-other-cpu $MAX_OTHER"; } > "$OUT/host.txt"

timed() { # label outfile cmd...; retries while the machine is not quiet enough to start
  local label=$1 out=$2; shift 2
  for try in 1 2 3 4 5 6 7 8 9 10; do
    if $ISO run --wait --cpus "$CPU" --max-other-cpu "$MAX_OTHER" --out "$OUT/isolation.jsonl" \
        --label "$label" -- "$@" > "$out"; then return 0; fi
    sleep 15
  done
  echo "gave up: $label" >&2; return 1
}

for round in 1 2; do
  R="$OUT/round$round"; mkdir -p "$R"
  if [ $round = 1 ]; then order="ours pari"; else order="pari ours"; fi
  for op in codomain maps chain; do
    for tool in $order; do
      case $tool in
        ours) timed "ours $op r$round" "$R/ours.$op.jsonl" $OURS $op ;;
        pari) printf 'OP="%s"; read("%s"); read("external/pari_families.gp")\n' "$op" "$OUT/cases.gp" > "$OUT/.gp_in"
              timed "pari $op r$round" "$R/pari.$op.jsonl" gp -q -D parisizemax=4G "$OUT/.gp_in" < /dev/null ;;
      esac
    done
  done
done
rm -f "$OUT/.gp_in"
echo done
