#!/usr/bin/env bash
# External comparison: isogeny_algos against FLINT and PARI/GP on identical inputs
# (examples/ext_compare.rs, external/flint_bench.c, external/pari_bench.gp). Every timed
# command runs through tools/isolated_bench.py on CPU $CPU; two rounds, the order of the
# three tools alternated. Writes results/p15-external/<round>/<tool>.<op>.<p>.jsonl and the
# isolation records to results/p15-external/isolation.jsonl.
#
#   external/run.sh            (from isogeny_algos/; needs gp, libflint-dev, a C compiler)
#   OUT=results/p16-external external/run.sh   (another output directory)
set -euo pipefail
cd "$(dirname "$0")/.."
CPU=${CPU:-3}
# other processes may use up to this many CPUs while a run settles and runs (the agent driving
# the run polls a little); each record still says how much they used and whether it was contended
MAX_OTHER=${MAX_OTHER:-0.5}
# OPS=poly reruns only the polynomial workloads (products, squares, x^p mod g)
OPS=${OPS:-all}
# appended to every isolation-record label (to tell reruns apart in isolation.jsonl)
TAGSFX=${TAGSFX:-}
OUT=${OUT:-results/p15-external}
ISO="python3 ../tools/isolated_bench.py"
mkdir -p "$OUT"
cargo build --release -q --example ext_compare
cc -O2 -o target/flint_bench external/flint_bench.c -lflint -lgmp
OURS=target/release/examples/ext_compare
FLINT=target/flint_bench
{ echo "host: $(uname -srm); $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs); $(nproc) logical CPUs"
  echo "rustc: $(rustc --version)"; echo "cc: $(cc --version | head -1)"
  echo "flint: $(grep -m1 'define FLINT_VERSION ' /usr/include/flint/flint.h || true)"
  echo "pari: $(echo 'print(version())' | gp -q)"; echo "head: $(git rev-parse HEAD)"
  echo "isolation: tools/isolated_bench.py run --cpus $CPU --max-other-cpu $MAX_OTHER"; } > "$OUT/host.txt"

timed() { # label outfile cmd...; retries while the machine is not quiet enough to start
  local label="$1$TAGSFX" out=$2; shift 2
  for try in 1 2 3 4 5 6 7 8 9 10; do
    if $ISO run --wait --cpus "$CPU" --max-other-cpu "$MAX_OTHER" --out "$OUT/isolation.jsonl" \
        --label "$label" -- "$@" > "$out"; then return 0; fi
    sleep 15
  done
  echo "gave up: $label" >&2; return 1
}
gpcmd() { # gp code -> run pari_bench.gp
  printf '%s; read("external/pari_bench.gp")\n' "$1" > "$OUT/.gp_in"
  timed "$2" "$3" gp -q -D parisizemax=4G "$OUT/.gp_in" < /dev/null
}

for round in 1 2; do
  R="$OUT/round$round"; mkdir -p "$R"
  if [ $round = 1 ]; then order="ours flint pari"; else order="pari flint ours"; fi
  for p in p61 p127 p255 p511; do
    for op in poly powx; do
      for tool in $order; do
        case $tool in
          ours)  timed "ours $op $p r$round" "$R/ours.$op.$p.jsonl" $OURS $op $p ;;
          flint) timed "flint $op $p r$round" "$R/flint.$op.$p.jsonl" $FLINT $op $p ;;
          pari)  gpcmd "OP=\"$op\"; TAG=\"$p\"" "pari $op $p r$round" "$R/pari.$op.$p.jsonl" ;;
        esac
      done
    done
  done
  [ "$OPS" = all ] || continue
  for tool in $order; do
    case $tool in
      ours) timed "ours phi p127 r$round" "$R/ours.phi.p127.jsonl" $OURS phi p127 89
            timed "ours sea p61 r$round" "$R/ours.sea.p61.jsonl" $OURS sea p61 20
            timed "ours sea p127 r$round" "$R/ours.sea.p127.jsonl" $OURS sea p127 10 ;;
      pari) gpcmd 'OP="phi"; TAG="p127"; MAXELL=89' "pari phi p127 r$round" "$R/pari.phi.p127.jsonl"
            gpcmd 'OP="sea"; TAG="p61"; CURVES=20' "pari sea p61 r$round" "$R/pari.sea.p61.jsonl"
            gpcmd 'OP="sea"; TAG="p127"; CURVES=10' "pari sea p127 r$round" "$R/pari.sea.p127.jsonl"
            gpcmd 'OP="sea"; TAG="p61"; CURVES=20; NOSEADATA=1' "pari-nosea sea p61 r$round" "$R/pari-nosea.sea.p61.jsonl"
            gpcmd 'OP="sea"; TAG="p127"; CURVES=10; NOSEADATA=1' "pari-nosea sea p127 r$round" "$R/pari-nosea.sea.p127.jsonl" ;;
      flint) ;;
    esac
  done
done
rm -f "$OUT/.gp_in"
echo done
