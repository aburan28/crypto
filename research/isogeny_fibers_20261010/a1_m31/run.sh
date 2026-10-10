#!/usr/bin/env bash
# Experiment A1 at m = 31 (RESEARCH_ISOGENY_FIBERS_ECDLP.md §7.A1) on the
# ledger cell icv1-f2m31-tm90707-c95f16f5 (icv1-f2m31-tm90707-c95f16f5), with the
# same seed, repeats and rho reference as research/koblitz_symmetrised_e2e_20260927.
#
# Arms (all on the divisor 0;1;2 subspace, dimension 11, unless noted):
#   x        koblitz-orbit        : the ledger's Frobenius-invariant x-frame base
#   u        koblitz-symmetrised  : F_u = V^{-1}(F_x), the +T-closed base (the
#                                   quotient-isogeny base of §4.4 / H1)
#   tz       koblitz-trace-zero   : divisor 1;2 (dimension 10), all points in [2]E
#   *-proj   the same base with columns merged by cofactor projection
#            (ColumnFold::ProjectedSignedFrobeniusOrbit): P and P + T share a column
# Oracles: mitm-frobenius at m = 2 and m = 3 (the combinatorial arm, exact and
# complete, so yield and column count are measured without a solver constant);
# symmetrised F4 at m = 3 under a wall cap (feasibility only).
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
OUT="$HERE/results"
IC="$ROOT/target/release/ic"
SEED=123212651130
REPEATS=3
RHO_RUNS=8
WALL="${WALL:-2h}"
mkdir -p "$OUT"
N=31; A=0; D="0;1;2"; DTZ="1;2"

declare -A ARM=(
  [mitm-x-m2]="--factor-base koblitz-orbit:divisor=$D --oracle mitm-frobenius:m=2"
  [mitm-x-m2-proj]="--factor-base koblitz-orbit:divisor=$D,projected=1 --oracle mitm-frobenius:m=2"
  [mitm-u-m2]="--factor-base koblitz-symmetrised:divisor=$D --oracle mitm-frobenius:m=2"
  [mitm-u-m2-proj]="--factor-base koblitz-symmetrised:divisor=$D,projected=1 --oracle mitm-frobenius:m=2"
  [mitm-x-m3]="--factor-base koblitz-orbit:divisor=$D --oracle mitm-frobenius:m=3"
  [mitm-x-m3-proj]="--factor-base koblitz-orbit:divisor=$D,projected=1 --oracle mitm-frobenius:m=3"
  [mitm-u-m3]="--factor-base koblitz-symmetrised:divisor=$D --oracle mitm-frobenius:m=3"
  [mitm-u-m3-proj]="--factor-base koblitz-symmetrised:divisor=$D,projected=1 --oracle mitm-frobenius:m=3"
  [mitm-tz-m3]="--factor-base koblitz-trace-zero:divisor=$DTZ --oracle mitm-frobenius:m=3"
  [mitm-tz-m3-proj]="--factor-base koblitz-trace-zero:divisor=$DTZ,projected=1 --oracle mitm-frobenius:m=3"
  [sym-m2-proj]="--factor-base koblitz-symmetrised:divisor=$D,projected=1 --oracle symmetrised:m=2"
  [sym-m3-proj]="--factor-base koblitz-symmetrised:divisor=$D,projected=1 --oracle symmetrised:m=3"
)
ORDER="mitm-x-m2 mitm-x-m2-proj mitm-u-m2 mitm-u-m2-proj mitm-x-m3 mitm-x-m3-proj mitm-u-m3 mitm-u-m3-proj mitm-tz-m3 mitm-tz-m3-proj sym-m2-proj sym-m3-proj"

one() {
  local arm="$1"
  local flags="${ARM[$arm]}"
  local report="$OUT/k0_31__${arm}.json"
  local fb="$OUT/k0_31__${arm}.fb.json"
  local log="$OUT/k0_31__${arm}.log"
  if [ -s "$report" ]; then echo "skip $arm (done)"; return; fi
  rm -f "$fb"
  echo "run  $arm"
  {
    echo "# $(date -u +%FT%TZ) ic bench --koblitz-degree $N --koblitz-a $A $flags --repeats $REPEATS --rho-runs $RHO_RUNS --seed $SEED"
    # shellcheck disable=SC2086
    timeout "$WALL" "$IC" --json --out "$report" bench \
      --koblitz-degree "$N" --koblitz-a "$A" $flags \
      --repeats "$REPEATS" --rho-runs "$RHO_RUNS" --seed "$SEED" \
      --factor-base-out "$fb" >/dev/null
    echo "# exit $? at $(date -u +%FT%TZ)"
  } >"$log" 2>&1
}

for arm in ${*:-$ORDER}; do one "$arm"; done

{
  echo "date: $(date -u +%FT%TZ)"
  echo "host: $(uname -srmo)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //') x$(nproc)"
  echo "rustc: $(rustc --version)"
  echo "commit: $(cd "$ROOT" && git rev-parse HEAD)"
  echo "ic sha256: $(sha256sum "$IC" | cut -d' ' -f1)"
} >"$OUT/host.txt"
