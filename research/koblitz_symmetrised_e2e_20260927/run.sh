#!/usr/bin/env bash
# The exact runs behind results/.  Every arm on an instance shares the
# seed (so the planted logarithms and the rho reference are matched),
# three repeats and eight counted rho runs.  A run whose report already
# exists is skipped, so the script is resumable; a run that fails or
# hits the wall cap leaves its log and no report.
#
#   bash research/koblitz_symmetrised_e2e_20260927/run.sh          # everything
#   bash research/koblitz_symmetrised_e2e_20260927/run.sh k1_17    # one instance
#
# Builds the release `ic` binary if it is missing.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
OUT="$HERE/results"
IC="$ROOT/target/release/ic"
SEED=123212651130
REPEATS=3
RHO_RUNS=8
WALL=3h
# Serial by default: the algebraic arms are priced by measured wall
# (their word-XOR count is partial, so the runner will not price by it),
# and parallel runs would share that price.
JOBS="${JOBS:-1}"
mkdir -p "$OUT"
[ -x "$IC" ] || (cd "$ROOT" && cargo build --release --bin ic)

# instance key  -> degree a divisor
declare -A INST=(
  [k1_17]="17 1 0;1"
  [k0_23]="23 0 0;1"
  [k1_23]="23 1 0;1"
  [k0_31]="31 0 0;1;2"
)
# arm -> extra bench flags (D is replaced by the divisor)
declare -A ARM=(
  [sym]="--factor-base koblitz-symmetrised:divisor=D --oracle symmetrised:m=2"
  [sym4]="--factor-base koblitz-symmetrised:divisor=D --oracle symmetrised:m=2,max_degree=4"
  [x]="--factor-base koblitz-orbit:divisor=D --oracle descent-algebraic:m=2 --solver inherited-f4"
  [mitm-x]="--factor-base koblitz-orbit:divisor=D --oracle mitm-frobenius:m=2"
  [mitm-u]="--factor-base koblitz-symmetrised:divisor=D --oracle mitm-frobenius:m=2"
  [sym-m3]="--factor-base koblitz-symmetrised:divisor=D --oracle symmetrised:m=3"
  [x-m3]="--factor-base koblitz-orbit:divisor=D --oracle descent-algebraic:m=3 --solver inherited-f4"
)
ORDER_ARMS="sym sym4 x mitm-x mitm-u sym-m3 x-m3"
M3_ONLY_AT="k1_17"

one() {
  local key="$1" arm="$2"
  read -r n a d <<<"${INST[$key]}"
  local flags="${ARM[$arm]//D/$d}"
  local report="$OUT/${key}__${arm}.json"
  local fb="$OUT/${key}__${arm}.fb.json"
  local log="$OUT/${key}__${arm}.log"
  if [ -s "$report" ]; then echo "skip $key $arm (done)"; return; fi
  rm -f "$fb"
  echo "run  $key $arm"
  {
    echo "# $(date -u +%FT%TZ) ic bench --koblitz-degree $n --koblitz-a $a $flags --repeats $REPEATS --rho-runs $RHO_RUNS --seed $SEED"
    # shellcheck disable=SC2086
    timeout "$WALL" "$IC" --json --out "$report" bench \
      --koblitz-degree "$n" --koblitz-a "$a" $flags \
      --repeats "$REPEATS" --rho-runs "$RHO_RUNS" --seed "$SEED" \
      --factor-base-out "$fb" >/dev/null
    echo "# exit $? at $(date -u +%FT%TZ)"
  } >"$log" 2>&1
}
export -f one
export OUT IC SEED REPEATS RHO_RUNS WALL
export INST_SER="$(declare -p INST)" ARM_SER="$(declare -p ARM)"

keys="${*:-k1_17 k0_23 k1_23 k0_31}"
jobs=()
late=()
for key in $keys; do
  for arm in $ORDER_ARMS; do
    case "$arm" in
      *-m3) [[ " $M3_ONLY_AT " == *" $key "* ]] || continue; late+=("$key $arm"); continue ;;
      sym4) [ "$key" = k0_31 ] && { late+=("$key $arm"); continue; } ;;
    esac
    jobs+=("$key $arm")
  done
done
jobs+=("${late[@]}")
printf '%s\n' "${jobs[@]}" | xargs -P "$JOBS" -I{} bash -c 'eval "$INST_SER"; eval "$ARM_SER"; one {}'

{
  echo "date: $(date -u +%FT%TZ)"
  echo "host: $(uname -srmo)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //') x$(nproc)"
  echo "rustc: $(rustc --version)"
  echo "commit: $(cd "$ROOT" && git rev-parse HEAD)"
  echo "ic sha256: $(sha256sum "$IC" | cut -d' ' -f1)"
} >"$OUT/host.txt"
