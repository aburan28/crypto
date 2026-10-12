#!/usr/bin/env bash
# Round 4, lead 16.1 (RESEARCH_ISOGENY_FIBERS_ECDLP.md §16.1): the feasibility
# crossover in m.  Two parts, both native Rust:
#
#  A. The per-attempt panel (examples/r4_16_1_m_crossover.rs): plain x-frame
#     and kernel-saturated u-frame decomposition arms at m = 2, 3, 4 on one
#     shared V ∋ 1 per (n, m) and one shared target list per curve, every
#     verdict gated against exact enumeration, one counted cap
#     (node_budget) for every solve, cost in the engine's word XORs.
#  B. End-to-end anchors through `ic bench` (the repository's index-calculus
#     engine, planted logarithm recovered and verified, matched counted rho)
#     on the Frobenius-invariant bases the earlier rounds used: the best
#     feasible m = 2 arms and the m = 3 arms under the same node budget and a
#     trial cap.
#
#   bash research/isogeny_fibers_20261010/round4/16_1_m_crossover/run.sh          # A and B
#   bash research/isogeny_fibers_20261010/round4/16_1_m_crossover/run.sh panel    # A only
#   bash research/isogeny_fibers_20261010/round4/16_1_m_crossover/run.sh anchors  # B only
#
# Resumable: a run whose output exists is skipped.  Single-threaded, serial.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../../../.." && pwd)"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-/tmp/claude-0/-home-user/717e4727-b5a7-5b25-924a-7b31edb43331/scratchpad/r4target}"
BIN="$CARGO_TARGET_DIR/release"
OUT="$HERE/results"
mkdir -p "$OUT"
BUDGET="${BUDGET:-4096}"
SEED=123212651130
WALL="${WALL:-10m}"

panel() {
  (cd "$ROOT" && cargo build --release --example r4_16_1_m_crossover >/dev/null 2>&1)
  local tag="$1"; shift
  local json="$OUT/panel_${tag}.json"
  if [ -s "$json" ] && grep -q '# exit 0' "$OUT/panel_${tag}.log" 2>/dev/null; then echo "skip panel $tag"; return; fi
  echo "run  panel $tag"
  {
    echo "# $(date -u +%FT%TZ) r4_16_1_m_crossover $*"
    "$BIN/examples/r4_16_1_m_crossover" "$@" --json "$json"
    echo "# exit $? at $(date -u +%FT%TZ)"
  } >"$OUT/panel_${tag}.log" 2>&1
}

# key -> "n a divisor-or-dim"
declare -A INST=(
  [k1_17]="17 1 0;1"
  [k0_19]="19 0 dim=10"
  [k0_23]="23 0 0;1"
)
# arm -> "repeats max_trials flags" (D = divisor)
declare -A ARM=(
  [mitm-x-m2]="3 2000000 --factor-base koblitz-orbit:divisor=D --oracle mitm-frobenius:m=2"
  [x-m2]="3 2000000 --factor-base koblitz-orbit:divisor=D --oracle descent-algebraic:m=2 --solver inherited-f4:node_budget=$BUDGET"
  [sym-m2]="3 2000000 --factor-base koblitz-symmetrised:divisor=D --oracle symmetrised:m=2,node_budget=$BUDGET"
  [x-m3]="1 300 --factor-base koblitz-orbit:divisor=D --oracle descent-algebraic:m=3 --solver inherited-f4:node_budget=$BUDGET"
  [sym-m3]="1 300 --factor-base koblitz-symmetrised:divisor=D --oracle symmetrised:m=3,node_budget=$BUDGET"
  # n = 19 has no Frobenius-invariant V ∋ 1 besides F_2: plain subspace base.
  [mitm-s-m2]="3 2000000 --factor-base binary-subspace:dimension=10 --oracle mitm:m=2"
  [xs-m2]="3 2000000 --factor-base binary-subspace:dimension=10 --oracle descent-algebraic:m=2 --solver inherited-f4:node_budget=$BUDGET"
)
declare -A PLAN=(
  [k1_17]="mitm-x-m2 x-m2 sym-m2 sym-m3 x-m3"
  [k0_19]="mitm-s-m2 xs-m2"
  [k0_23]="mitm-x-m2 x-m2 sym-m2 sym-m3 x-m3"
)

anchor() {
  local key="$1" arm="$2"
  read -r n a d <<<"${INST[$key]}"
  read -r reps trials flags <<<"${ARM[$arm]}"
  flags="${flags//D/$d}"
  local report="$OUT/${key}__${arm}.json" log="$OUT/${key}__${arm}.log"
  if [ -s "$report" ]; then echo "skip $key $arm"; return; fi
  echo "run  $key $arm"
  {
    echo "# $(date -u +%FT%TZ) ic bench --koblitz-degree $n --koblitz-a $a $flags --repeats $reps --max-trials $trials --rho-runs 8 --seed $SEED"
    # shellcheck disable=SC2086
    timeout "$WALL" "$BIN/ic" --json --out "$report" bench \
      --koblitz-degree "$n" --koblitz-a "$a" $flags \
      --repeats "$reps" --max-trials "$trials" --rho-runs 8 --seed "$SEED" >/dev/null
    echo "# exit $? at $(date -u +%FT%TZ)"
  } >"$log" 2>&1
}

anchors() {
  (cd "$ROOT" && cargo build --release --bin ic >/dev/null 2>&1)
  for key in k1_17 k0_19 k0_23; do
    for arm in ${PLAN[$key]}; do anchor "$key" "$arm"; done
  done
}

what="${1:-all}"
if [ "$what" = all ] || [ "$what" = panel ]; then
  # One process per curve; the wall guard (contended seconds per cell) only
  # stops a cell whose single solves run for minutes, and is reported.
  panel "n17_b${BUDGET}" --curves 17:1 --ms 2,3,4 --v-spec 17:3:div=0.1 \
    --targets 24 --node-budget "$BUDGET" --cell-wall-s 300
  panel "n19_b${BUDGET}" --curves 19:0 --ms 2,3,4 \
    --targets 24 --node-budget "$BUDGET" --cell-wall-s 300
  # x-direct at m = 3 is skipped at n = 23: one degree-6 solve in 27 unknowns
  # at n = 17 ran 905 s to its node budget, the first attempt at n = 23
  # (27 unknowns) did not return before the 30-minute process limit, and the
  # ic anchor k0_23__x-m3 (36 unknowns) returned nothing in 10 minutes.
  panel "n23_b${BUDGET}" --curves 23:0 --ms 2,3,4 --v-spec 23:3:div=0.1 --skip x-direct@3 \
    --targets 24 --node-budget "$BUDGET" --cell-wall-s 300
fi
if [ "$what" = all ] || [ "$what" = anchors ]; then anchors; fi

{
  echo "date: $(date -u +%FT%TZ)"
  echo "host: $(uname -srmo)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //') x$(nproc) (shared by other agents: contended)"
  echo "rustc: $(rustc --version)"
  echo "commit: $(cd "$ROOT" && git rev-parse HEAD) (worktree has uncommitted round-4 files)"
  [ -x "$BIN/ic" ] && echo "ic sha256: $(sha256sum "$BIN/ic" | cut -d' ' -f1)"
  [ -x "$BIN/examples/r4_16_1_m_crossover" ] && echo "example sha256: $(sha256sum "$BIN/examples/r4_16_1_m_crossover" | cut -d' ' -f1)"
} >"$OUT/host.txt"
