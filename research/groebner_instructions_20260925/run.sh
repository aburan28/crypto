#!/usr/bin/env bash
# Registered runs of research/notes/ecc2k130/RESEARCH_GROEBNER_STAGE_INSTRUCTIONS.md §2.
# Arms A0 (before #690), A1 (#690) and A2 (#712, the default) on one binary.
#   stage/       every rung of the five ladders, instructions inside groebner_decompose
#   determinism/ three rungs per arm, run a second time (gate G2)
#   e2e/         whole logarithms, instructions of the whole ic run process
#   rho/         the counted rho reference on the same seeds, same unit
#   add/         instructions per group addition (the conversion to GAE)
# Nothing is overwritten: a finished job is skipped.
#
#   cargo build --release --example groebner_stage_bench --example ir_calibration --bin ic
#   research/groebner_instructions_20260925/run.sh
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
export BENCH="${BENCH:-$PWD/target/release/examples/groebner_stage_bench}"
export IC="${IC:-$PWD/target/release/ic}"
export CAL="${CAL:-$PWD/target/release/examples/ir_calibration}"
export HERE="$here"
jobs="${JOBS:-4}"

arm_env() { # A0 | A1 | A2
  case $1 in
    A0) echo "KIC_CHAIN_ORDER=layout KIC_LINEAR_ELIM=0 KIC_F4_MULTIPLIERS=occurring KIC_F4_DROP=complete" ;;
    A1) echo "KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_MULTIPLIERS=occurring KIC_F4_DROP=complete" ;;
    A2) echo "KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_MULTIPLIERS=support KIC_F4_DROP=complete" ;;
  esac
}
export -f arm_env

stage_one() { # KIND SUITE RUNG ARM  (KIND = stage | determinism)
  local kind=$1 suite=$2 rung=$3 arm=$4
  local d="$HERE/$kind/$suite/$arm/rung$rung"
  [ -e "$d/summary.json" ] && return 0
  mkdir -p "$d"
  local ladder=(); [ "$suite" != frozen ] && ladder=(--ladder "$suite")
  env $(arm_env "$arm") valgrind --tool=callgrind --toggle-collect='*groebner_decompose*' \
    --callgrind-out-file="$d/callgrind.out" "$BENCH" --label "$arm" "${ladder[@]}" --rung "$rung" \
    --out "$d/stage" > "$d/stdout" 2> "$d/valgrind.log"
  python3 "$HERE/summarise.py" stage "$d" && gzip -9 -f "$d/callgrind.out"
}
export -f stage_one

e2e_one() { # N SEED ARM
  local n=$1 seed=$2 arm=$3
  local d="$HERE/e2e/$arm/K0_2^${n}_seed$seed"
  [ -e "$d/summary.json" ] && return 0
  mkdir -p "$d"
  env $(arm_env "$arm") valgrind --tool=callgrind --callgrind-out-file="$d/callgrind.out" \
    "$IC" run --degree "$n" --curve-a 0 --summands 3 --solver groebner --random-target \
    --seed "$seed" --batch 1 --json > "$d/run.json" 2> "$d/valgrind.log"
  python3 "$HERE/summarise.py" e2e "$d" && gzip -9 -f "$d/callgrind.out"
}
export -f e2e_one

rho_one() { # N SEED
  local n=$1 seed=$2
  local d="$HERE/rho/K0_2^${n}_seed$seed"
  [ -e "$d/summary.json" ] && return 0
  mkdir -p "$d"
  valgrind --tool=callgrind --toggle-collect='*rho_solve*' --callgrind-out-file="$d/callgrind.out" \
    "$CAL" rho 0 "$n" "$seed" > "$d/run.json" 2> "$d/valgrind.log"
  python3 "$HERE/summarise.py" rho "$d" && gzip -9 -f "$d/callgrind.out"
}
export -f rho_one

add_one() { # N
  local n=$1
  local d="$HERE/add/K0_2^$n"
  [ -e "$d/summary.json" ] && return 0
  mkdir -p "$d"
  valgrind --tool=callgrind --toggle-collect='*add_loop*' --callgrind-out-file="$d/callgrind.out" \
    "$CAL" add 0 "$n" 200000 > "$d/run.json" 2> "$d/valgrind.log"
  python3 "$HERE/summarise.py" add "$d" && gzip -9 -f "$d/callgrind.out"
}
export -f add_one

{
  for n in 13 9; do echo "add_one $n"; done
  for arm in A0 A1 A2; do
    for r in 0 1 2 3 4 5; do echo "stage_one stage frozen $r $arm"; done
    for r in 0 1 4; do echo "stage_one stage chain $r $arm"; done
    for r in 6 7 8 9; do echo "stage_one stage chain-holdout $r $arm"; done
    for r in $(seq 0 8); do echo "stage_one stage chain-holdout-2 $r $arm"; done
    for r in $(seq 0 18); do echo "stage_one stage r2-holdout $r $arm"; done
    echo "stage_one determinism frozen 2 $arm"
    echo "stage_one determinism chain-holdout-2 7 $arm"
    echo "stage_one determinism r2-holdout 5 $arm"
    for s in 201 202 203 204 205 206 207 208 209 210 301 302 303 304 305; do echo "e2e_one 13 $s $arm"; done
    for s in 201 202 203 204 205; do echo "e2e_one 9 $s $arm"; done
  done
  for s in 201 202 203 204 205 206 207 208 209 210 301 302 303 304 305; do echo "rho_one 13 $s"; done
  for s in 201 202 203 204 205; do echo "rho_one 9 $s"; done
} | xargs -P "$jobs" -I{} bash -c '{} && echo "done: {}" || echo "FAILED: {}"'
