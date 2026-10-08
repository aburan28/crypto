#!/usr/bin/env bash
# Head-engine rerun of the m = 4 exponent audit: the nine registered cells of
# ../ic_m4_exponent_audit_20260928 (§8), on the head engine (PREREGISTRATION.md here).
#
#   BIN=$WORK/bin/m4_exponent_audit-4ff512f2 \
#     research/ic_m4_head_engine_20260929/run_audit.sh [--print] [RUNS_DIR]
#
# --print only prints the per-cell commands.  Cells are independent processes, one
# thread each; the Semaev/enumeration/null cells run four at a time, then the degree
# cells one at a time (memory).  A cell's CPU-seconds limit (ulimit -t) is machine
# protection: a killed cell keeps every line it wrote, and its unwritten targets are
# censored -- never negative evidence.  Outputs are never overwritten.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
print_only=0
if [ "${1:-}" = "--print" ]; then print_only=1; shift; fi
runs="${1:-$here/runs}"
bin="${BIN:?set BIN to the m4_exponent_audit-4ff512f2 binary that build.sh produced}"
seed=20260928
engine_env="KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete"
# Degree arm only: the dreg_ladder caps, so a degree is built in full or not at all.
degree_env="$engine_env F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000"

# Registered cells (a n): the survey's K_0, K_1 x {9, 11, 13, 15, 17, 19} minus the three
# that do not exist in the tooling (K_0/2^11, K_0/2^17, K_1/2^13; §1).
cells="0 9|0 13|0 15|0 19|1 9|1 11|1 15|1 17|1 19"
degree_cells="0 9|1 9|1 11|0 13"

semaev_cpu() { case $1 in 17) echo 300;; 19) echo 600;; *) echo 120;; esac; }

cmd() { # arm a n targets cpu_limit env extra...
  local arm=$1 a=$2 n=$3 t=$4 cpu=$5 envs=$6; shift 6
  local ell=$(( (n + 2) / 4 ))
  local cell="K${a}n${n}l${ell}"
  local mem=""
  [ "$arm" = degree ] && mem="ulimit -v 10000000; "
  echo "mkdir -p $runs/$arm && ( ulimit -t $cpu; ${mem}exec env -i PATH=\"\$PATH\" $envs $bin --arm $arm --a $a --n $n --targets $t --first 0 --seed $seed --label audit $* --out $runs/$arm/$cell.jsonl ) > /dev/null 2> $runs/$arm/$cell.stderr; echo \$? > $runs/$arm/$cell.exit"
}

parallel_cmds=()
IFS='|' read -ra list <<<"$cells"
for c in "${list[@]}"; do
  read -r a n <<<"$c"
  parallel_cmds+=("$(cmd semaev "$a" "$n" 16 "$(semaev_cpu "$n")" "$engine_env" --node-budget 20000)")
  parallel_cmds+=("$(cmd enumerate "$a" "$n" 16 60 "$engine_env")")
  parallel_cmds+=("$(cmd null "$a" "$n" 8 120 "$engine_env" --node-budget 20000)")
done
serial_cmds=()
IFS='|' read -ra dlist <<<"$degree_cells"
for c in "${dlist[@]}"; do
  read -r a n <<<"$c"
  serial_cmds+=("$(cmd degree "$a" "$n" 16 300 "$degree_env" --node-budget 20000 --d-max 6 --max-unsat 4)")
done

if [ $print_only = 1 ]; then
  printf '%s\n' "${parallel_cmds[@]}" "${serial_cmds[@]}"
  exit 0
fi
[ -x "$bin" ] || { echo "no binary at $bin; run build.sh" >&2; exit 1; }
[ -e "$runs" ] && { echo "$runs exists; never overwritten" >&2; exit 1; }
mkdir -p "$runs"
sha256sum "$bin" > "$runs/binary.sha256"
git -C "$here" rev-parse HEAD > "$runs/branch_head.txt"
date -u +%FT%TZ > "$runs/started_utc.txt"
# Pinning (AGENTS.md section 10): three parallel slots, each pinned to its own CPU
# (1, 2, 3); CPU 0 is left to the system and the driver.  Degree cells run one at
# a time on CPU 3.  The counted units are deterministic; pinning keeps each cell's
# CPU-seconds limit, which decides censoring, free of time-slicing by its siblings.
printf '%s\n' "${parallel_cmds[@]}" | xargs -d '\n' -P 3 --process-slot-var=SLOT -I{} \
  bash -c 'exec taskset -c $((SLOT + 1)) bash -c "$1"' _ {}
printf '%s\n' "${serial_cmds[@]}" | xargs -d '\n' -P 1 -I{} taskset -c 3 bash -c '{}'
date -u +%FT%TZ > "$runs/finished_utc.txt"
python3 "$here/../ic_m4_exponent_audit_20260928/analyze.py" "$runs" | tee "$runs/readout.txt"
