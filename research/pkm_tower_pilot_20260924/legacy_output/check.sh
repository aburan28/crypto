#!/bin/bash
# The native checker (examples/pkm_tower_check.rs) must print what the legacy
# Python scripts printed, line for line, with the same exit status (note
# section 16.3). The reference is CI's output on 2026-09-30: workflow run
# 36689898077, job 109804351090, ubuntu-24.04 with its Python 3.12, at commit
# a8b1e7ad, whose rows and scripts are those on main. It covers verify.py on
# each round's rows and analyze.py on rounds 1-5, as that workflow ran them.
# From the repository root:
#
#     cargo build --release --example pkm_tower_check
#     research/pkm_tower_pilot_20260924/legacy_output/check.sh
set -u
export LC_ALL=C
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../../.." && pwd)
bin=${BIN:-$root/target/release/examples/pkm_tower_check}
cd "$root" || exit 2
status=0
check() {
  name=$1
  shift
  out=$("$bin" "$@")
  rc=$?
  if [ "$rc" -eq 0 ] && [ "$out" == "$(cat "$here/$name.txt")" ]; then
    echo "$name: identical ($(wc -l < "$here/$name.txt") lines), exit 0"
  else
    echo "$name: DIFFERS (exit $rc)"
    diff <(printf '%s\n' "$out") "$here/$name.txt" | head -40
    status=1
  fi
}
r=research
check verify-pilot verify $r/pkm_tower_pilot_20260924/runs/*.jsonl
check verify-round2 verify $r/pkm_tower_round2_20260925/runs/*.jsonl
check verify-round3 verify $r/pkm_tower_round3_20260925/runs/*.jsonl
check verify-round4 verify $r/pkm_tower_round4_20260926/runs/*.jsonl
check verify-round5 verify $r/pkm_tower_round5_20260926/runs/*.jsonl
check verify-round6 verify $r/pkm_tower_round6_20260930/replay/*.jsonl \
  $r/pkm_tower_round6_20260930/off/*.jsonl
check analyze-rounds1-5 analyze $r/pkm_tower_pilot_20260924/runs/*.jsonl \
  $r/pkm_tower_round2_20260925/runs/*.jsonl $r/pkm_tower_round3_20260925/runs/*.jsonl \
  $r/pkm_tower_round4_20260926/runs/*.jsonl $r/pkm_tower_round5_20260926/runs/*.jsonl
exit $status
