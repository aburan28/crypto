#!/bin/sh
# Thin orchestration of the frozen protocol (PROTOCOL.md).  Usage:
#   sh research/joux_vitse_n1_20261015/run.sh ladder|dlp11|dlp31
set -eu
D=research/joux_vitse_n1_20261015/data
BIN=./target/release/examples/joux_vitse_ladder
mkdir -p "$D/ladder" "$D/dlp"
case "$1" in
ladder)
  for p in 11 31 41 61 71 101 131 151 181 191 211 241; do
    for s in 1 2; do
      out="$D/ladder/p${p}_s${s}.jsonl"; rm -f "$out"
      $BIN trials --p $p --seed $s --f4-trials 32 --mitm-trials $((500 * p)) --json "$out" > /dev/null
    done
  done ;;
dlp11)
  for s in 1 2 3; do
    out="$D/dlp/p11_s${s}.jsonl"; rm -f "$out"
    $BIN dlp --p 11 --seed $s --min-ell-bits 10 --max-trials 200000 --json "$out" > /dev/null
  done ;;
dlp31)
  out="$D/dlp/p31_s1.jsonl"; rm -f "$out"
  timeout 36000 $BIN dlp --p 31 --seed 1 --min-ell-bits 16 --max-trials 200000 --json "$out" > /dev/null ;;
esac
