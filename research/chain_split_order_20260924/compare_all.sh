#!/usr/bin/env bash
# The registered comparisons of RESEARCH_CHAIN_SPLIT_ORDER.md §2–§3, on the
# runs run.sh wrote.  Same-tree pairs (the degree-drop rule D against its
# `complete` twin) go through the frozen suite's compare.py; every pair that
# changes the order O or linear elimination L goes through the cross-tree
# script.  A refusal is saved like an acceptance: its log is the evidence.
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="${1:-$here}"
same=research/groebner_stage_20260915/compare.py
cross=research/inherited_f4_20260922/compare_cross_tree.py
mkdir -p "$root/comparisons"
for suite in frozen chain chain-holdout; do
  s="$root/$suite"
  for pair in "D reference" "LD L" "OD O" "OLD candidate"; do
    read -r cand base <<<"$pair"
    tag="$suite.same_tree.$cand-vs-$base"
    python3 "$same" "$s/$base" "$s/$cand" --output "$root/comparisons/$tag.json" \
      > "$root/comparisons/$tag.log" 2>&1
    echo "$tag: exit $?"
  done
  for cand in candidate L O LD OD OLD; do
    tag="$suite.cross_tree.$cand-vs-reference"
    python3 "$cross" "$s/reference" "$s/$cand" --output "$root/comparisons/$tag.json" \
      > "$root/comparisons/$tag.log" 2>&1
    echo "$tag: exit $?"
  done
done
