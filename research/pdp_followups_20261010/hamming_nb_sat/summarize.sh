#!/bin/sh
# Markdown table of the panel's frozen rows, one line per (file, solver).
# Thin orchestration over jq; the JSONL in results/ is the evidence.
#   sh summarize.sh > results/summary.md
set -e
here=$(cd "$(dirname "$0")" && pwd)
echo "| n | w | m | base size F_w | T_enum | oracle max wall (ms) | solver | cells | found | unsat | timeout / budget | oracle: decomposable | disagree | conflicts on found (min..max) | conflicts on unsat (min..max) | median wall found (ms) | blocked models |"
echo "|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|:--|:--|--:|--:|"
for f in "$here"/results/n*.jsonl; do
  for s in cryptominisat5 cadical; do
    jq -rs --arg s "$s" '
      (map(select(.kind=="header"))[0]) as $h
      | map(select(.kind=="cell" and .solver==$s)) as $c
      | if ($c|length)==0 then empty else
        ($c|map(select(.verdict=="found"))) as $f
        | ($c|map(select(.verdict=="unsat"))) as $u
        | def rng(a): (a|map(.conflicts)|map(select(.!=null))) as $x
            | if ($x|length)==0 then "-" else "\($x|min)..\($x|max)" end;
          def med(a): (a|map(.wall_ms)|sort) as $w
            | if ($w|length)==0 then "-" else ($w[(($w|length)-1)/2|floor]|tostring) end;
        "| \($h.n) | \($h.w) | \($h.m) | \($h.base_points) | \($c[0].enumeration_floor) | \(($c|map(.ref_wall_us)|max)/1000|floor) | \($s) | \($c|length) | \($f|length) | \($u|length) | \($c|map(select(.verdict!="found" and .verdict!="unsat"))|length) | \($c|map(select(.truth_decomposable))|length) | \($c|map(select(.agrees|not))|length) | \(rng($f)) | \(rng($u)) | \(med($f)) | \($c|map(.blocked_models)|add) |"
        end' "$f"
  done
done
