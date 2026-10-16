#!/bin/bash
# Lay out q1/structure.jsonl: per (t, family_a, family_b) the pooled counts over
# primes (sums only), and per (p, t) the supersingular-vs-torus rows.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
fam='def fam: if test("^ordinary-") then "ordinary" else . end;'
{
  printf 't\tfamily_a\tfamily_b\tcells\tsets_a\tequivalent\tfraction\tcells_all\tcells_none\n'
  jq -s -r "$fam"' map(select(.family_a != .family_b) | .fa = (.family_a|fam) | .fb = (.family_b|fam))
    | group_by([.t, .fa, .fb])[] | . as $g
    | [$g[0].t, $g[0].fa, $g[0].fb, ($g|length), ($g|map(.sets_a)|add), ($g|map(.sets_a_equivalent_to_some_b)|add),
       (($g|map(.sets_a_equivalent_to_some_b)|add) / ($g|map(.sets_a)|add)),
       ($g|map(select(.sets_a_equivalent_to_some_b == .sets_a))|length),
       ($g|map(select(.sets_a_equivalent_to_some_b == 0))|length)] | @tsv' structure.jsonl
} > pooled.tsv
{
  printf 'p\tt\tv2\tfamily_a\tsets_a\tequiv_to_torus\tclasses_a\tbrute_checked\tbrute_agree\n'
  jq -r 'select(.family_b == "torus" and .family_a != "torus") | [.p, .t, .v2_p_plus_1, .family_a, .sets_a, .sets_a_equivalent_to_some_b, .classes_a, (.brute_checked // ""), (.brute_agree // "")] | @tsv' structure.jsonl
} > vs_torus.tsv
{
  printf 'p\tt\tfamily_a\tfamily_b\tsets_a\tequivalent\n'
  jq -r 'select(.family_a == "ss-j1728-iso16" and .family_b == "ss-j1728") | [.p, .t, .family_a, .family_b, .sets_a, .sets_a_equivalent_to_some_b] | @tsv' structure.jsonl
} > positive_control.tsv
