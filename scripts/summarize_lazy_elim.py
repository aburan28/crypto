#!/usr/bin/env python3
"""Section 11.14's table: the demand-driven elimination against the full one,
paired cell by cell on section 11.5's protocol.

`research/notes/index-calculus/RESEARCH_RESIDUAL_WALKS.md` §11.12 registered
the border basis with a falsifier; §11.13 bounded it from the solver's own
counters; §11.14 writes it (`GAUDRY_LAZY_ELIM`) and measures it.  This reads
the two protocol runs and prints the one table.

**Why paired, and why the pairing is checked rather than assumed.**  Both runs
use `--cross-check` and the same seeds, so the residual stream is the same by
construction -- *provided* each solve consumes the random stream identically.
The lazy path draws from it only in root finding, after `M_{e1}` is built, and
a fallback that reached root finding and then gave up would draw twice and
desynchronise every later residual.  So each cell asserts that residuals,
decompositions, independent relations and the charpoly, eigenvector and root
counters -- none of which the elimination touches -- are identical in the two
arms before any cost is compared.

**What is charged.**  `oracle_fp_muls / residuals` is `C3`, everything the
solve spends per residual, including the lazy attempts that fell back and the
full elimination that then ran.  `total_ops` is the whole method, cold, in
section 11's unit, and the speedup is AGENTS.md section 8's
`baseline_total_operations / candidate_total_operations`.

Usage:
    python3 scripts/summarize_lazy_elim.py BASELINE.json LAZY.json
"""

from __future__ import annotations

import argparse
import json

#: Counters the elimination does not touch.  Equal in both arms, or the
#: residual streams diverged and nothing below is a paired comparison.
UNTOUCHED = ("residuals", "decompositions", "relations_independent", "correct")
UNTOUCHED_SOLVE = ("solves", "charpoly_muls", "eigenvector_muls", "roots_muls", "unsolved")


def cells(path: str) -> list:
    return [(r["gaudry"], r["rho"]) for r in json.load(open(path))]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("baseline")
    ap.add_argument("lazy")
    args = ap.parse_args()
    base, lazy = cells(args.baseline), cells(args.lazy)
    assert len(base) == len(lazy), (len(base), len(lazy))

    print(
        "| p | log2 n | base C3 | lazy C3 | C3 ratio | elim/res base | elim/res lazy "
        "| lazy solved | fallbacks | needless | M_e1 mismatches | xcheck | speedup | S/rho |"
    )
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    tot_b = tot_l = 0.0
    tot_solved = tot_fb = tot_needless = tot_mm = tot_x = tot_xc = 0
    for (gb, rb), (gl, rl) in zip(base, lazy):
        assert gb["p"] == gl["p"] and gb["base"] == gl["base"], (gb["p"], gl["p"])
        for k in UNTOUCHED:
            assert gb[k] == gl[k], (gb["p"], k, gb[k], gl[k])
        sb, sl = gb["solve_stats"], gl["solve_stats"]
        for k in UNTOUCHED_SOLVE:
            assert sb[k] == sl[k], (gb["p"], k, sb[k], sl[k])
        assert gb["correct"] is True and rb["correct"] is True
        c3b = gb["oracle_fp_muls"] / gb["residuals"]
        c3l = gl["oracle_fp_muls"] / gl["residuals"]
        speed = gb["total_ops"] / gl["total_ops"]
        s_rho = gl["s"] / rl["s"]
        tot_b += gb["total_ops"]
        tot_l += gl["total_ops"]
        tot_solved += sl["lazy_solves"]
        tot_fb += sl["lazy_fallbacks"]
        tot_needless += sl.get("lazy_needless_fallbacks", 0)
        tot_mm += sl["lazy_verify_mismatches"]
        tot_x += gl["cross_check_mismatches"] + gb["cross_check_mismatches"]
        tot_xc += gl["cross_checked"]
        print(
            f"| {gb['p']} | {gb['bits']:.1f} | {c3b:,.0f} | {c3l:,.0f} | {c3l / c3b:.4f} "
            f"| {sb['echelon_muls'] / sb['solves']:,.0f} | {sl['echelon_muls'] / sl['solves']:,.0f} "
            f"| {sl['lazy_solves']} / {sl['lazy_verified']} verified | {sl['lazy_fallbacks']} "
            f"| {sl.get('lazy_needless_fallbacks', 0)} | {sl['lazy_verify_mismatches']} "
            f"| {gl['cross_check_mismatches']} / {gl['cross_checked']} | {speed:.4f} | {s_rho:,.0f} |"
        )
    print()
    print(f"pooled speedup (sum of base total_ops / sum of lazy): {tot_b / tot_l:.5f}")
    print(
        f"lazy solves {tot_solved}, fallbacks {tot_fb} "
        f"({100 * tot_fb / max(tot_solved + tot_fb, 1):.2f}% of lazy attempts), "
        f"needless fallbacks {tot_needless}, M_e1 mismatches {tot_mm}, "
        f"cross-check mismatches {tot_x} over {tot_xc} residuals per arm"
    )


if __name__ == "__main__":
    main()
