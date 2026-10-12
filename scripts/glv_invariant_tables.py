#!/usr/bin/env python3
"""Print the tables of RESEARCH_GLV_INVARIANT_FACTOR_BASES.md from the
frozen pilot files.

    python3 scripts/glv_invariant_tables.py experiments/23_glv_invariant_pilot.json \
        experiments/23_glv_invariant_gls_pilot.json

Every number in the note's §5 is printed by this script; the note cites,
it does not compute.
"""
import json
import sys
from collections import defaultdict


def load(paths):
    rows, type_c = [], []
    for p in paths:
        d = json.load(open(p))
        rows.extend(d["rows"])
        type_c.extend(d.get("type_c", []))
    return rows, type_c


def fold_table(rows):
    print(
        "| family | log2 r | h | oracle | signed points | columns fold | columns control | "
        "column ratio | relations fold | relations control | trials fold | trials control | "
        "S fold | S control | S control / S fold | rho S (A = 2) | S fold / rho S | correct |"
    )
    print("|:--|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|")
    for r in sorted(rows, key=lambda r: (r["family"], r["oracle"], r["log2_r"])):
        f, c = r["folded"], r["control"]
        print(
            f"| {r['family']} | {r['log2_r']:.1f} | {r['cofactor']} | {r['oracle']} | "
            f"{f['factor_base']['signed_points']} | {f['factor_base']['columns']} | "
            f"{c['factor_base']['columns']} | {r['column_ratio']:.2f} | "
            f"{f['decomposition']['relations_found']} | {c['decomposition']['relations_found']} | "
            f"{f['decomposition']['targets_tried']} | {c['decomposition']['targets_tried']} | "
            f"{f['s']:.1f} | {c['s']:.1f} | {r['s_ratio']:.2f} | {r['rho_s_mean']:.2f} | "
            f"{f['s'] / r['rho_s_mean']:.0f}× | "
            f"{'yes' if f['verified'] and c['verified'] else 'NO'} |"
        )


def phase_table(rows):
    print(
        "| family | log2 r | oracle | arm | base build | oracle setup | relations | "
        "linear algebra | verify | total GAE | S |"
    )
    print("|:--|--:|:--|:--|--:|--:|--:|--:|--:|--:|--:|")
    for r in sorted(rows, key=lambda r: (r["family"], r["oracle"], r["log2_r"])):
        for arm in ("folded", "control"):
            x = r[arm]
            print(
                f"| {r['family']} | {r['log2_r']:.1f} | {r['oracle']} | {arm} | "
                f"{x['factor_base']['cost']['gae']:.0f} | {x['decomposition']['setup']['gae']:.0f} | "
                f"{x['decomposition']['cost']['gae']:.0f} | {x['linear_algebra']['cost']['gae']:.0f} | "
                f"{x['verify']['gae']:.0f} | {x['total_gae']:.0f} | {x['s']:.1f} |"
            )


def summary_table(rows):
    """Per family and oracle: the range of the column, trial and S ratios."""
    groups = defaultdict(list)
    for r in rows:
        groups[(r["family"], r["oracle"])].append(r)
    print("| family | oracle | instances | column ratio | trial ratio (min–max) | S ratio (min–max) | S fold / rho S (min–max) |")
    print("|:--|:--|--:|--:|--:|--:|--:|")
    for (fam, orc), rs in sorted(groups.items()):
        col = {round(r["column_ratio"], 2) for r in rs}
        tr = [r["trial_ratio"] for r in rs]
        sr = [r["s_ratio"] for r in rs]
        vr = [r["folded"]["s"] / r["rho_s_mean"] for r in rs]
        print(
            f"| {fam} | {orc} | {len(rs)} | {', '.join(str(c) for c in sorted(col))} | "
            f"{min(tr):.2f}–{max(tr):.2f} | {min(sr):.2f}–{max(sr):.2f} | {min(vr):.0f}×–{max(vr):.0f}× |"
        )


def type_c_table(type_c):
    print("| family | log2 r | rational degree-2 maps | ord_r(λ) | images in base | base points | chance fraction | verified |")
    print("|:--|--:|--:|--:|--:|--:|--:|:--|")
    for t in sorted(type_c, key=lambda t: (t["family"], t["r"])):
        o = t.get("overlap")
        if o is None:
            print(f"| {t['family']} | {t['r'].bit_length()} | {t['endomorphisms_found']} | — | — | — | — | {t['verified']} |")
            continue
        import math

        print(
            f"| {t['family']} | {math.log2(t['r']):.1f} | {t['endomorphisms_found']} | {o['eigenvalue_order']} | "
            f"{o['images_in_base']} | {o['base_points']} | {o['chance_fraction']:.5f} | {t['verified']} |"
        )


def main():
    rows, type_c = load(sys.argv[1:])
    print("## Fold against control, every row\n")
    fold_table(rows)
    print("\n## Summary by family and oracle\n")
    summary_table(rows)
    print("\n## Phase costs (group-addition equivalents)\n")
    phase_table(rows)
    if type_c:
        print("\n## Type C: degree-2 CM endomorphisms\n")
        type_c_table(type_c)


if __name__ == "__main__":
    main()
