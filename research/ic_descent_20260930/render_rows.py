#!/usr/bin/env python3
"""Rows of ledger §22's table, rendered from analysis.json.

The scoreboard cites and never computes, so its §22 table is printed
from the frozen analysis rather than typed: one HTML row per size, and
the same rows as Markdown for the ledger note.  Rounding is the only
arithmetic here.

    python3 render_rows.py html > rows.html
    python3 render_rows.py md   > rows.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
A = json.loads((HERE / "analysis.json").read_text())


def times(v: float, digits: int = 2, html: bool = True) -> str:
    return f"{v:.{digits}f}{'&times;' if html else '×'}"


def ci(c: dict, digits: int = 2) -> str:
    return f"[{c['lo']:.{digits}f}, {c['hi']:.{digits}f}]"


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "html"
    for r in A["sizes"]:
        sp, own, dsc = r["speedup"], r["speedup_own_unit"], r["descent_speedup_stage_diagnostic"]
        share = r["descent_share_before"] * 100
        if mode == "html":
            curve = f"<code>K_{r['a']}/GF(2<sup>{r['n']}</sup>)</code> &middot; r = 2<sup>{r['log2_r']:.1f}</sup>"
            print(
                f"        <tr><td>{curve}</td>"
                f"<td class=\"n\">{r['s_after']:.3f} <small>was {r['s_before']:.3f}</small></td>"
                f"<td class=\"n\"><b>{times(r['ratio_after'])}</b> <small>was {times(r['ratio_before'])}</small>"
                f"<br><small>{ci(r['ratio_after_ci'])}</small></td>"
                f"<td class=\"n\"><b>{times(sp['geomean'], 3)}</b><br><small>{ci(sp, 3)}</small></td>"
                f"<td class=\"n\">{times(own['geomean'], 3)}</td>"
                f"<td class=\"n\">{times(dsc['geomean'], 2)}<br><small>{ci(dsc, 2)}</small></td>"
                f"<td class=\"n\">{share:.0f}</td></tr>"
            )
        else:
            print(
                f"| `{r['curve']}` | {r['log2_r']:.1f} | {r['s_before']:.3f} → {r['s_after']:.3f} "
                f"| {times(r['ratio_before'], 2, False)} → **{times(r['ratio_after'], 2, False)}** "
                f"{ci(r['ratio_after_ci'])} "
                f"| **{times(sp['geomean'], 3, False)}** {ci(sp, 3)} "
                f"| {times(own['geomean'], 3, False)} "
                f"| {times(dsc['geomean'], 2, False)} {ci(dsc, 2)} "
                f"| {share:.0f}% |"
            )


if __name__ == "__main__":
    main()
