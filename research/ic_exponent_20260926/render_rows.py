#!/usr/bin/env python3
"""Rows of ledger §20's tables, rendered from analysis.json.

The scoreboard cites and never computes, so its §20 table is printed
from the frozen analysis rather than typed: one HTML row per size, and
the same rows as Markdown for the ledger note (which alone carries the
cold column; the page gives the cold range in its legend).  Rounding
is the only arithmetic here.

    python3 render_rows.py html > rows.html
    python3 render_rows.py md   > rows.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
A = json.loads((HERE / "analysis.json").read_text())


def curve_html(row: dict) -> str:
    a, n = row["a"], row["n"]
    return f"<code>K_{a}/GF(2<sup>{n}</sup>)</code> &middot; r = 2<sup>{row['log2_r']:.1f}</sup>"


def x(v: float, digits: int = 2) -> str:
    return f"{v:.{digits}f}&times;" if v < 100 else f"{v:,.0f}&times;"


def xm(v: float, digits: int = 2) -> str:
    return f"{v:.{digits}f}×" if v < 100 else f"{v:,.0f}×"


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "html"
    rows = [r for r in A["sizes"] if "ratio" in r]
    for r in rows:
        ci = r["ratio_ci"]
        ch = r["chosen"]
        pts = r["sets"][0]["points"]
        sh = r["shares"]
        cells = {
            "recipe": f"{ch['columns']} &middot; {ch['descent_summands']} &middot; {pts:,}",
            "s_ic": f"{r['s_ic']:.3f}",
            "s_rho": f"{r['s_rho']:.4f}",
            "ratio": f"<b>{x(r['ratio'])}</b><br><small>[{ci['lo']:.2f}, {ci['hi']:.2f}]</small>",
            "law": x(r["law"]),
            "model": x(r["model"]),
            "work": x(r["work_ratio_stage_diagnostic"]),
            "shares": f"{sh['work'] * 100:.0f} / {sh['constructions'] * 100:.0f} / {sh['verification'] * 100:.0f}",
            "cold": x(r["cold_ratio_M1"], 1) if "cold_ratio_M1" in r else "&mdash;",
        }
        if mode == "html":
            print(
                f"        <tr><td>{curve_html(r)}</td><td class=\"n\">{cells['recipe']}</td>"
                f"<td class=\"n\">{cells['s_ic']}</td><td class=\"n\">{cells['s_rho']}</td>"
                f"<td class=\"n\">{cells['ratio']}</td><td class=\"n\">{cells['law']}</td>"
                f"<td class=\"n\">{cells['model']}</td><td class=\"n\">{cells['work']}</td>"
                f"<td class=\"n\">{cells['shares']}</td></tr>"
            )
        else:
            print(
                f"| `{r['curve']}` | {r['log2_r']:.1f} | {ch['columns']}, {ch['descent_summands']}, {pts:,} "
                f"| {r['s_ic']:.3f} | {r['s_rho']:.4f} | **{xm(r['ratio'])}** [{ci['lo']:.2f}, {ci['hi']:.2f}] "
                f"| {xm(r['law'])} | {xm(r['model'])} | {xm(r['work_ratio_stage_diagnostic'])} "
                f"| {sh['work']:.0%} / {sh['constructions']:.0%} / {sh['verification']:.0%} "
                f"| {xm(r['cold_ratio_M1'], 1) if 'cold_ratio_M1' in r else '—'} |"
            )


if __name__ == "__main__":
    main()
