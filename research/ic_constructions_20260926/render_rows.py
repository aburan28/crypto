#!/usr/bin/env python3
"""Rows of ledger §21's table, rendered from analysis.json.

The scoreboard cites and never computes, so its §21 table is printed
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
    sign = "&times;" if html else "×"
    return f"{v:.{digits}f}{sign}"


def ci(c: dict, digits: int = 2) -> str:
    return f"[{c['lo']:.{digits}f}, {c['hi']:.{digits}f}]"


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "html"
    for r in A["sizes"]:
        sp, cs = r["speedup"], r["constructions_speedup_stage_diagnostic"]
        share_b = r["constructions_share"]["baseline"] * 100
        share_c = r["constructions_share"]["candidate"] * 100
        if mode == "html":
            curve = f"<code>K_{r['a']}/GF(2<sup>{r['n']}</sup>)</code> &middot; r = 2<sup>{r['log2_r']:.1f}</sup>"
            print(
                f"        <tr><td>{curve}</td>"
                f"<td class=\"n\">{r['s_after']:.3f} <small>was {r['s_before']:.3f}</small></td>"
                f"<td class=\"n\"><b>{times(r['ratio_after'])}</b> <small>was {times(r['ratio_before'])}</small>"
                f"<br><small>{ci(r['ratio_after_ci'])}</small></td>"
                f"<td class=\"n\"><b>{times(sp['geomean'], 3)}</b><br><small>{ci(sp, 3)}</small></td>"
                f"<td class=\"n\">{times(cs['geomean'], 1)}<br><small>{ci(cs, 1)}</small></td>"
                f"<td class=\"n\">{share_c:.1f} <small>was {share_b:.1f}</small></td></tr>"
            )
        else:
            print(
                f"| `{r['curve']}` | {r['log2_r']:.1f} | {r['s_before']:.3f} → {r['s_after']:.3f} "
                f"| {times(r['ratio_before'], 2, False)} → **{times(r['ratio_after'], 2, False)}** "
                f"{ci(r['ratio_after_ci'])} "
                f"| **{times(sp['geomean'], 3, False)}** {ci(sp, 3)} "
                f"| {times(cs['geomean'], 1, False)} {ci(cs, 1)} "
                f"| {share_b:.1f}% → {share_c:.1f}% |"
            )


if __name__ == "__main__":
    main()
