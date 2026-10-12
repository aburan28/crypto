#!/usr/bin/env python3
"""Ledger §23's table rows, from analysis.json: the note's Markdown rows
and the page's HTML rows, so that neither is typed by hand.

    python3 render_rows.py [analysis.json]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def curve_md(s: dict) -> str:
    return f"`K_{s['a']}/GF(2^{s['n']})`"


def curve_html(s: dict) -> str:
    return f"<code>K_{s['a']}/GF(2<sup>{s['n']}</sup>)</code> &middot; r = 2<sup>{s['log2_r']:.1f}</sup>"


def num(v: float) -> str:
    """Two decimals below 10, one above: enough to tell an interval apart."""
    return f"{v:.2f}" if v < 10 else f"{v:.1f}"


def ci(v: list[float], fmt: str = "") -> str:
    return f"[{num(v[0])}, {num(v[1])}]"


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "analysis.json"
    sizes = json.loads(path.read_text())["sizes"]
    print("| curve | log₂ r | S, index calculus online | S, rho online | online speedup [95%] "
          "| at the canonical step [95%] | cold ratio [95%] | S, setup | break-even | IC online over the BL model |")
    print("|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
    for s in sizes:
        o, m, c, b = (s["online_speedup"], s["online_speedup_rho_model"], s["cold_ratio_ic_over_rho"],
                      s["break_even_targets"])
        print(f"| {curve_md(s)} | {s['log2_r']:.1f} | {s['s_ic_online_mean']:.4f} | {s['s_rho_online_mean']:.3f} "
              f"| **{num(o['mean_ratio'])}×** {ci(o['ci95'])} | {num(m['mean_ratio'])}× {ci(m['ci95'])} "
              f"| {num(c['value'])}× {ci(c['ci95'])} | {s['s_setup']:.2f} | {b['value']:.0f} "
              f"| {s['precomputation_boundary_model']['ic_online_over_bl']:.1f}× |")
    print()
    for s in sizes:
        o, m, c, b = (s["online_speedup"], s["online_speedup_rho_model"], s["cold_ratio_ic_over_rho"],
                      s["break_even_targets"])
        aa = [a["online_speedup_R2_over_R1"] for a in s["aa"]]
        print(f"        <tr><td>{curve_html(s)}</td>"
              f"<td class=\"n\">{s['s_ic_online_mean']:.4f}</td>"
              f"<td class=\"n\">{s['s_rho_online_mean']:.3f}</td>"
              f"<td class=\"n\"><b>{num(o['mean_ratio'])}&times;</b><br><small>{ci(o['ci95'])}</small></td>"
              f"<td class=\"n\">{num(m['mean_ratio'])}&times;<br><small>{ci(m['ci95'])}</small></td>"
              f"<td class=\"n\">{num(c['value'])}&times;<br><small>{ci(c['ci95'])}</small></td>"
              f"<td class=\"n\">{s['s_setup']:.2f}</td>"
              f"<td class=\"n\">{b['value']:.0f}</td>"
              f"<td class=\"n\">{s['precomputation_boundary_model']['ic_online_over_bl']:.1f}&times;</td>"
              f"<td class=\"n\">{min(aa):.2f}&ndash;{max(aa):.2f}</td></tr>")


if __name__ == "__main__":
    main()
