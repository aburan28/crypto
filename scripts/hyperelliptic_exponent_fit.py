#!/usr/bin/env python3
"""Fit the genus exponent law, and correct for a drifting rho reference.

`RESEARCH_HYPERELLIPTIC_IC_RHO.md` round five pre-registers

    S_ic / S_rho  ∝  g! · p^(1 − g/2)

so the slope of `log2(S_ic/S_rho)` against `log2(p)` should be `1 − g/2`.

It also pre-registers a falsification condition that fired: `S_walk`, rho's
walk cost alone, must stay inside `1.0 – 1.5` against the `sqrt(pi/2) =
1.2533` ideal, and at genus 3 it drifts upward with `p` (1.16 → 1.82).  A
reference that gets more expensive as the group grows depresses the measured
ratio, which flatters the algorithm under study -- the wrong direction for a
reference to be wrong in.

So two fits are reported, never one: the raw ratio, and a ratio with the walk
term renormalised to the ideal,

    S_rho_corrected = S_rho - S_walk + sqrt(pi/2)

which charges rho what a correct walk would cost and leaves its (measured,
not modelled) precomputation alone.  If the trend only exists in the raw
column it is an artefact of the reference; if it survives the correction it
is a property of the algorithms.

    python3 scripts/hyperelliptic_exponent_fit.py [LOG]
"""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

IDEAL_WALK = math.sqrt(math.pi / 2)          # 1.25331...
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "experiments/hyperelliptic_exponent_fit.json"

ROW = re.compile(
    r"^\s*(\d+)\s+(baseline|optimised)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+"
    r"([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$")


def parse(text: str):
    genus = None
    rows = []
    for line in text.splitlines():
        m = re.match(r"^---\s*genus\s+(\d+)\s*---", line.strip())
        if m:
            genus = int(m.group(1))
            continue
        m = ROW.match(line)
        if not m or genus is None:
            continue
        p, arm = int(m.group(1)), m.group(2)
        n, fb = int(m.group(3)), int(m.group(4))
        s_ic, s_rho, ratio, s_walk = (float(m.group(7)), float(m.group(8)),
                                      float(m.group(9)), float(m.group(10)))
        rows.append({
            "genus": genus, "p": p, "arm": arm, "N": n, "factor_base": fb,
            "S_ic": s_ic, "S_rho": s_rho, "ratio_raw": ratio,
            "S_walk": s_walk,
        })
    return rows


def correct(row):
    """Renormalise rho's walk to the ideal, leaving its precomputation alone."""
    s_rho_corr = row["S_rho"] - row["S_walk"] + IDEAL_WALK
    row["S_rho_corrected"] = round(s_rho_corr, 4)
    row["ratio_corrected"] = (round(row["S_ic"] / s_rho_corr, 4)
                              if s_rho_corr > 0 else None)
    row["walk_inside_band"] = 1.0 <= row["S_walk"] <= 1.5
    return row


def fit(xs, ys):
    """Least-squares slope and intercept of y against x, plus R^2."""
    n = len(xs)
    if n < 2:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    intercept = my - slope * mx
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return {"slope": round(slope, 4), "intercept": round(intercept, 4),
            "r2": round(r2, 4), "points": n}


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/hyper5.log")
    rows = [correct(r) for r in parse(path.read_text())]
    if not rows:
        print("no rows parsed", file=sys.stderr)
        return 1

    out = {
        "source_log": str(path),
        "ideal_walk": round(IDEAL_WALK, 5),
        "predicted_slope": {g: 1 - g / 2 for g in (2, 3, 4)},
        "rows": rows,
        "fits": {},
    }

    print(f"{'g':>2} {'p':>5} {'N':>10} {'S_ic':>7} {'S_rho':>7} {'raw':>6} "
          f"{'S_walk':>7} {'S_rho*':>7} {'corr':>6} band")
    print("-" * 74)
    for r in rows:
        if r["arm"] != "optimised":
            continue
        print(f"{r['genus']:>2} {r['p']:>5} {r['N']:>10} {r['S_ic']:>7.2f} "
              f"{r['S_rho']:>7.2f} {r['ratio_raw']:>6.2f} {r['S_walk']:>7.2f} "
              f"{r['S_rho_corrected']:>7.2f} {r['ratio_corrected']:>6.2f} "
              f"{'ok' if r['walk_inside_band'] else 'OUT'}")

    print()
    print(f"{'g':>2} {'pred':>6} {'raw slope':>10} {'r2':>6} "
          f"{'corr slope':>11} {'r2':>6} {'n':>3}  verdict")
    print("-" * 68)
    for g in sorted({r["genus"] for r in rows}):
        sel = [r for r in rows if r["genus"] == g and r["arm"] == "optimised"]
        xs = [math.log2(r["p"]) for r in sel]
        raw = fit(xs, [math.log2(r["ratio_raw"]) for r in sel])
        cor = fit(xs, [math.log2(r["ratio_corrected"]) for r in sel])
        pred = 1 - g / 2
        out["fits"][str(g)] = {"predicted_slope": pred, "raw": raw,
                              "corrected": cor,
                              "walk_out_of_band": [r["p"] for r in sel
                                                   if not r["walk_inside_band"]]}
        if raw and cor:
            ok = abs(cor["slope"] - pred) <= 0.25
            print(f"{g:>2} {pred:>6.2f} {raw['slope']:>10.3f} {raw['r2']:>6.2f} "
                  f"{cor['slope']:>11.3f} {cor['r2']:>6.2f} {raw['points']:>3}  "
                  f"{'within +-0.25' if ok else 'OUTSIDE +-0.25'}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
