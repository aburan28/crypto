#!/usr/bin/env python3
"""Score the surplus control by PREREGISTRATION.md's rules, from committed runs only.

    python3 research/dreg_surplus_control_20260925/score.py

Reads this study's runs/cell-*.jsonl and, for context only, the committed
ell = 4 and ell = 5 cells of the ladder and the grid.  A cell with fewer than
three values is not testable; nothing here infers or fills in a missing cell.
"""
import json
import statistics
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCES = [HERE / "runs",
           HERE.parent / "dreg_ell_grid_20260925" / "runs",
           HERE.parent / "dreg_fixed_surplus_20260923" / "runs"]


def load():
    draws = defaultdict(list)
    for src in SOURCES:
        for f in sorted(src.glob("cell-*.jsonl")):
            for line in f.read_text().splitlines():
                if line.strip() and '"control"' not in line:
                    r = json.loads(line)
                    draws[r["cell"]].append(r)
    return draws


def values(rows):
    """(value, exact?) per unsatisfiable draw; caps-hit draws are excluded."""
    out = []
    for r in rows:
        o = r["outcome"]
        if o["kind"] == "resolved":
            out.append((o["degree"], True))
        elif o["kind"] == "at_least":
            out.append((o["degree"], False))
    return out


def verdict(vs, above, level, below):
    """Median against 6: >= 7 (a bound counts), == 6 exact, < 6 exact."""
    if len(vs) < 3:
        return "not testable (fewer than 3 values)"
    m = statistics.median(v for v, _ in vs)
    if m >= 7:
        return above
    if m == 6 and all(e for v, e in vs if v == 6):
        return level
    if m < 6 and all(e for v, e in vs if v == m):
        return below
    return "not testable (a bound decides it)"


def fmt(vs):
    return " ".join(f"{v}" if e else f"≥{v}" for v, e in vs) or "—"


def main():
    d = load()
    q3 = verdict(values(d["n7l4"]), "rises at S = −5", "no rise at S = −5", "falls at S = −5")
    q4 = verdict(values(d["n11l5"]), "replicates", "does not replicate", "falls")
    print("| cell | ℓ | N | S | values | role |")
    print("|---|--:|--:|--:|---|---|")
    for cell, role in [("n7l4", "**Q3**"), ("n8l4", "secondary"), ("n11l4", "committed"), ("n13l4", "committed"),
                       ("n10l5", "committed (grid Q2)"), ("n11l5", "**Q4**")]:
        rows = d.get(cell, [])
        head = rows[0] if rows else None
        meta = f"{head['ell']} | {head['n_vars']} | {head['surplus']:+d}" if head else "— | — | —"
        print(f"| {cell} | {meta} | {fmt(values(rows))} | {role} |")
    print(f"\nQ3, (7, 4) at S = −5, by the registered rule: {q3}")
    print(f"Q4, (11, 5) at S = −4, by the registered rule: {q4}")
    if q3 == "no rise at S = −5" and q4 == "replicates":
        joint = "the rise belongs to ℓ at these sizes"
    elif q3 == "rises at S = −5":
        joint = "confounded: the equation count is a live explanation"
    elif q4 == "does not replicate":
        joint = "the (10, 5) rise is specific to S = −5"
    else:
        joint = "no registered joint reading applies"
    print(f"Joint reading, as registered: {joint}")
    print("\nSecondary (not decided on)")
    for cell in ("n7l4", "n8l4", "n11l5"):
        rows = d.get(cell, [])
        uns = [r for r in rows if r["outcome"]["kind"] != "satisfiable"]
        sat = sum(1 for r in rows if r["outcome"]["kind"] == "satisfiable")
        if not uns:
            print(f"- {cell}: no unsatisfiable draw committed")
            continue
        secs = [r["secs"] for r in uns]
        drawn = f"satisfiable {sat}/{len(rows)}" if "unsat_index" not in rows[0] else \
            f"unsatisfiable at draw indices {[r['draw'] for r in uns]}"
        print(f"- {cell}: FFD {[r['ffd'] for r in uns]} {drawn} secs {min(secs):.1f}–{max(secs):.1f}")


if __name__ == "__main__":
    main()
