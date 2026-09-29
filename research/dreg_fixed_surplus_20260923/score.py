#!/usr/bin/env python3
"""Score the fixed-surplus ladder by PREREGISTRATION.md's rules, from committed runs only.

    python3 research/dreg_fixed_surplus_20260923/score.py

Reads runs/cell-*.jsonl.  A cell with no committed draws has no values;
nothing here infers, extrapolates or fills in a missing cell.
"""
import json
import statistics
from collections import defaultdict
from pathlib import Path

RUNS = Path(__file__).resolve().parent / "runs"
PAIRS = [  # (label, surplus, small cell, large cell, large d_max), as registered
    ("primary", -2, "n7l3", "n13l5", 6),
    ("", -1, "n5l2", "n11l4", 6),
    ("", 0, "n9l3", "n15l5", 6),
    ("", +1, "n7l2", "n13l4", 6),
]


def load():
    draws, controls = defaultdict(list), {}
    for f in sorted(RUNS.glob("cell-*.jsonl")):
        for line in f.read_text().splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if "control" in r:
                controls[r["cell"]] = r["outcome"]
            else:
                draws[r["cell"]].append(r)
    return draws, controls


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


def main():
    draws, controls = load()
    print("| pair | S | small | values | median | large | values | median | verdict |")
    print("|---|--:|---|---|--:|---|---|--:|---|")
    verdicts = {}
    for label, s, small, large, d_max in PAIRS:
        vs, vl = values(draws.get(small, [])), values(draws.get(large, []))
        ms = statistics.median(v for v, _ in vs) if vs else None
        ml = statistics.median(v for v, _ in vl) if vl else None
        testable = (sum(e for _, e in vs) >= 3 and ms is not None and ms <= d_max and len(vl) >= 3)
        if not testable:
            why = "no committed draws" if not vl else "too few values"
            verdict = f"not testable ({large}: {why})"
        elif ml > ms:
            verdict = "grows"
        elif ml == ms and all(e for v, e in vl if v == ml):
            verdict = "flat"
        elif ml < ms:
            verdict = "falls"
        else:
            verdict = "not testable (median of L is a bound)"
        verdicts[s] = verdict
        fmt = lambda v: " ".join(f"{x}" if e else f"≥{x}" for x, e in v) or "—"
        print(f"| {label} | {s:+d} | {small} | {fmt(vs)} | {ms if ms is not None else '—'} | {large} | {fmt(vl)} "
              f"| {ml if ml is not None else '—'} | {verdict} |")
    primary = verdicts[-2]
    others = [v for s, v in verdicts.items() if s != -2 and not v.startswith("not testable")]
    if primary.startswith("not testable"):
        overall = "inconclusive (the primary pair is not testable)"
    elif primary == "grows" and all(v == "grows" for v in others):
        overall = "grows at fixed surplus"
    elif primary == "flat" and all(v == "flat" for v in others):
        overall = "flat"
    else:
        overall = "mixed"
    print(f"\nOverall, by the registered rule: {overall}")
    print("Testable pairs other than the primary:", {s: v for s, v in verdicts.items() if s != -2})
    print("\nSecondary (not decided on): FFD, controls, satisfiable fraction")
    for cell, rows in sorted(draws.items(), key=lambda kv: (kv[1][0]["n"], kv[1][0]["ell"])):
        uns = [r for r in rows if r["outcome"]["kind"] != "satisfiable"]
        sat = len(rows) - len(uns)
        print(f"  {cell}: ℓ={rows[0]['ell']} S={rows[0]['surplus']:+d} FFD {[r['ffd'] for r in uns]} "
              f"satisfiable {sat}/{len(rows)} control {controls.get(cell, 'unfinished')}")


if __name__ == "__main__":
    main()
