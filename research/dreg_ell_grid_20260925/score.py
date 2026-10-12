#!/usr/bin/env python3
"""Score the (n, ell) grid by PREREGISTRATION.md's rules, from committed runs only.

    python3 research/dreg_ell_grid_20260925/score.py

Reads this grid's runs/*.jsonl and the ladder's committed runs for the cells
it reuses.  A cell with fewer than three values is not testable; nothing here
infers, extrapolates or fills in a missing cell.
"""
import json
import statistics
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCES = [HERE / "runs", HERE.parent / "dreg_fixed_surplus_20260923" / "runs"]
Q1_PAIRS = [  # (N, ell = 2 cell, ell = 3 cell), as registered
    (13, "n7l2", "n4l3"),
    (14, "n8l2", "n5l3"),
    (16, "n10l2", "n7l3"),
    (18, "n12l2", "n9l3"),
]
Q2 = ("n13l4", "n10l5")  # both at N = 25


def load():
    draws = defaultdict(list)
    for src in SOURCES:
        for f in sorted(src.glob("cell-*.jsonl")):  # not the ladder's identity-check replays
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


def median(vs):
    return statistics.median(v for v, _ in vs) if len(vs) >= 3 else None


def exact_at(vs, m):
    return all(e for v, e in vs if v == m)


def fmt(vs):
    return " ".join(f"{v}" if e else f"≥{v}" for v, e in vs) or "—"


def main():
    draws = load()
    print("## Q1: at matched unknown count N, does ℓ = 3 resolve above ℓ = 2?\n")
    print("| N | ℓ = 2 cell | values | median | ℓ = 3 cell | values | median | pair |")
    print("|--:|---|---|--:|---|---|--:|---|")
    verdicts = []
    for N, c2, c3 in Q1_PAIRS:
        v2, v3 = values(draws.get(c2, [])), values(draws.get(c3, []))
        m2, m3 = median(v2), median(v3)
        if m2 is None or m3 is None:
            pair = "not testable (fewer than 3 values)"
        elif m3 > m2 and exact_at(v2, m2):
            pair = "ℓ-step"
        elif m3 == m2 and exact_at(v2, m2) and exact_at(v3, m3):
            pair = "no ℓ-step"
        elif m3 < m2 and exact_at(v3, m3):
            pair = "reversed"
        else:
            pair = "not testable (a bound decides it)"
        verdicts.append(pair)
        show = lambda m: "—" if m is None else f"{m:g}"
        print(f"| {N} | {c2} | {fmt(v2)} | {show(m2)} | {c3} | {fmt(v3)} | {show(m3)} | {pair} |")
    testable = [v for v in verdicts if not v.startswith("not testable")]
    if len(testable) < 3:
        q1 = f"inconclusive ({len(testable)} testable pairs; 3 are required)"
    elif all(v == "ℓ-step" for v in testable):
        q1 = "tracks ℓ, not size"
    elif all(v == "no ℓ-step" for v in testable):
        q1 = "tracks size, not ℓ"
    else:
        q1 = "mixed"
    print(f"\nQ1, by the registered rule: {q1}\n")

    print("## Q2: at N = 25, does ℓ = 5 resolve above ℓ = 4?\n")
    v4, v5 = values(draws.get(Q2[0], [])), values(draws.get(Q2[1], []))
    m5 = median(v5)
    if m5 is None:
        q2 = "not testable (fewer than 3 values)"
    elif m5 >= 7:
        q2 = "rises at ℓ = 5"
    elif m5 == 6 and exact_at(v5, 6):
        q2 = "flat through ℓ = 5"
    elif m5 < 6:
        q2 = "falls"
    else:
        q2 = "not testable (a bound decides it)"
    print(f"{Q2[0]}: {fmt(v4)}; {Q2[1]}: {fmt(v5)}")
    print(f"\nQ2, by the registered rule: {q2}\n")

    print("## Secondary (not decided on)\n")
    for cell, rows in sorted(draws.items(), key=lambda kv: (kv[1][0]["ell"], kv[1][0]["n_vars"])):
        uns = [r for r in rows if r["outcome"]["kind"] != "satisfiable"]
        sat = sum(1 for r in rows if r["outcome"]["kind"] == "satisfiable")
        secs = [r["secs"] for r in uns]
        drawn = f"satisfiable {sat}/{len(rows)}" if "unsat_index" not in rows[0] else "one process a draw"
        print(f"- {cell}: ℓ={rows[0]['ell']} N={rows[0]['n_vars']} S={rows[0]['surplus']:+d} "
              f"FFD {[r['ffd'] for r in uns]} {drawn} "
              f"secs {min(secs):.1f}–{max(secs):.1f}" if secs else f"- {cell}: no unsatisfiable draw")


if __name__ == "__main__":
    main()
