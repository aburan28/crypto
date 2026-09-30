#!/usr/bin/env python3
"""Score the degree-7 study at ℓ = 5 by the rules in PREREGISTRATION.md.

    python3 research/dreg_degree7_ell5_20260930/score.py            # score runs/
    python3 research/dreg_degree7_ell5_20260930/score.py --selftest # the rules on made-up values

Every degree-7 row is first matched to its committed degree-6 row: the
same draw index, subspace and target, and a degree-6 outcome of `≥7`.  A
row that fails either check is void, because `--d-min 7` is only valid on
a draw already refuted by no degree below 7.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = HERE.parent
# cell -> (surplus, committed degree-6 run files, role)
CELLS = {
    "10-5-7": (-5, "dreg_ell_grid_20260925/runs/cell-10-5-6.u{k}.jsonl", "Q6"),
    "13-5-7": (-2, "dreg_fixed_surplus_20260923/runs/cell-13-5-6.u{k}.jsonl", "Q5"),
    "11-5-7": (-4, "dreg_surplus_control_20260925/runs/cell-11-5-6.u{k}.jsonl", "secondary"),
}
# the committed ℓ = 3 and ℓ = 4 values at each surplus, for the growth column
BELOW = {-5: ("5 5 5 5", "7 7 6 7"), -4: ("6 6 6 6", "6 7 6 7"), -2: ("6 6 6 6", "—")}


def value(outcome):
    """A draw's value: (degree, exact).  None if it is not evidence."""
    if outcome["kind"] == "resolved":
        return (outcome["degree"], True)
    if outcome["kind"] == "at_least":
        return (outcome["degree"], False)
    return None  # caps_hit, satisfiable (impossible here): excluded


def fmt(v):
    return f"{v[0]}" if v[1] else f"≥{v[0]}"


def cell_reading(values):
    """The registered cell reading from its values (a bound counts as its bound)."""
    vals = [v for v in values if v is not None]
    if len(vals) < 3:
        return "not testable"
    xs = sorted(v[0] for v in vals)
    mid = len(xs) // 2
    med = xs[mid] if len(xs) % 2 else (xs[mid - 1] + xs[mid]) / 2
    if med == 7:
        return "7"
    if med >= 8:
        return "≥8"
    if med == 7.5:
        return "split 7 / ≥8"
    return f"median {med}: no registered reading"


def verdict(role, reading):
    if reading == "7":
        return f"{role}: one degree above 6"
    if reading == "≥8":
        return f"{role}: at least two degrees above 6"
    if reading.startswith("split"):
        return f"{role}: split between one and at least two above 6"
    return f"{role}: {reading}"


def score(runs):
    table = []
    for cell, (surplus, ref_pat, role) in CELLS.items():
        vals, shown, notes = [], [], []
        for k in range(4):
            f = runs / f"cell-{cell}.u{k}.jsonl"
            if not f.exists() or not f.read_text().strip():
                shown.append("·")
                continue
            got = json.loads(f.read_text().splitlines()[-1])
            ref = json.loads((R / ref_pat.format(k=k)).read_text().splitlines()[-1])
            same = all(got[x] == ref[x] for x in ("draw", "v_basis", "x_r", "n", "ell"))
            premise = ref["outcome"] == {"kind": "at_least", "degree": 7}
            if not (same and premise and got.get("d_min") == 7):
                shown.append("void")
                notes.append(f"u{k} void: same draw {same}, degree-6 premise {premise}, d_min {got.get('d_min')}")
                continue
            v = value(got["outcome"])
            vals.append(v)
            shown.append(fmt(v) if v else json.dumps(got["outcome"]))
            notes.append(f"u{k} draw {got['draw']}: {shown[-1]} ffd {got['ffd']} secs {got['secs']:.0f}")
        reading = cell_reading(vals)
        table.append((cell, surplus, role, " ".join(shown), reading, notes))
    print("| cell | S | ℓ = 3 | ℓ = 4 | ℓ = 5 at degree 7 | reading | verdict |")
    print("|---|--:|---|---|---|---|---|")
    for cell, surplus, role, shown, reading, _ in table:
        l3, l4 = BELOW[surplus]
        n, l, _ = cell.split("-")
        print(f"| ({n}, {l}) | {surplus} | {l3} | {l4} | {shown} | {reading} | {verdict(role, reading)} |")
    for cell, *_, notes in table:
        for line in notes:
            print(f"  {cell} {line}")


def selftest():
    e, b = (lambda d: (d, True)), (lambda d: (d, False))
    cases = [
        ([e(7)] * 4, "7"),
        ([b(8)] * 4, "≥8"),
        ([e(7), e(7), b(8), b(8)], "split 7 / ≥8"),
        ([e(7), e(7), e(7), b(8)], "7"),
        ([e(7), b(8), b(8), b(8)], "≥8"),
        ([e(7), e(7), None, None], "not testable"),
        ([e(7), e(7), b(8), None], "7"),
    ]
    for vals, want in cases:
        got = cell_reading(vals)
        assert got == want, (vals, got, want)
    print(f"selftest: {len(cases)} cases pass")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    else:
        score(Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "runs")
