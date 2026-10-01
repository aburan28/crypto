#!/usr/bin/env python3
"""Ledger §23's declared prediction, computed before anything in §23 ran.

§20's model and law, unchanged, carried to §23's six sizes: the four
largest of §20's ladder and the only two other Koblitz curves with
`n ≤ MAX_N = 63` and `r ≥ 2^36`.  The model's constants are §20's frozen
ones (`research/ic_exponent_20260926/predict.py`); nothing here is fitted
to a §23 measurement.  Each size's sweep grid is §20's rule: the model's
optimum column count times `2^{j/2}`, `j = −2 … 4`.

    python3 research/ic_exponent_top_20260930/predict.py > research/ic_exponent_top_20260930/prediction.json
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("s20", HERE.parent / "ic_exponent_20260926" / "predict.py")
S20 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S20)

# (a, n, r): r is the subgroup order KoblitzCurve::new builds, read with
# examples/koblitz_curve_records.rs (curve_records.json).
SIZES = [
    (1, 47, 106781081677),
    (0, 57, 275295876199),
    (0, 41, 549756390943),
    (0, 53, 21044858204113),
    (1, 59, 25179555920633),
    (0, 61, 162888033982417),
]
NEW = {(0, 57), (1, 59)}


def fit(rows: list[dict]) -> float:
    xs = [math.log(row["r"]) for row in rows]
    ys = [math.log(row["model_ratio"] * math.sqrt(row["n"])) for row in rows]
    xm, ym = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / sum((x - xm) ** 2 for x in xs)


def main() -> None:
    rows = []
    for a, n, r in SIZES:
        s, cols, m = S20.optimum(r, n)
        F = 2 * n * cols
        phases = S20.ic_phases(F, r, n, m)
        total = sum(phases.values())
        grid = sorted({max(2, round(cols * 2 ** (j / 2))) for j in range(-2, 5)})
        rows.append({
            "curve": f"K_{a}/GF(2^{n})", "a": a, "n": n, "r": r, "log2_r": round(math.log2(r), 2),
            "new_in_23": (a, n) in NEW,
            "model_optimum": {"columns": cols, "points": F, "descent_summands": m},
            "model_S_ic": round(s, 4),
            "model_S_batch_rho": round(S20.rho_s(n), 4),
            "model_ratio": round(s / S20.rho_s(n), 2),
            "model_phase_shares": {k: round(v / total, 3) for k, v in phases.items()},
            "law_ratio": round(S20.law_ratio(r, n), 2),
            "sweep_grid_columns": grid,
        })
    top4 = sorted(rows, key=lambda row: row["r"])[-4:]
    report = {
        "what_this_is": "Ledger §23's declared prediction: §20's frozen model and law at §23's six sizes, computed before any §23 measurement.",
        "source_model": "research/ic_exponent_20260926/predict.py (constants unchanged)",
        "k": S20.K,
        "law_exponent": 1 / 6,
        "model_local_exponent_six_sizes_n_corrected": round(fit(rows), 3),
        "model_local_exponent_top_four_n_corrected": round(fit(top4), 3),
        "rows": rows,
    }
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
