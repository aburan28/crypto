#!/usr/bin/env python3
"""Ledger §20's declared prediction, computed before anything in §20 ran.

Every constant below is read from a frozen file that predates §20, and each
one names its source.  Nothing here is fitted to a §20 measurement: this
script and its output, prediction.json, are committed with the protocol
and before the pricer or any run exists.

Two predictions are made for the same question: at k = 32 targets, the
Koblitz collection thread's m = 3 method, with every phase priced, against
batch rho at the same k.

1. The page's law.  The scoreboard derives an r^(1/6) rise for this family,
   and §19.5 closes at 6.38x batch rho at r = 2^39 (n = 41), both sides
   priced.  Carried along r^(1/6) n^(-1/2) from that anchor, the ratio
   reaches one near r = 2^23.
2. The phase model.  The same shared costs, plus a descent priced at the
   probes it actually scans.  Each m = 3 descent trial is a full scan of the
   base in blocks of 1,024, and each m = 2 trial is one walked probe in
   lockstep rounds of 64.  The frozen ledger priced the descent per trial
   (54 units for 32 targets at the headline).  Priced per probe, the
   descent's cost per target does not amortise over k and grows as r
   falls, so the model predicts a minimum ratio rather than a crossing.

    python3 research/ic_exponent_20260926/predict.py > research/ic_exponent_20260926/prediction.json
"""
from __future__ import annotations

import json
import math

K = 32

# Sizes declared in PROTOCOL.md (a, n, r), r from examples/koblitz_degree_census.rs.
SIZES = [
    (1, 19, 262543),
    (1, 23, 4196903),
    (1, 45, 29264761),
    (0, 37, 230603167),
    (1, 43, 4644189029),
    (1, 47, 106781081677),
    (0, 41, 549756390943),
    (0, 53, 21044858204113),
    (0, 61, 162888033982417),
]

CONSTANTS = {
    "scan_units_per_summand": {
        "value": 2.85,
        "source": "docs/ic/runs/koblitz-collection-aim-20260922.json conversions_measured.scan_adds_per_summand.by_base_then_window['15744']['256']",
    },
    "summands_per_relation_times_F2_over_r": {
        "value": 883200 / 197 * 15744**2 / 549756390943,
        "source": "same run, arm aimed_at_least_mentioned: 883,200 summands scanned for 197 relations at |F| = 15,744, r = 549,756,390,943",
    },
    "relations_over_counting_floor": {
        "value": 197 / 192,
        "source": "same arm: 197 relations for 192 columns",
    },
    "build_units_per_stored_pair": {
        "value": 6.14,
        "source": "ledger §19.4, n = 41, |F| = 15,744, at commit 3e8dd352",
    },
    "select_units_per_point": {
        "value": 37.3,
        "source": "docs/ic/runs/koblitz-select-packed-20260922.json measurements.selection_at_the_ledger_width.after_adds_per_point",
    },
    "descent_units_per_probe_m3": {
        "value": 2.67,
        "source": "docs/ic/runs/koblitz-collection-aim-20260922.json scan_adds_per_summand '15744' 'full' (the descent's m = 3 trial is a full scan)",
    },
    "descent_units_per_probe_m2": {
        "value": 1.19,
        "source": "docs/ic/runs/koblitz-collection-aim-20260922.json conversions_measured.descent_probe_adds (one walked probe)",
    },
    "rho_step_units": {
        "value": 2.83,
        "source": "ledger §19.4 canonical step, 2.74 / 2.83 / 2.92 at n = 41 / 53 / 61; the middle value",
    },
    "m3_block": {"value": 1024, "source": "PairSumTable::witnesses_fast_inner, BLOCK"},
    "m2_round": {"value": 64, "source": "IndividualLogSolver::solve_by_walking, WALKS"},
    "anchor_ratio": {
        "value": 6.38,
        "source": "ledger §19.5 headline, both corrections, canonical step: 6.38x batch rho at k = 32, n = 41",
    },
}


def c(name: str) -> float:
    return CONSTANTS[name]["value"]


def batch_law(k: int) -> float:
    """Kuhn–Struik: the i-th target's expected share, averaged over k."""
    return sum(math.comb(2 * i, i) / 4**i for i in range(k)) / k


def descent_probes(mu: float, granule: float) -> float:
    """Probes to the first witness when probes are paid `granule` at a
    time and the witness arrives after an exponential number with mean
    `mu`: `granule / (1 - exp(-granule / mu))`."""
    return granule / -math.expm1(-granule / mu)


def ic_phases(F: int, r: int, n: int, m: int) -> dict[str, float]:
    columns = F / (2 * n)
    mu = c("summands_per_relation_times_F2_over_r") * r / F**2
    collect = c("relations_over_counting_floor") * columns * mu * c("scan_units_per_summand")
    build = (F * F / (4 * n) + F) * c("build_units_per_stored_pair")
    select = F * c("select_units_per_point")
    if m == 3:
        per_target = descent_probes(mu, min(c("m3_block"), F)) * c("descent_units_per_probe_m3")
    else:
        per_target = descent_probes(mu, c("m2_round")) * c("descent_units_per_probe_m2")
    return {"collect": collect, "build": build, "select": select, "descent": K * per_target}


def ic_s(F: int, r: int, n: int, m: int) -> float:
    return sum(ic_phases(F, r, n, m).values()) / (K * math.sqrt(r))


def rho_s(n: int) -> float:
    return batch_law(K) * math.sqrt(math.pi / (4 * n)) * c("rho_step_units")


def optimum(r: int, n: int) -> tuple[float, int, int]:
    best = None
    for m in (2, 3):
        for cols in range(2, 20000):
            s = ic_s(2 * n * cols, r, n, m)
            if best is None or s < best[0]:
                best = (s, cols, m)
    return best


def law_ratio(r: int, n: int) -> float:
    return c("anchor_ratio") * (r / 549756390943) ** (1 / 6) * (n / 41) ** -0.5


def main() -> None:
    rows = []
    for a, n, r in SIZES:
        s, cols, m = optimum(r, n)
        F = 2 * n * cols
        phases = ic_phases(F, r, n, m)
        total = sum(phases.values())
        grid = sorted({max(2, round(cols * 2 ** (j / 2))) for j in range(-2, 5)})
        rows.append({
            "curve": f"K_{a}/GF(2^{n})", "a": a, "n": n, "r": r, "log2_r": round(math.log2(r), 2),
            "model_optimum": {"columns": cols, "points": F, "descent_summands": m},
            "model_S_ic": round(s, 4),
            "model_S_batch_rho": round(rho_s(n), 4),
            "model_ratio": round(s / rho_s(n), 2),
            "model_phase_shares": {k: round(v / total, 3) for k, v in phases.items()},
            "law_ratio": round(law_ratio(r, n), 2),
            "sweep_grid_columns": grid,
        })
    ratios = [(row["log2_r"], row["model_ratio"]) for row in rows]
    lo = min(ratios, key=lambda t: t[1])
    top = [row for row in rows if row["log2_r"] >= 36]
    xs = [math.log(row["r"]) for row in top]
    ys = [math.log(row["model_ratio"] * math.sqrt(row["n"])) for row in top]
    xm, ym = sum(xs) / len(xs), sum(ys) / len(ys)
    beta = sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / sum((x - xm) ** 2 for x in xs)
    report = {
        "what_this_is": "Ledger §20's declared prediction, computed from frozen constants before any §20 measurement.",
        "k": K,
        "batch_law_k32": batch_law(K),
        "constants": CONSTANTS,
        "model": {
            "S_ic": "[collect + build + select] / (k sqrt r) + descent_per_target / sqrt r, minimised over |F| = 2n x columns (columns >= 2) and descent summands in {2, 3}",
            "collect": "relations_over_counting_floor x |F|/2n x mu x scan_units_per_summand, mu = summands_per_relation_times_F2_over_r x r / |F|^2",
            "build": "(|F|^2 / 4n + |F|) x build_units_per_stored_pair",
            "select": "|F| x select_units_per_point",
            "descent_per_target": "g / (1 - exp(-g / mu)) probes, g = min(1024, |F|) for m = 3 and 64 for m = 2, times that descent's units per probe",
            "S_batch_rho": "batch_law(32) x sqrt(pi / 4n) x rho_step_units",
            "not_modelled": "linear algebra (0.06% of the frozen headline), verification, setup",
        },
        "law": {
            "ratio": "anchor_ratio x (r / 2^39.0)^(1/6) x (n / 41)^(-1/2)",
            "crossing_log2_r_at_n41": round(math.log2(549756390943) - 6 * math.log2(c("anchor_ratio")), 2),
        },
        "model_minimum": {"log2_r": lo[0], "ratio": lo[1]},
        "model_local_exponent_top_four_n_corrected": round(beta, 3),
        "rows": rows,
    }
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
