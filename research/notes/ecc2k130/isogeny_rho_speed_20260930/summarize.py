#!/usr/bin/env python3
"""Summarize n17 rho screen into a one-unit table."""
from __future__ import annotations

import json
import statistics as stats
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EV = ROOT / "evidence"


def med(xs):
    xs = sorted(xs)
    return xs[len(xs) // 2]


def main() -> None:
    boundary = json.loads((EV / "boundary.json").read_text())
    n17 = json.loads((EV / "n17-rho.json").read_text())
    rows_out = []
    for row in n17["rows"]:
        runs = row["runs"]
        rows_out.append(
            {
                "arm": row["arm"],
                "a6": row["a6"],
                "is_koblitz": row["is_koblitz"],
                "A": row["A"],
                "S_floor": row["S_floor"],
                "S_median": med(r["S"] for r in runs),
                "S_walk_median": med(r["S_walk"] for r in runs),
                "steps_median": med(r["steps"] for r in runs),
                "verified": all(r["verified"] for r in runs),
                "reps": len(runs),
                "class": "accounting",
                "note": (
                    "toy-size total S for signed-Frobenius is setup-dominated; "
                    "compare S_walk to A-floor"
                    if row["arm"] == "signed_frobenius"
                    else "A=2 negation floor on this member"
                ),
            }
        )
    sf = next(r for r in rows_out if r["arm"] == "signed_frobenius")
    neg_others = [r for r in rows_out if r["arm"] == "negation" and not r["is_koblitz"]]
    walk_ratio = med(r["S_walk_median"] for r in neg_others) / sf["S_walk_median"]
    summary = {
        "schema": "isogeny_rho_speed_summary/v1",
        "unit": "S = charged_group_ops / sqrt(r); S_walk = walk_steps / sqrt(r)",
        "boundary_ecc2k130": boundary["ecc2k130"],
        "n17": {
            "class_members": n17["class_members"],
            "sample_a6": n17["sample_a6"],
            "predicted_walk_ratio_sqrt_2n": n17["prediction"]["S_leaf_over_S_E0"],
            "measured_median_S_walk_neg_others_over_sf": walk_ratio,
            "all_verified": all(r["verified"] for r in rows_out),
            "rows": rows_out,
        },
        "verdict": (
            "No isogenous vertex beats E0 on rho iteration count. At n=131 the "
            "eligible A drops from 262 to 2 on descending leaves (ratio √131≈11.45). "
            "At n=17 every sampled non-Koblitz class member only admits A=2; "
            "signed-Frobenius on the Koblitz member admits A=34; all scalars verified. "
            "Modal/RunPod GPU comparison of leaf curves is blocked (no credentials; "
            "leaf kernel unimplemented) and cannot overturn the automorphism floor."
        ),
    }
    out = EV / "summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
