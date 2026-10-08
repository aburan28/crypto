#!/usr/bin/env python3
"""Derive automorphism floors for isogenous ECC2K-130 rho (no measurement)."""
from __future__ import annotations

import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent / "evidence" / "boundary.json"


def s_for(a: float) -> float:
    return math.sqrt(math.pi / (2.0 * a))


def main() -> None:
    rows = []
    for n in (13, 17, 19, 23, 31, 41, 53, 83, 131):
        a_e0 = 2 * n
        a_leaf = 2
        rows.append(
            {
                "n": n,
                "A_E0": a_e0,
                "A_leaf": a_leaf,
                "S_E0": s_for(a_e0),
                "S_leaf": s_for(a_leaf),
                "S_leaf_over_S_E0": math.sqrt(a_e0 / a_leaf),
                "prime_extension": True,  # all listed n are prime
                "intermediate_F2_subfields": [],
                "note": "only F_2-rational ordinary model in class is Koblitz when n is prime",
            }
        )
    n131 = next(r for r in rows if r["n"] == 131)
    payload = {
        "schema": "isogeny_rho_speed_boundary/v1",
        "unit": "S = charged_group_ops / sqrt(r)",
        "plain_rho_S": s_for(1.0),
        "negation_only_S": s_for(2.0),
        "ecc2k130": {
            "n": 131,
            "r_bits_approx": 129,
            "A_E0": 262,
            "S_E0": n131["S_E0"],
            "S_leaf": n131["S_leaf"],
            "iteration_ratio_leaf_over_E0": n131["S_leaf_over_S_E0"],
            "class_number_OK": 1,
            "horizontal_263_loops_to_j": 1,
            "descending_263_neighbors": 262,
            "derived_verdict": (
                "No reachable isogenous leaf beats E0 on rho iteration count: "
                "eligible A drops from 262 to 2, ratio sqrt(131)≈11.45."
            ),
        },
        "sizes": rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload["ecc2k130"], indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
