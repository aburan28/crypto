#!/usr/bin/env python3
"""Write one `ic bench` sweep file per (degree, columns) cell of the
double-large-prime measurement: LP in {0, 1, 2} x |D|/|F| in {3, 10, 30}
on the compact-orbit base with the three-summand Frobenius table.

    python3 make_sweeps.py            # all cells
    python3 make_sweeps.py 23 31      # these degrees only
"""
import json, sys, pathlib

here = pathlib.Path(__file__).resolve().parent
degrees = [int(a) for a in sys.argv[1:]] or [23, 31, 41]
# Repeats per column count: enough targets that every cell holds >= 32
# complete relations (a target pins after about columns + 1 rows).
repeats = {5: 8, 10: 4, 20: 4}
# Raw-x caps generous enough for the base and the largest domain.
raw_x_cap = {23: 4_000_000, 31: 400_000_000, 41: 400_000_000}
factors = [3, 10, 30]
for n in degrees:
    for columns in (5, 10, 20):
        oracles = ["mitm-frobenius:m=3"]
        for lp in (1, 2):
            for c in factors:
                oracles.append(
                    f"mitm-frobenius:m=3,large_primes={lp},lp_domain_factor={c},lp_raw_x_cap={raw_x_cap[n]}"
                )
        sweep = {
            "_comment": (
                f"Double large prime, n={n}, {columns} compact-orbit columns: the classic m=3 "
                f"Frobenius table against the large-prime collector with 1 or 2 residuals in a "
                f"domain of 3, 10 or 30 x |F| signed points. Operation counts are the metric."
            ),
            "instance": {"regime": "koblitz", "degree": n},
            "repeats": repeats[columns],
            "seed": 20261012,
            "max_trials": 40_000_000,
            "matrix": {
                "factor_base": [f"compact-orbit-scan:columns={columns},raw_x_cap={raw_x_cap[n]}"],
                "oracle": oracles,
                "targets": ["walk"],
            },
        }
        path = here / "sweeps" / f"dlp_n{n}_c{columns}.json"
        path.write_text(json.dumps(sweep, indent=1) + "\n")
        print(path)

# Equal-envelope control at n = 41: the classic m=3 table on a base as
# large as the LP cells' envelope F ∪ D (10 columns + |D|/|F| x 10 orbit
# columns: 40 columns for |D|/|F| = 3, 110 for 10), so the large-prime
# graph is compared with simply spending the same table on more columns.
if 41 in degrees:
    control = {
        "_comment": (
            "Equal-envelope control for the n=41, 10-column double-large-prime cells: the classic "
            "m=3 Frobenius table over 40 and 110 compact-orbit columns (the size of F ∪ D at "
            "|D|/|F| = 3 and 10), same seed and repeats."
        ),
        "instance": {"regime": "koblitz", "degree": 41},
        "repeats": repeats[10],
        "seed": 20261012,
        "max_trials": 40_000_000,
        "matrix": {
            "factor_base": [
                f"compact-orbit-scan:columns={c},raw_x_cap={raw_x_cap[41]}" for c in (40, 110)
            ],
            "oracle": ["mitm-frobenius:m=3"],
            "targets": ["walk"],
        },
    }
    path = here / "sweeps" / "dlp_n41_control_envelope.json"
    path.write_text(json.dumps(control, indent=1) + "\n")
    print(path)
