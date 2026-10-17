"""Conic-evaluation screen for the 16 unit-coefficient Type I divisors."""

import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, matrix, prod, vector
import type_i_divisor_transfer as dt


HERE = Path(__file__).resolve().parent
input_path = HERE / "type_i_half_fiber_probe.json"
plane_path = HERE / "type_i_plane_quartic.json"
source = json.loads(input_path.read_text())
plane = json.loads(plane_path.read_text())
p, L, F = dt.p, dt.L, dt.gm.F
R = PolynomialRing(F, 3, names=("Z0", "Z1", "Z2"))
Z = R.gens()
quartic = sum(F(int(c)) * prod(Z[i] ** e[i] for i in range(3))
              for key, c in plane["quartic_coefficients"].items()
              for e in [tuple(map(int, key.split(",")))])


def decode(encoded):
    return tuple(L([int(c) for c in a]) for a in encoded)


def coordinates(point):
    values = dt.check_model(point)
    assert all(y != 0 for y in values)
    result = dt.Vandermonde.solve_right(vector(L, [1 / y for y in values]))
    assert quartic(*result) == 0
    return result


def screen(seed):
    points = [tuple(a ** (p ** i) for a in decode(encoded))
              for i in range(3) for encoded in seed]
    rows = []
    for point in points:
        u, v, w = coordinates(point)
        rows.append([u * u, v * v, w * w, u * v, u * w, v * w])
    evaluation = matrix(L, rows)
    determinant = evaluation.det()
    return {"rank": int(evaluation.rank()),
            "determinant_coefficients": [str(c) for c in determinant.polynomial().list()],
            "determinant_nonzero": determinant != 0}


start = time.monotonic_ns()
classes = []
for item in source["classes"]:
    assert item["unit_divisor_verified"]
    positive = screen(item["target_fiber"])
    negative = screen(item["half_fiber"])
    assert positive["rank"] == negative["rank"] == 6
    assert positive["determinant_nonzero"] and negative["determinant_nonzero"]
    classes.append({"base_index": item["base_index"],
                    "positive": positive, "negative": negative})
    print(f"conic rank checked {item['base_index'] + 1}/16", flush=True)

out = {
    "status": "verified_full_conic_evaluation_rank_for_16_type_i_divisors",
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "half_fiber_receipt_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
    "plane_receipt_sha256": hashlib.sha256(plane_path.read_bytes()).hexdigest(),
    "field_prime": str(p),
    "monomial_order": ["Z0^2", "Z1^2", "Z2^2", "Z0*Z1", "Z0*Z2", "Z1*Z2"],
    "classes": classes,
    "full_rank_positive": sum(c["positive"]["rank"] == 6 for c in classes),
    "full_rank_negative": sum(c["negative"]["rank"] == 6 for c in classes),
    "elapsed_ns": time.monotonic_ns() - start,
    "boundary": "No nonzero conic vanishes on either complete six-point support set for these 16 divisors; other Riemann-Roch reductions remain open.",
}
(HERE / "type_i_conic_screen.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps({k: v for k, v in out.items() if k != "classes"}, sort_keys=True))
