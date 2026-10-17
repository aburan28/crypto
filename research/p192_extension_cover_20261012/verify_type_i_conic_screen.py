"""Replay every P-192 conic-evaluation determinant from frozen half fibers."""

import hashlib
import json
from pathlib import Path

from sage.all import matrix, vector
import type_i_divisor_transfer as dt


HERE = Path(__file__).resolve().parent
receipt_path = HERE / "type_i_conic_screen.json"
half_path = HERE / "type_i_half_fiber_probe.json"
receipt = json.loads(receipt_path.read_text())
half = json.loads(half_path.read_text())
assert receipt["source_sha256"] == hashlib.sha256(
    (HERE / "type_i_conic_screen.py").read_bytes()
).hexdigest()
assert receipt["half_fiber_receipt_sha256"] == hashlib.sha256(half_path.read_bytes()).hexdigest()
assert receipt["field_prime"] == str(dt.p)
assert len(receipt["classes"]) == len(half["classes"]) == 16


def decode(encoded):
    return tuple(dt.L([int(c) for c in coordinate]) for coordinate in encoded)


def determinant(encoded_fiber):
    rows = []
    for power in range(3):
        for encoded in encoded_fiber:
            point = tuple(a ** (dt.p ** power) for a in decode(encoded))
            ys = dt.check_model(point)
            u, v, w = dt.Vandermonde.solve_right(vector(dt.L, [1 / y for y in ys]))
            rows.append((u * u, v * v, w * w, u * v, u * w, v * w))
    evaluation = matrix(dt.L, rows)
    return evaluation.det(), int(evaluation.rank())


for item, source in zip(receipt["classes"], half["classes"]):
    assert item["base_index"] == source["base_index"]
    for side, fiber_key in (("positive", "target_fiber"),
                            ("negative", "half_fiber")):
        value, rank = determinant(source[fiber_key])
        reported = item[side]
        assert rank == reported["rank"] == 6
        assert value != 0 and reported["determinant_nonzero"]
        assert [str(c) for c in value.polynomial().list()] == reported[
            "determinant_coefficients"
        ]

out = {
    "status": "verified_full_conic_evaluation_rank_for_16_type_i_divisors",
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "producer_receipt_sha256": hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
    "classes_replayed": 16,
    "nonzero_determinants_replayed": 32,
    "boundary": "Six-point conic interpolation only; no general Jacobian reduction certificate.",
}
(HERE / "type_i_conic_verification.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps(out, indent=2, sort_keys=True))
