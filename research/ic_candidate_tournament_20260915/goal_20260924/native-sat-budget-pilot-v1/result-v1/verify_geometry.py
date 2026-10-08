"""Postexecution independent Sage geometry check for the disclosed SAT pilot.

Run with /Volumes/SSD990/cryptanalysis/sage -python. This does not run a SAT
solver, estimate natural yield, or establish complete target admission.
"""

import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
PREPARATION = ROOT / "native-sat-control-registration-v1/result-v1/data/preparation.json"
RECORD = json.loads(PREPARATION.read_text())["record"]

F2 = GF(2)
t = PolynomialRing(F2, "t").gen()
field = GF(2**17, "z", modulus=t**17 + t**3 + 1)
z = field.gen()
curve = EllipticCurve(field, [1, 1, 0, 0, 1])


def element(bits):
    return sum(((bits >> i) & 1) * z**i for i in range(17))


def point(coords):
    return curve(element(coords[0]), element(coords[1]))


base = [point(coords) for coords in RECORD["factor_base"]["points"]]
assert len(base) == 63
pairs = {a + b for a in base for b in base}
checks = []
for trial, expected in ((0, False), (15, True)):
    result = json.loads((Path(__file__).parent / f"query-{trial:03d}/result.json").read_text())
    query = point(result["public_point"])
    possible = any(query - c in pairs for c in base)
    assert possible is expected
    assert result["geometric_class_valid"] is True
    assert result["prestarted_stdin_compatible"] is True
    for mode in ("cold", "stdin"):
        row = result[mode]
        assert row["timed_out"] is False
        assert row["status"] == ("SAT_MODEL" if expected else "SOURCE_UNSAT")
        if expected:
            model = row["model_check"]
            assert model["source_model_valid"] is True
            indices = model["full_point_witness_indices"]
            assert len(indices) == 3 and all(0 <= i < len(base) for i in indices)
            assert sum((base[i] for i in indices), curve(0)) == query
    checks.append({"trial": trial, "public_point": result["public_point"],
                   "geometric_three_sum_exists": possible,
                   "cold_status": result["cold"]["status"],
                   "stdin_status": result["stdin"]["status"]})

print(json.dumps({"schema_version": 1,
                  "status": "PASS_INDEPENDENT_SAGE_GEOMETRY_REPLAY",
                  "preparation_curve_id": RECORD["curve"]["curve"]["curve_id"],
                  "geometric_point_count": len(base), "checks": checks}, sort_keys=True))
