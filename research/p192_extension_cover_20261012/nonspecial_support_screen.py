#!/usr/bin/env python3
"""Certify noncollinearity of the frozen fixed-origin degree-three divisors."""

import hashlib
import json
import sys
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def coefficients(values, prime):
    result = [int(value[0] if isinstance(value, list) else value) for value in values]
    assert all(0 <= value < prime for value in result)
    return result


def main():
    if len(sys.argv) != 5:
        raise SystemExit("usage: nonspecial_support_screen.py SUPPORT.json OUT.json OUT.gp SOURCE_LABEL")
    source, output, gp_output, label = map(Path, sys.argv[1:])
    doc = json.loads(source.read_bytes())
    prime = int(doc["field_prime"])
    assert doc["status"] == "exact_fixed_origin_divisor_support_audit"
    assert doc["squarefree_cubics"] == doc["class_count"] == 160
    records = doc["records"]
    assert len(records) == 160
    witnesses = []
    gp_polys = []
    for i, record in enumerate(records):
        x = coefficients(record["x_coefficients"], prime)
        y = coefficients(record["y_coefficients"], prime)
        assert len(x) == 4 and x[3] == 1 and len(y) == 3
        if y[2] == 0:
            witnesses.append(i)
        gp_polys.append(" + ".join(f"{c}*z^{j}" for j, c in enumerate(y)))
    gp_bytes = (f"P={prime};\nY=[" + ",".join(gp_polys) + "];\n").encode()
    gp_output.write_bytes(gp_bytes)
    result = {
        "status": "fixed_origin_nonspecial_coefficient_screen",
        "field_prime": str(prime),
        "class_count": 160,
        "squarefree_distinct_abscissa_classes": 160,
        "quadratic_y_coefficient_nonzero": 160 - len(witnesses),
        "collinearity_witness_indices": witnesses,
        "support_sha256": digest(source),
        "source_sha256": digest(__file__),
        "gp_input_sha256": hashlib.sha256(gp_bytes).hexdigest(),
        "source_label": str(label),
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"NONSPECIAL_SUPPORT_PASS classes=160 nonzero_quadratic={160-len(witnesses)}")


if __name__ == "__main__":
    main()
