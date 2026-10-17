#!/usr/bin/env python3
"""Exact support-polynomial audit for 160 fixed-origin Type I divisors."""

import collections
import hashlib
import json
import sys
from pathlib import Path


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def trim(poly):
    poly = list(poly)
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def remainder(poly, modulus, prime):
    poly = trim(poly)
    inverse = pow(modulus[-1], -1, prime)
    while len(poly) >= len(modulus):
        shift = len(poly) - len(modulus)
        scale = poly[-1] * inverse % prime
        for i, value in enumerate(modulus):
            poly[shift + i] = (poly[shift + i] - scale * value) % prime
        poly = trim(poly)
    return poly


def gcd(left, right, prime):
    while right:
        left, right = right, remainder(left, right, prime)
    if not left:
        return []
    inverse = pow(left[-1], -1, prime)
    return [value * inverse % prime for value in left]


def product(left, right, modulus, prime):
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] = (out[i + j] + a * b) % prime
    return remainder(out, modulus, prime)


def power(base, exponent, modulus, prime):
    result = [1]
    while exponent:
        if exponent & 1:
            result = product(result, base, modulus, prime)
        base = product(base, base, modulus, prime)
        exponent >>= 1
    return result


def coefficients(values, prime):
    result = []
    for value in values:
        if isinstance(value, list):
            assert len(value) == 1
            value = value[0]
        assert isinstance(value, str) and value.isdecimal()
        result.append(int(value))
    assert all(0 <= value < prime for value in result)
    return result


def main():
    if len(sys.argv) != 6:
        raise SystemExit("usage: divisor_support_screen.py POS NEG GUIDED_NEG OUT.json OUT.gp")
    raw = [Path(path).read_bytes() for path in sys.argv[1:4]]
    docs = [json.loads(value) for value in raw]
    prime = int(docs[0]["field_prime"])
    assert all(int(doc["field_prime"]) == prime for doc in docs)
    positive, negative, guided_negative = [doc["representatives"] for doc in docs]
    assert (len(positive), len(negative), len(guided_negative)) == (80, 16, 64)
    grouped = [
        ("base_positive", positive[:16]),
        ("guided_positive", positive[16:]),
        ("base_negative", negative),
        ("guided_negative", guided_negative),
    ]
    records = []
    cubic_vectors = {"numerator": set(), "denominator": set()}
    for group, rows in grouped:
        for index, row in enumerate(rows):
            x = coefficients(row["representative_x_polynomial"], prime)
            y = coefficients(row["representative_y_polynomial"], prime)
            assert len(x) == 4 and x[-1] == 1 and len(y) == 3
            for label in cubic_vectors:
                cubic_vectors[label].add(tuple(coefficients(row[f"{label}_cubic_coefficients"], prime)))
            derivative = [x[i] * i % prime for i in range(1, 4)]
            assert len(gcd(x, trim(derivative), prime)) == 1
            xp = power([0, 1], prime, x, prime)
            frobenius_difference = trim([
                ((xp[i] if i < len(xp) else 0) - (1 if i == 1 else 0)) % prime
                for i in range(3)
            ])
            linear_roots = len(gcd(x, frobenius_difference, prime)) - 1
            assert linear_roots in (0, 1, 3)
            records.append({"group": group, "index": index,
                            "linear_roots": linear_roots,
                            "x_coefficients": [str(value) for value in x],
                            "y_coefficients": [str(value) for value in y]})
    assert len(records) == 160
    degree_counts = collections.Counter()
    first_shared = []
    for i, left in enumerate(records):
        x_left = list(map(int, left["x_coefficients"]))
        for j in range(i + 1, len(records)):
            degree = len(gcd(x_left, list(map(int, records[j]["x_coefficients"])), prime)) - 1
            degree_counts[degree] += 1
            if degree and len(first_shared) < 20:
                first_shared.append([i, j, degree])
    assert sum(degree_counts.values()) == 160 * 159 // 2
    counts = {group: dict(sorted(collections.Counter(
        row["linear_roots"] for row in records if row["group"] == group).items()))
        for group, _ in grouped}
    gp_polynomials = []
    for row in records:
        x = row["x_coefficients"]
        gp_polynomials.append(" + ".join(f"{value}*x^{i}" for i, value in enumerate(x)))
    gp_lines = [f"P={prime};", "X=[" + ",".join(gp_polynomials) + "];", ""]
    gp_bytes = "\n".join(gp_lines).encode()
    Path(sys.argv[5]).write_bytes(gp_bytes)
    output = {
        "status": "exact_fixed_origin_divisor_support_audit",
        "field_prime": str(prime),
        "class_count": len(records),
        "group_counts": {group: len(rows) for group, rows in grouped},
        "linear_root_count_by_group": counts,
        "squarefree_cubics": len(records),
        "distinct_numerator_cubic_sections": len(cubic_vectors["numerator"]),
        "distinct_denominator_cubic_sections": len(cubic_vectors["denominator"]),
        "pairwise_x_gcd_degree_counts": dict(sorted(degree_counts.items())),
        "first_shared_x_witnesses": first_shared,
        "records": records,
        "input_sha256": list(map(sha, raw)),
        "source_sha256": sha(Path(__file__).read_bytes()),
        "gp_input_sha256": sha(gp_bytes),
    }
    Path(sys.argv[4]).write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print("DIVISOR_SUPPORT_PASS classes=160 shared_x_pairs="
          f"{sum(n for d, n in degree_counts.items() if d)} "
          f"split={sum(row['linear_roots'] == 3 for row in records)}")


if __name__ == "__main__":
    main()
