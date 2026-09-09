"""Independent GF(2^19) certificate for the full-universe selected base."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent
N = 19
FIELD = (1 << N) | (1 << 5) | (1 << 2) | (1 << 1) | 1


def remainder(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def multiply(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a


def field_mul(a, b):
    return remainder(multiply(a, b), FIELD)


def field_pow(a, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = field_mul(result, a)
        a = field_mul(a, a)
        exponent >>= 1
    return result


def trace(a):
    result = 0
    for _ in range(N):
        result ^= a
        a = field_mul(a, a)
    assert result in (0, 1)
    return result


def irreducible_degree_19(polynomial):
    assert polynomial.bit_length() - 1 == N
    x = 2
    for _ in range(N):
        x = remainder(multiply(x, x), polynomial)
    return x == 2 and gcd(polynomial, 0b110) == 1


def evaluate(polynomial, x):
    result = 0
    for i in reversed(range(polynomial.bit_length())):
        result = field_mul(result, x) ^ ((polynomial >> i) & 1)
    return result


def orbit(x):
    result = []
    for _ in range(N):
        result.append(x)
        x = field_mul(x, x)
    assert x == result[0] and len(set(result)) == N
    return result


def minimal_polynomial(x):
    coefficients = [1]
    for root in orbit(x):
        product = [0] * (len(coefficients) + 1)
        for i, coefficient in enumerate(coefficients):
            product[i] ^= field_mul(coefficient, root)
            product[i + 1] ^= coefficient
        coefficients = product
    assert all(coefficient in (0, 1) for coefficient in coefficients)
    polynomial = sum(coefficient << i for i, coefficient in enumerate(coefficients))
    assert irreducible_degree_19(polynomial)
    assert all(evaluate(polynomial, root) == 0 for root in orbit(x))
    return polynomial


assert len(sys.argv) == 2, "usage: full_orbit_predicate_certificate.py RESULT.jsonl"
assert irreducible_degree_19(FIELD)
result_path = Path(sys.argv[1]).resolve()
rows = [json.loads(line) for line in result_path.read_text().splitlines()]
selected = next(row for row in rows if row.get("kind") == "full_orbit_search_result")
verification = next(row for row in rows if row.get("kind") == "affine_quotient_base")
assert selected["candidate_orbits"] == 6909
assert selected["covered_target_orbits"] == 6300
assert selected["selected_indices"] == [200, 1913, 2481, 5643]
assert selected["selected_scalar_representatives"] == [203, 2143, 2901, 9853]
representatives = selected["selected_point_representatives"]
assert representatives == [[16795, 144921], [1315, 90425], [8461, 38471], [6685, 369649]]
assert verification["identity"]["point_representatives"] == representatives
assert verification["hits_at_most_three"] == 6300
assert verification["hits_exactly_three"] == 6256
assert verification["quotient_points"] == 152
assert verification["projected_signed_orbits"] == 4

roots = set()
factors = []
product = 1
for x, y in representatives:
    assert field_mul(y, y) ^ field_mul(x, y) == field_mul(field_mul(x, x), x) ^ field_mul(x, x) ^ 1
    new_roots = set(orbit(x))
    assert not roots & new_roots
    roots |= new_roots
    factor = minimal_polynomial(x)
    factors.append(
        {
            "u_representative": x,
            "binary_polynomial_mask": factor,
            "exponents": [i for i in range(N + 1) if (factor >> i) & 1],
        }
    )
    product = multiply(product, factor)

assert len(roots) == 76
assert product.bit_length() - 1 == 76
assert all(evaluate(product, root) == 0 for root in roots)
assert all(trace(root) == 1 and trace(field_pow(root, (1 << N) - 2)) == 0 for root in roots)
derivative = sum(1 << (i - 1) for i in range(1, product.bit_length(), 2) if (product >> i) & 1)
assert gcd(product, derivative) == 1
exponents = [i for i in range(product.bit_length()) if (product >> i) & 1]
certificate = {
    "kind": "full_universe_factor_base_membership_polynomial_certificate",
    "field_degree": N,
    "field_modulus_exponents": [0, 1, 2, 5, 19],
    "representative_points": representatives,
    "scalar_representatives": selected["selected_scalar_representatives"],
    "target_orbit_indices": selected["selected_indices"],
    "minimal_polynomials": factors,
    "membership_polynomial_degree": 76,
    "membership_polynomial_hex": hex(product),
    "membership_polynomial_exponents": exponents,
    "exact_field_root_count": len(roots),
    "root_representations": sorted(roots),
    "rational_quotient_point_count": 152,
    "signed_frobenius_orbit_count": 4,
    "covered_target_orbits_at_most_three": verification["hits_at_most_three"],
    "covered_target_orbits_exactly_three": verification["hits_exactly_three"],
    "target_orbits": verification["targets"],
    "target_orbit_size": verification["target_orbit_size"],
    "squarefree": True,
    "trace_u": 1,
    "trace_inverse_u": 0,
    "boolean_degree_upper_bound": max(i.bit_count() for i in exponents),
    "bound_scope": "syntactic degree bound for the descended membership equation; not solving degree or a complexity estimate",
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "selected_result_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
}
output = ROOT / "full-orbit-selected-predicate-certificate.json"
output.write_text(json.dumps(certificate, indent=2) + "\n")
print(
    json.dumps(
        {
            key: certificate[key]
            for key in [
                "minimal_polynomials",
                "membership_polynomial_hex",
                "membership_polynomial_exponents",
                "boolean_degree_upper_bound",
                "covered_target_orbits_at_most_three",
                "covered_target_orbits_exactly_three",
            ]
        }
    )
)
