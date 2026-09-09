"""Independent GF(2) / GF(2^19) certificate for the selected factor-base predicate."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
N = 19
FIELD = (1 << N) | (1 << 5) | (1 << 2) | (1 << 1) | 1

def remainder(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a

def multiply(a, b):
    z = 0
    while b:
        if b & 1:
            z ^= a
        a <<= 1
        b >>= 1
    return z

def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a

def field_mul(a, b):
    return remainder(multiply(a, b), FIELD)

def field_pow(a, e):
    z = 1
    while e:
        if e & 1:
            z = field_mul(z, a)
        a = field_mul(a, a)
        e >>= 1
    return z

def trace(a):
    z = 0
    for _ in range(N):
        z ^= a
        a = field_mul(a, a)
    assert z in (0, 1)
    return z

def irreducible_degree_19(f):
    assert f.bit_length() - 1 == N
    x = 2
    for _ in range(N):
        x = remainder(multiply(x, x), f)
    return x == 2 and gcd(f, 0b110) == 1

def evaluate(f, x):
    z = 0
    for i in reversed(range(f.bit_length())):
        z = field_mul(z, x) ^ ((f >> i) & 1)
    return z

def orbit(x):
    result = []
    for _ in range(N):
        result.append(x)
        x = field_mul(x, x)
    assert x == result[0] and len(set(result)) == N
    return result

def minimal_polynomial(x):
    # Coefficients initially lie in GF(2^19); invariance must reduce them to GF(2).
    coefficients = [1]
    for r in orbit(x):
        product = [0] * (len(coefficients) + 1)
        for i, c in enumerate(coefficients):
            product[i] ^= field_mul(c, r)
            product[i + 1] ^= c
        coefficients = product
    assert all(c in (0, 1) for c in coefficients)
    f = sum(c << i for i, c in enumerate(coefficients))
    assert irreducible_degree_19(f)
    assert all(evaluate(f, r) == 0 for r in orbit(x))
    return f

assert irreducible_degree_19(FIELD)
rows = [json.loads(line) for line in (ROOT / 'orbit-global-n19.jsonl').read_text().splitlines()]
finalists = [r for r in rows if r.get('exhaustive_via_symmetry')]
assert len(finalists) == 1
selected = max(finalists, key=lambda r: r['hits_at_most_three'])
representatives = selected['identity']['representatives']
assert [p[0] for p in representatives] == [6685, 8461, 25103, 32881]
roots = set()
factors = []
product = 1
for x, y in representatives:
    # Independent verification of the supplied curve point.
    assert field_mul(y, y) ^ field_mul(x, y) == field_mul(field_mul(x, x), x) ^ field_mul(x, x) ^ 1
    new_roots = set(orbit(x))
    assert not roots & new_roots
    roots |= new_roots
    f = minimal_polynomial(x)
    factors.append({'u_representative': x, 'binary_polynomial_mask': f,
                    'exponents': [i for i in range(N + 1) if (f >> i) & 1]})
    product = multiply(product, f)
assert len(roots) == 76
assert product.bit_length() - 1 == 76
assert all(evaluate(product, x) == 0 for x in roots)
assert all(trace(x) == 1 and trace(field_pow(x, (1 << N) - 2)) == 0 for x in roots)
derivative = sum(1 << (i - 1) for i in range(1, product.bit_length(), 2) if (product >> i) & 1)
assert gcd(product, derivative) == 1
exponents = [i for i in range(product.bit_length()) if (product >> i) & 1]
certificate = {
    'kind': 'factor_base_membership_polynomial_certificate',
    'field_degree': N,
    'field_modulus_exponents': [0, 1, 2, 5, 19],
    'representative_points': representatives,
    'minimal_polynomials': factors,
    'membership_polynomial_degree': 76,
    'membership_polynomial_hex': hex(product),
    'membership_polynomial_exponents': exponents,
    'exact_field_root_count': len(roots),
    'root_representations': sorted(roots),
    'rational_quotient_point_count': 152,
    'signed_frobenius_orbit_count': 4,
    'squarefree': True,
    'trace_u': 1,
    'trace_inverse_u': 0,
    'boolean_degree_upper_bound': max(i.bit_count() for i in exponents),
    'bound_scope': 'syntactic degree bound for the descended membership equation; not solving degree or a complexity estimate',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'selected_result_sha256': hashlib.sha256((ROOT / 'orbit-global-n19.jsonl').read_bytes()).hexdigest(),
}
(ROOT / 'global-selected-predicate-certificate.json').write_text(json.dumps(certificate, indent=2) + '\n')
print(json.dumps({k: certificate[k] for k in ['minimal_polynomials', 'membership_polynomial_hex', 'membership_polynomial_exponents', 'boolean_degree_upper_bound']}))
