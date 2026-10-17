"""Reduce the P-192 Type I conorm and its cubic Frobenius norm by elliptic quotients.

The three degree-two quotient pushforwards give a faithful coordinate system
on the odd-order subgroup of the genus-three Jacobian.  This script checks
their values on the exact four-point conorm and its three Frobenius translates.
"""
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve
import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
prior = json.loads((HERE / 'type_i_divisor_transfer.json').read_text())
assert prior['status'] == 'verified_type_i_point_level_conorm_norm'
assert prior['source_sha256'] == hashlib.sha256(
    (HERE / 'type_i_divisor_transfer.py').read_bytes()).hexdigest()

p, q, r = dt.p, dt.q, dt.r
E = EllipticCurve(dt.L, [dt.L(dt.gm.a), dt.L(dt.gm.b)])
G = E(dt.L(dt.gm.gx), dt.L(dt.gm.gy))
zero = E(0)
assert r.is_prime() and r % 2 == 1
assert r*G == zero and G != zero
target = -2*G
assert target != zero and r*target == zero

def quotient_image(point, i):
    """The p^i-conjugate quartic quotient, mapped to standard P-192."""
    assert i in (0, 1, 2)
    X = point[0]
    y = dt.check_model(point)[i]
    rho = dt.inf**(p**i)
    C = dt.c0**(p**i)
    D = dt.d0**(p**i)
    kappa = dt.kappa**(p**i)
    assert X != rho
    image = E(C+D/(X-rho), kappa*y/(X-rho)**2)
    assert image.curve() == E
    return image

def pushforward_vector(positive, negative):
    assert len(positive) == len(negative) == 2
    return [sum((quotient_image(P, i) for P in positive), zero)
            - sum((quotient_image(P, i) for P in negative), zero)
            for i in range(3)]

vectors = []
for j in range(3):
    power = q**j
    pos = [dt.frobenius(P, power) for P in dt.fiber_G]
    neg = [dt.frobenius(P, power) for P in dt.fiber_G2]
    vector = pushforward_vector(pos, neg)
    expected = [zero, zero, zero]
    expected[(-j) % 3] = target
    assert vector == expected, (j, vector, expected)
    vectors.append(vector)

norm_vector = [sum((vectors[j][i] for j in range(3)), zero)
               for i in range(3)]
assert norm_vector == [target, target, target]
assert all(v != zero and r*v == zero for v in norm_vector)
assert [dt.frobenius(P, q**3) for P in dt.fiber_G] == dt.fiber_G
assert [dt.frobenius(P, q**3) for P in dt.fiber_G2] == dt.fiber_G2

def encode_point(P):
    if P == zero:
        return None
    return [[str(c) for c in P[k].polynomial().list()] for k in (0, 1)]

out = {
    'status': 'verified_type_i_norm_nonzero_via_three_elliptic_quotients',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'divisor_transfer_source_sha256': prior['source_sha256'],
    'divisor_fiber_sha256': prior['fiber_sha256'],
    'base_prime_decimal': str(p),
    'subgroup_order_decimal': str(r),
    'three_quotient_pushforward_vectors':
        [[encode_point(P) for P in row] for row in vectors],
    'norm_pushforward_vector': [encode_point(P) for P in norm_vector],
    'expected_nonzero_component': encode_point(target),
    'norm_vector_equals_diagonal_minus_2G': True,
    'norm_class_fixed_by_q_frobenius': True,
    'norm_class_nonzero_and_exact_order_r': True,
    'reduction_model': 'three elliptic quotient pushforwards; faithful on odd torsion by the Klein-four push-pull identity',
    'boundary': 'r-primary isogeny-coordinate reduction and norm-order check; no generic plane-quartic divisor reduction or P-192 relation collector',
}
(HERE / 'type_i_norm_quotients.json').write_text(
    json.dumps(out, indent=2, sort_keys=True) + '\n')
print(json.dumps(out, indent=2, sort_keys=True))
