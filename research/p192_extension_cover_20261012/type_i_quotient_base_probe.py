"""Bounded P-192 Type I base and known-scalar three-sum probe.

All arithmetic is in the three exact elliptic quotient coordinates of the
descended genus-three cover. The receipt labels these as quotient-coordinate
classes; a generic Jacobian reduction/inverse-isogeny implementation remains
a separate gate.
"""
import hashlib
import json
import random
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, matrix, prod, vector
import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
p, q, r, L = dt.p, dt.q, dt.r, dt.L
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
G = E(L(dt.gm.gx), L(dt.gm.gy))
O = E(0)
assert G != O and r.is_prime() and r*G == O
trace = p + 1 - r
s0, s1 = Integer(2), trace
for _ in range(2, 7):
    s0, s1 = s1, trace*s1 - p*s0
extension_order = p**6 + 1 - s1
assert extension_order % r == 0
cofactor = extension_order // r
assert cofactor % r != 0

def cover_point(x):
    """One F_p-descended affine point above x, with the norm-component sign."""
    X = L(x)
    U = (X-dt.zeta)*(X-dt.beta)
    y = U**((q+1)//2)
    assert y*y == prod(X-root for root in dt.quartic_roots)
    norm_roots = prod(X-root**(p**i)
                      for root in (dt.zeta, dt.beta) for i in range(3))
    if y*y**p*y**(p**2) != norm_roots:
        y = -y
    assert y*y**p*y**(p**2) == norm_roots
    values = vector(L, [y, y**p, y**(p**2)])
    coeffs = dt.Vandermonde.solve_right(values)
    point = (X, coeffs[0], coeffs[1], coeffs[2])
    assert all(c**p == c for c in point)
    dt.check_model(point)
    return point

def quotient_image(point, i):
    X = point[0]
    y = dt.check_model(point)[i]
    rho, C, D, kappa = [v**(p**i) for v in
                         (dt.inf, dt.c0, dt.d0, dt.kappa)]
    assert X != rho
    return E(C+D/(X-rho), kappa*y/(X-rho)**2)

def vector_key(points):
    return tuple(str(P) for P in points)

def add_vectors(left, right):
    return tuple(a+b for a, b in zip(left, right))

def sub_vectors(left, right):
    return tuple(a-b for a, b in zip(left, right))

seed = 0x19220261013
rng = random.Random(seed)
xs = [Integer(rng.randrange(int(p))) for _ in range(256)]
pilot_digest = hashlib.sha256(','.join(map(str, xs)).encode()).hexdigest()
assert pilot_digest == '92166e6477a7738fdae1381531c864499523d404f0ac6e3ba3cadb0cdc8cfaf8'
anchor = cover_point(xs[0])
anchor_images = tuple(quotient_image(anchor, i) for i in range(3))

base = []
selection = []
start = time.monotonic_ns()
for index, x in enumerate(xs[1:], 1):
    if len(base) >= 16:
        break
    point = cover_point(x)
    images = tuple(quotient_image(point, i) for i in range(3))
    delta = sub_vectors(images, anchor_images)
    projected = tuple(cofactor*P for P in delta)
    assert all(r*P == O for P in projected)
    if projected == (O, O, O) or vector_key(projected) in [vector_key(v) for v in base]:
        continue
    base.append(projected)
    selection.append({'pilot_index': index, 'x': str(x),
                      'diagonal_after_projection': projected[0] == projected[1] == projected[2]})
    print(f'base {len(base)}/16, pilot index {index}', flush=True)
base_ns = time.monotonic_ns() - start
assert len(base) == 16

start = time.monotonic_ns()
pairs = {}
for i in range(16):
    for j in range(i, 16):
        pairs.setdefault(vector_key(add_vectors(base[i], base[j])), []).append((i, j))
assert sum(map(len, pairs.values())) == 136
pair_ns = time.monotonic_ns() - start

N = (-2*G, -2*G, -2*G)
assert all(P != O and r*P == O for P in N)
query_seed = 0x19220261019
query_rng = random.Random(query_seed)
scalars = [Integer(query_rng.randrange(1, int(r))) for _ in range(64)]
query_digest = hashlib.sha256(','.join(map(str, scalars)).encode()).hexdigest()
hits = []
start = time.monotonic_ns()
for query_index, scalar in enumerate(scalars):
    target = tuple(scalar*P for P in N)
    for k in range(16):
        remainder = sub_vectors(target, base[k])
        for i, j in pairs.get(vector_key(remainder), []):
            assert add_vectors(add_vectors(base[i], base[j]), base[k]) == target
            hits.append({'query_index': query_index, 'scalar': str(scalar),
                         'indices': [i, j, k]})
    if (query_index+1) % 8 == 0:
        print(f'checked {query_index+1}/64 known-scalar quotient queries', flush=True)
query_ns = time.monotonic_ns() - start
rows = []
for hit in hits:
    row = [0]*16
    for index in hit['indices']:
        row[index] += 1
    rows.append(row)
rank = int(matrix(GF(r), rows).rank()) if rows else 0

out = {
    'status': 'bounded_type_i_quotient_coordinate_base_probe',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'divisor_transfer_source_sha256': hashlib.sha256((HERE/'type_i_divisor_transfer.py').read_bytes()).hexdigest(),
    'proposal_id': 'Q-P192-TypeI-Fp-quotient-16',
    'candidate_id': None,
    'p': str(p), 'r': str(r),
    'extension_elliptic_order': str(extension_order),
    'candidate_cofactor': str(cofactor),
    'pilot_seed': seed,
    'pilot_x_sha256_decimal_csv': pilot_digest,
    'anchor_pilot_index': 0,
    'base_policy': 'first 16 distinct nonzero projected three-quotient vectors from frozen Fp abscissae after index 0 anchor',
    'base_selection': selection,
    'base_classes': 16,
    'diagonal_projected_classes': sum(v['diagonal_after_projection'] for v in selection),
    'pair_index_entries': 136,
    'query_policy': '64 independent frozen uniform nonzero scalars times the checked Type I three-step norm quotient vector (-2G,-2G,-2G)',
    'query_seed': query_seed,
    'query_scalar_sha256_decimal_csv': query_digest,
    'ordinary_known_scalar_queries': 64,
    'verified_quotient_three_sum_matches': len(hits),
    'relation_witnesses': hits,
    'quotient_incidence_rank': rank,
    'base_build_ns': base_ns,
    'pair_index_ns': pair_ns,
    'query_and_match_ns': query_ns,
    'boundary': 'exact elliptic quotient-coordinate r-vectors; generic genus-three Jacobian representative, inverse-isogeny, and complete relation collector remain separate',
}
(HERE/'type_i_quotient_base_probe.json').write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
print(json.dumps(out, indent=2, sort_keys=True))
