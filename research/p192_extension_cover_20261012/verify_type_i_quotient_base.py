"""Independent exhaustive replay of the P-192 Type I quotient-base probe.

This verifier rebuilds the cover points and projected quotient coordinates,
then enumerates all 816 unordered triples directly instead of using the
producer's pair-index lookup.
"""
import hashlib
import itertools
import json
import random
from pathlib import Path

from sage.all import EllipticCurve, Integer, matrix, prod, vector
import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
receipt = json.loads((HERE/'type_i_quotient_base_probe.json').read_text())
assert receipt['status'] == 'bounded_type_i_quotient_coordinate_base_probe'
assert receipt['source_sha256'] == hashlib.sha256((HERE/'type_i_quotient_base_probe.py').read_bytes()).hexdigest()
p, r, L = dt.p, dt.r, dt.L
assert str(p) == receipt['p'] and str(r) == receipt['r']
E = EllipticCurve(L, [L(dt.gm.a), L(dt.gm.b)])
G = E(L(dt.gm.gx), L(dt.gm.gy))
O = E(0)
assert r*G == O and G != O
h = Integer(receipt['candidate_cofactor'])
assert h*r == Integer(receipt['extension_elliptic_order']) and h % r != 0

rng = random.Random(receipt['pilot_seed'])
xs = [Integer(rng.randrange(int(p))) for _ in range(256)]
assert hashlib.sha256(','.join(map(str, xs)).encode()).hexdigest() == receipt['pilot_x_sha256_decimal_csv']

def point_from_x(x):
    X = L(x)
    U = (X-dt.zeta)*(X-dt.beta)
    y = U**((dt.q+1)//2)
    six_roots = prod((X-root**(p**j)) for root in (dt.zeta, dt.beta) for j in range(3))
    if y*y**p*y**(p*p) == -six_roots:
        y = -y
    assert y*y**p*y**(p*p) == six_roots
    assert y*y == prod(X-z for z in dt.quartic_roots)
    coord = dt.Vandermonde.solve_right(vector(L,[y,y**p,y**(p*p)]))
    P = (X, *coord)
    assert all(c**p == c for c in P)
    dt.check_model(P)
    return P

def image(P, j):
    X = P[0]
    yj = dt.check_model(P)[j]
    rho = dt.inf**(p**j)
    C = dt.c0**(p**j)
    D = dt.d0**(p**j)
    kappa = dt.kappa**(p**j)
    return E(C + D/(X-rho), kappa*yj/(X-rho)**2)

anchor = point_from_x(xs[receipt['anchor_pilot_index']])
anchor_images = [image(anchor,j) for j in range(3)]
base = []
for item in receipt['base_selection']:
    idx = item['pilot_index']
    assert str(xs[idx]) == item['x']
    P = point_from_x(xs[idx])
    projected = tuple(h*(image(P,j)-anchor_images[j]) for j in range(3))
    assert all(r*Q == O for Q in projected)
    assert projected != (O,O,O)
    assert projected not in base
    assert (projected[0] == projected[1] == projected[2]) == item['diagonal_after_projection']
    base.append(projected)
assert len(base) == receipt['base_classes'] == 16

triples = list(itertools.combinations_with_replacement(range(16),3))
assert len(triples) == 816
triple_values = [tuple(base[i][j]+base[k][j]+base[l][j] for j in range(3))
                 for i,k,l in triples]
query_rng = random.Random(receipt['query_seed'])
scalars = [Integer(query_rng.randrange(1,int(r))) for _ in range(64)]
assert hashlib.sha256(','.join(map(str,scalars)).encode()).hexdigest() == receipt['query_scalar_sha256_decimal_csv']
N = (-2*G,-2*G,-2*G)
matches = []
for q,a in enumerate(scalars):
    target = tuple(a*Q for Q in N)
    for indices, candidate in zip(triples,triple_values):
        if candidate == target:
            matches.append({'query_index':q,'scalar':str(a),'indices':list(indices)})
assert len(matches) == receipt['verified_quotient_three_sum_matches']
assert matches == receipt['relation_witnesses']

out = {
    'status':'verified_type_i_quotient_base_and_exhaustive_triples',
    'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'producer_sha256':receipt['source_sha256'],
    'base_classes_replayed':len(base),
    'all_base_vectors_diagonal':all(v[0]==v[1]==v[2] for v in base),
    'all_base_vectors_nonzero_r_torsion':True,
    'unordered_triples_enumerated':len(triples),
    'known_scalar_queries_replayed':len(scalars),
    'exact_matches':len(matches),
}
(HERE/'type_i_quotient_base_verification.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
