"""Explicit point-level conorm and Frobenius norm for the P-192 Type I cover.

The output is a formal degree-zero divisor on the descended affine model.
Genus-three Jacobian reduction is a separate arithmetic gate.
"""
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, matrix, prod, vector
import type_i_generator_map as gm

HERE = Path(__file__).resolve().parent
prior = json.loads((HERE / 'type_i_generator_map.json').read_text())
assert prior['status'] == 'verified_type_i_quartic_generator_coordinates'
L, p = gm.L, gm.p
zeta, beta = gm.z, gm.be
quartic_roots = gm.rl
inf, c0, d0, kappa = gm.inf, gm.c0, gm.d0, gm.kappa
r = gm.Integer('FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831', 16)
q = p**2
E = EllipticCurve(gm.F, [gm.F(gm.a), gm.F(gm.b)])
G = E(gm.F(gm.gx), gm.F(gm.gy))
assert r * G == E(0) and G != E(0)
G2 = 2*G
assert G2 != E(0)

conjugates = [zeta**(p**i) for i in range(3)]
Vandermonde = matrix(L, [[L(1), c, c**2] for c in conjugates])
assert Vandermonde.det() != 0

def check_model(point):
    X, Y0, Y1, Y2 = point
    values = [Y0+c*Y1+c**2*Y2 for c in conjugates]
    for i in range(3):
        rhs = prod(X-v**(p**i) for v in quartic_roots)
        assert values[i]**2 == rhs
    norm_roots = prod(X-v**(p**i) for v in (zeta,beta) for i in range(3))
    assert prod(values) == norm_roots
    return values

def fiber(P):
    x, y = L(P[0]), L(P[1])
    assert x != c0
    X = inf+d0/(x-c0)
    y0 = y*(X-inf)**2/kappa
    assert y0**2 == prod(X-v for v in quartic_roots)
    rhs1 = prod(X-v**p for v in quartic_roots)
    assert rhs1.is_square() and y0 != 0
    root1 = rhs1.sqrt()
    norm_roots = prod(X-v**(p**i) for v in (zeta,beta) for i in range(3))
    points = []
    for sign in (1,-1):
        y1 = sign*root1
        assert y1 != 0
        y2 = norm_roots/(y0*y1)
        coords = Vandermonde.solve_right(vector(L, [y0,y1,y2]))
        point = (X,coords[0],coords[1],coords[2])
        check_model(point)
        points.append(point)
    assert points[0] != points[1]
    return points

print('building two rational fibers and six Frobenius pairs', flush=True)
fiber_G, fiber_G2 = fiber(G), fiber(G2)
frobenius = lambda point, power: tuple(c**power for c in point)
for initial in fiber_G+fiber_G2:
    orbit = [frobenius(initial, q**i) for i in range(3)]
    for point in orbit:
        check_model(point)
    assert frobenius(orbit[-1], q) == initial

def encode(point):
    return [[str(c) for c in coordinate.polynomial().list()]
            for coordinate in point]

positive = [encode(P) for P in fiber_G]
negative = [encode(P) for P in fiber_G2]
raw = json.dumps({'G':positive,'2G':negative},sort_keys=True,separators=(',',':'))
out = {
    'status': 'verified_type_i_point_level_conorm_norm',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'generator_map_source_sha256': prior['source_sha256'],
    'group_relation': 'G - 2G = -G, order r on P-192',
    'fiber_G_Fp6': positive,
    'fiber_2G_Fp6': negative,
    'fiber_sha256': hashlib.sha256(raw.encode()).hexdigest(),
    'fiber_cardinality_each': 2,
    'all_fiber_points_on_descended_model': True,
    'frobenius_q_orbits_return_after_three_steps': True,
    'formal_conorm': 'Pplus(G)+Pminus(G)-Pplus(2G)-Pminus(2G)',
    'formal_norm_to_Fp2': 'sum of q^i-Frobenius translates of the formal conorm for i=0,1,2; q=p^2',
    'boundary': 'exact point-level divisor and Frobenius orbit; genus-three Jacobian class reduction, nonzero norm class, and subgroup order need a separate checker',
}
(HERE / 'type_i_divisor_transfer.json').write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('fiber_G_Fp6','fiber_2G_Fp6')}, indent=2, sort_keys=True))
