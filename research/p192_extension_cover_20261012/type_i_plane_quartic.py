"""Exact descended Ciani plane quartic for the P-192 Type I cover.

The coordinates are the Frobenius-descent basis of the three elliptic
quotient differentials dx/y_i. This checks the polynomial and known points;
Jacobian class reduction is a separate arithmetic operation.
"""
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, ProjectiveSpace, prod, vector
from sage.schemes.curves.projective_curve import IntegralProjectivePlaneCurve_finite_field
import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
p, F, L = dt.p, dt.gm.F, dt.L
R = PolynomialRing(L, 3, names=('u','v','w'))
u,v,w = R.gens()
z,beta = dt.zeta,dt.beta
s = [z**(p**i)+beta**(p**i) for i in range(3)]
n = [z**(p**i)*beta**(p**i) for i in range(3)]
s0,s1,s2 = s
n0,n1,n2 = n
L0 = (n0-n2)*u**2 + (n1-n0)*v**2 + (n2-n1)*w**2
L1 = (s2-s0)*u**2 + (s0-s1)*v**2 + (s1-s2)*w**2
L2 = (s0*n2-n0*s2)*u**2 + (s1*n0-n1*s0)*v**2 + (n1*s2-s1*n2)*w**2
ciani = L0**2-L1*L2
assert ciani != 0 and ciani.total_degree() == 4

# Direct Sylvester resultant of t0*u^2-t1*w^2 and t0*v^2-t2*w^2.
a2,a1,a0 = u**2-w**2,-s0*u**2+s1*w**2,n0*u**2-n1*w**2
b2,b1,b0 = v**2-w**2,-s0*v**2+s2*w**2,n0*v**2-n2*w**2
resultant = (a2*b0-a0*b2)**2-(a2*b1-a1*b2)*(a1*b0-a0*b1)
assert resultant == w**4*ciani

T = PolynomialRing(L, 3, names=('Z0','Z1','Z2'))
Z0,Z1,Z2 = T.gens()
forms = [Z0+z**(p**i)*Z1+z**(2*p**i)*Z2 for i in range(3)]
descended_L = ciani(*forms)
assert descended_L != 0 and descended_L.total_degree() == 4
assert all(c**p == c for c in descended_L.dict().values())
TF = PolynomialRing(F, 3, names=('Z0','Z1','Z2'))
quartic = TF({exp:F(c.polynomial()[0]) for exp,c in descended_L.dict().items()})
assert T(quartic) == descended_L
P2 = ProjectiveSpace(F,2,names=('Z0','Z1','Z2'))
Z = P2.coordinate_ring().gens()
plane_equation = sum(F(c)*prod(Z[i]**exp[i] for i in range(3))
                     for exp,c in quartic.dict().items())
# The direct constructor avoids a Sage factorization precheck. The subsequent
# smoothness check independently certifies integrality of this plane quartic.
curve = IntegralProjectivePlaneCurve_finite_field(P2,plane_equation)
assert curve.degree() == 4 and curve.is_smooth()
jacobian_error = None
try:
    curve.jacobian(model='km_small')
except Exception as exc:
    jacobian_error = {'type':type(exc).__name__,
                      'message':str(exc).splitlines()[0]}

point_checks = []
for index, P in enumerate(dt.fiber_G+dt.fiber_G2):
    ys = dt.check_model(P)
    assert all(y != 0 for y in ys)
    inv = vector(L,[1/y for y in ys])
    coord = dt.Vandermonde.solve_right(inv)
    assert quartic(*coord) == 0
    point_checks.append({'source':'generator_fiber','index':index,
                         'coordinate_encoding':[[str(c) for c in a.polynomial().list()] for a in coord]})

# The first 16 F_p base points and the anchor provide rational plane points.
base_receipt = json.loads((HERE/'type_i_quotient_base_probe.json').read_text())
def cover_point(x):
    X = L(x)
    U = (X-dt.zeta)*(X-dt.beta)
    y = U**((dt.q+1)//2)
    roots = prod(X-root**(p**j) for root in (dt.zeta,dt.beta) for j in range(3))
    if y*y**p*y**(p**2) != roots:
        y = -y
    assert y*y**p*y**(p**2) == roots
    coord = dt.Vandermonde.solve_right(vector(L,[y,y**p,y**(p**2)]))
    P = (X,*coord)
    assert all(c**p == c for c in P)
    dt.check_model(P)
    return P
xs = [int(base_receipt['base_selection'][i]['x']) for i in range(16)]
for index,x in enumerate(xs):
    P = cover_point(x)
    ys = dt.check_model(P)
    assert all(y != 0 for y in ys)
    coord = dt.Vandermonde.solve_right(vector(L,[1/y for y in ys]))
    assert all(c**p == c for c in coord)
    assert quartic(*coord) == 0
    point_checks.append({'source':'base','index':index,
                         'coordinates_Fp_decimal':[str(c.polynomial()[0]) for c in coord]})

out = {'status':'verified_type_i_descended_plane_quartic_points',
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'p':str(p), 'equation':'L0^2-L1*L2 in conjugate differential coordinates, then Frobenius-linear descent',
       'quartic_coefficients':{','.join(map(str,e)):str(c) for e,c in sorted(quartic.dict().items())},
       'quartic_degree':int(quartic.total_degree()),
       'smooth_plane_quartic_verified':True,
       'geometric_genus_from_smooth_plane_degree':3,
       'resultant_factorization_verified':True,
       'coefficients_in_Fp_verified':True,
       'known_points_verified':point_checks,
       'sage_km_small_jacobian_error':jacobian_error,
       'reduced_jacobian_classes':None}
(HERE/'type_i_plane_quartic.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('quartic_coefficients','known_points_verified')},indent=2,sort_keys=True))
