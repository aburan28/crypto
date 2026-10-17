"""Exact primitive genus-three equation and bounded Sage reduction attempt."""
import hashlib
import json
from pathlib import Path

from sage.all import FunctionField, PolynomialRing, prod
import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
p, F, L = dt.p, dt.gm.F, dt.L
RL = PolynomialRing(L,'X')
X = RL.gen()
f = [prod(X-z**(p**i) for z in dt.quartic_roots) for i in range(3)]
A = sum(f)
B = f[0]*f[1]+f[0]*f[2]+f[1]*f[2]
H = prod(X-z**(p**i) for z in (dt.zeta,dt.beta) for i in range(3))
assert [A.degree(),B.degree(),H.degree()] == [4,8,6]
assert all(c**p == c for poly in (A,B,H) for c in poly.list())

def descend(poly):
    RF = PolynomialRing(F,'X')
    return RF([F(c.polynomial()[0]) for c in poly.list()])

a,b,h = map(descend,(A,B,H))
assert all(RL(poly) == original for poly,original in zip((a,b,h),(A,B,H)))
for point in dt.fiber_G + dt.fiber_G2:
    x = point[0]
    ys = dt.check_model(point)
    s = sum(ys)
    assert s**4-2*A(x)*s**2-8*H(x)*s+A(x)**2-4*B(x) == 0

K = FunctionField(F,'x')
x = K.gen()
S = PolynomialRing(K,'s').gen()
quartic = S**4-2*K(a(x))*S**2-8*K(h(x))*S+K(a(x))**2-4*K(b(x))
failure = None
try:
    field = K.extension(quartic,'s')
    genus = int(field.genus())
    jacobian = field.jacobian(model='hess')
    status = 'jacobian_constructed'
except Exception as exc:
    status = 'jacobian_construction_failed'
    failure = {'type':type(exc).__name__,'message':str(exc)}
    genus = None

out = {'status':status,
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'p':str(p), 'primitive_degree':4,
       'A_degree':int(A.degree()),'B_degree':int(B.degree()),'H_degree':int(H.degree()),
       'all_coefficients_in_Fp':True,
       'fiber_points_checked':len(dt.fiber_G)+len(dt.fiber_G2),
       'quartic_equation':'s^4 - 2*A(x)*s^2 - 8*H(x)*s + A(x)^2 - 4*B(x)',
       'A_coefficients_decimal':list(map(str,a.list())),
       'B_coefficients_decimal':list(map(str,b.list())),
       'H_coefficients_decimal':list(map(str,h.list())),
       'sage_jacobian_model':'hess', 'sage_genus':genus,
       'failure':failure,
       'reduced_P192_jacobian_representatives':None}
(HERE/'type_i_jacobian_model_probe.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in out.items() if not k.endswith('coefficients_decimal')},indent=2,sort_keys=True))
