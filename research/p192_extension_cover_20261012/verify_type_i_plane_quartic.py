"""Independent point and smoothness replay of the descended plane quartic."""
import hashlib
import json
from pathlib import Path

from sage.all import GF, ProjectiveSpace, prod, vector
from sage.schemes.curves.projective_curve import IntegralProjectivePlaneCurve_finite_field
import type_i_divisor_transfer as dt

HERE = Path(__file__).resolve().parent
j = json.loads((HERE/'type_i_plane_quartic.json').read_text())
assert j['source_sha256'] == hashlib.sha256((HERE/'type_i_plane_quartic.py').read_bytes()).hexdigest()
p,F,L = dt.p,dt.gm.F,dt.L
assert j['p'] == str(p)
P2 = ProjectiveSpace(F,2,names=('Z0','Z1','Z2'))
Z = P2.coordinate_ring().gens()
f = sum(F(int(c))*prod(Z[i]**exp[i] for i in range(3))
        for key,c in j['quartic_coefficients'].items()
        for exp in [tuple(map(int,key.split(',')))])
assert f.total_degree() == 4
C = IntegralProjectivePlaneCurve_finite_field(P2,f)
assert C.is_smooth()

points = j['known_points_verified']
assert len(points) == 20
for item in points:
    if item['source'] == 'generator_fiber':
        P = (dt.fiber_G+dt.fiber_G2)[item['index']]
        y = dt.check_model(P)
        coord = dt.Vandermonde.solve_right(vector(L,[1/a for a in y]))
        encoded = [[str(c) for c in a.polynomial().list()] for a in coord]
        assert encoded == item['coordinate_encoding']
        assert f(*coord) == 0
    else:
        coord = [F(int(v)) for v in item['coordinates_Fp_decimal']]
        assert f(*coord) == 0
        assert any(coord)
        i = item['index']
        x = int(json.loads((HERE/'type_i_quotient_base_probe.json').read_text())['base_selection'][i]['x'])
        U = L(x)
        y = ((U-dt.zeta)*(U-dt.beta))**((dt.q+1)//2)
        H = prod(U-z**(p**k) for z in (dt.zeta,dt.beta) for k in range(3))
        if y*y**p*y**(p**2) != H:
            y = -y
        assert y*y**p*y**(p**2) == H
        Y = (y,y**p,y**(p**2))
        inverse = [sum(L(coord[k])*dt.zeta**(p**m*k) for k in range(3))
                   for m in range(3)]
        assert all(inverse[m]*Y[m] == 1 for m in range(3))

out = {'status':'verified_type_i_plane_quartic_smoothness_and_points',
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'producer_sha256':j['source_sha256'],
       'smooth':True,'degree':4,'genus_by_plane_formula':3,
       'generator_fiber_points_replayed':4,'rational_base_points_replayed':16}
(HERE/'type_i_plane_quartic_verification.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
