"""Independent replay of Type I point-level inverse-isogeny divisors."""
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, prod, vector
import type_i_divisor_transfer as dt
import verify_type_i_quotient_base as vb

HERE=Path(__file__).resolve().parent
j=json.loads((HERE/'type_i_inverse_divisors.json').read_text())
assert j['source_sha256']==hashlib.sha256((HERE/'type_i_inverse_divisors.py').read_bytes()).hexdigest()
assert j['base_receipt_sha256']==hashlib.sha256((HERE/'type_i_quotient_base_probe.json').read_bytes()).hexdigest()
assert j['plane_receipt_sha256']==hashlib.sha256((HERE/'type_i_plane_quartic.json').read_bytes()).hexdigest()
p,r,L,F=dt.p,dt.r,dt.L,dt.gm.F
assert j['p']==str(p) and j['r']==str(r)
half=int(j['inverse_of_two_mod_r'])
assert (2*half)%int(r)==1
E=EllipticCurve(L,[L(dt.gm.a),L(dt.gm.b)])
O=E(0)

plane=json.loads((HERE/'type_i_plane_quartic.json').read_text())
R=PolynomialRing(F,3,names=('Z0','Z1','Z2'))
Z=R.gens()
quartic=sum(F(int(c))*prod(Z[i]**e[i] for i in range(3))
            for key,c in plane['quartic_coefficients'].items()
            for e in [tuple(map(int,key.split(',')))])

def decode(point):
    return tuple(L([int(c) for c in a]) for a in point)

def image(P,k):
    X=P[0]
    values=dt.check_model(P)
    rho=dt.inf**(p**k)
    C=dt.c0**(p**k)
    D=dt.d0**(p**k)
    scale=dt.kappa**(p**k)
    assert X!=rho
    return E(C+D/(X-rho),scale*values[k]/(X-rho)**2)

def check_plane(P):
    y=dt.check_model(P)
    coord=dt.Vandermonde.solve_right(vector(L,[1/a for a in y]))
    assert quartic(*coord)==0

assert len(j['classes'])==len(vb.base)==16
for item in j['classes']:
    i=item['base_index']
    q=vb.base[i][0]
    assert vb.base[i]==(q,q,q) and r*q==O and q!=O
    neg=[decode(v) for v in item['fiber_q']]
    pos=[decode(v) for v in item['fiber_2q']]
    assert len(pos)==len(neg)==2
    assert all(image(P,0)==q for P in neg)
    assert all(image(P,0)==2*q for P in pos)
    totals=[O,O,O]
    for shift in range(3):
        positive=[tuple(a**(p**shift) for a in P) for P in pos]
        negative=[tuple(a**(p**shift) for a in P) for P in neg]
        for P in positive+negative:
            check_plane(P)
        row=[sum((image(P,k) for P in positive),O)-
             sum((image(P,k) for P in negative),O) for k in range(3)]
        assert row==[2*q if k==shift else O for k in range(3)]
        totals=[totals[k]+row[k] for k in range(3)]
    assert tuple(half*P for P in totals)==vb.base[i]
    assert item['positive_points']==item['negative_points']==6

out={'status':'verified_type_i_inverse_isogeny_formal_divisors',
     'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'producer_sha256':j['source_sha256'],
     'classes_replayed':16,'formal_points_per_class':12,
     'all_points_on_smooth_plane_quartic':True,
     'all_scaled_quotient_vectors_match':True}
(HERE/'type_i_inverse_divisors_verification.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
