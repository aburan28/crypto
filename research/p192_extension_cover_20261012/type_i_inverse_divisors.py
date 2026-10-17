"""Point-level inverse-isogeny divisors for 16 P-192 Type I quotient classes.

This realizes each exact quotient vector by a degree-zero formal divisor on
the smooth plane quartic. Coefficients are not reduced to degree-three
Jacobian normal form.
"""
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, prod, vector
import type_i_divisor_transfer as dt
import verify_type_i_quotient_base as vb

HERE=Path(__file__).resolve().parent
p,r,L,F=dt.p,dt.r,dt.L,dt.gm.F
E=EllipticCurve(L,[L(dt.gm.a),L(dt.gm.b)])
O=E(0)
assert len(vb.base)==16 and r.is_prime() and r%2==1
half=(r+1)//2
plane=json.loads((HERE/'type_i_plane_quartic.json').read_text())
R=PolynomialRing(F,3,names=('Z0','Z1','Z2'))
Z=R.gens()
quartic=sum(F(int(c))*prod(Z[i]**exp[i] for i in range(3))
            for key,c in plane['quartic_coefficients'].items()
            for exp in [tuple(map(int,key.split(',')))])

def quotient(P,i):
    X=P[0]
    y=dt.check_model(P)[i]
    rho,C,D,kappa=[a**(p**i) for a in (dt.inf,dt.c0,dt.d0,dt.kappa)]
    assert X!=rho
    return E(C+D/(X-rho),kappa*y/(X-rho)**2)

def frob(P,i):
    Q=tuple(a**(p**i) for a in P)
    dt.check_model(Q)
    return Q

def on_plane(P):
    ys=dt.check_model(P)
    assert all(y!=0 for y in ys)
    coord=dt.Vandermonde.solve_right(vector(L,[1/y for y in ys]))
    assert quartic(*coord)==0

def encode(P):
    return [[str(c) for c in a.polynomial().list()] for a in P]

def vec(positive,negative):
    return tuple(sum((quotient(P,i) for P in positive),O)-
                 sum((quotient(P,i) for P in negative),O)
                 for i in range(3))

classes=[]
for index,target in enumerate(vb.base):
    assert target[0]==target[1]==target[2]
    q=target[0]
    assert q!=O and r*q==O and all(c**p==c for c in (q[0],q[1]))
    fq=dt.fiber(q)
    f2q=dt.fiber(2*q)
    assert len(fq)==len(f2q)==2
    for P in fq+f2q:
        on_plane(P)
    component_vectors=[]
    positive=[]
    negative=[]
    for i in range(3):
        pos=[frob(P,i) for P in f2q]
        neg=[frob(P,i) for P in fq]
        for P in pos+neg:
            on_plane(P)
        component=vec(pos,neg)
        expected=tuple(2*q if j==i else O for j in range(3))
        assert component==expected,(index,i,component,expected)
        component_vectors.append(component)
        positive.extend(pos)
        negative.extend(neg)
    unscaled=vec(positive,negative)
    assert unscaled==(2*q,2*q,2*q)
    recovered=tuple(half*Q for Q in unscaled)
    assert recovered==target
    classes.append({'base_index':index,
                    'fiber_q':list(map(encode,fq)),
                    'fiber_2q':list(map(encode,f2q)),
                    'formal_divisor':'half_mod_r * sum_i (Frob_p^i fiber(2q) - Frob_p^i fiber(q))',
                    'positive_points':6,'negative_points':6,
                    'each_component_pushforward_twice_q':True,
                    'scaled_quotient_vector_matches_base':True,
                    'all_points_on_smooth_plane_quartic':True})
    print(f'inverse divisor {index+1}/16 verified',flush=True)

out={'status':'verified_type_i_inverse_isogeny_formal_divisors',
     'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'base_receipt_sha256':hashlib.sha256((HERE/'type_i_quotient_base_probe.json').read_bytes()).hexdigest(),
     'plane_receipt_sha256':hashlib.sha256((HERE/'type_i_plane_quartic.json').read_bytes()).hexdigest(),
     'p':str(p),'r':str(r),'inverse_of_two_mod_r':str(half),
     'classes':classes,'classes_verified':len(classes),
     'formal_points_per_class':12,
     'reduced_degree_three_representatives':None,
     'boundary':'point-level degree-zero inverse-isogeny divisors with exact quotient images; generic Jacobian reduction remains open'}
(HERE/'type_i_inverse_divisors.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='classes'},indent=2,sort_keys=True))
