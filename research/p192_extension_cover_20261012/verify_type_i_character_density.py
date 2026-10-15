"""Independent polynomial and integer replay of the P-192 density inputs."""
import contextlib
import hashlib
import io
import json
from fractions import Fraction
from math import isqrt
from pathlib import Path
import random

with contextlib.redirect_stdout(io.StringIO()):
    import type_i_divisor_transfer as dt
from sage.all import PolynomialRing, prod

HERE=Path(__file__).resolve().parent
j=json.loads((HERE/'type_i_character_density.json').read_text())
assert j['source_sha256']==hashlib.sha256((HERE/'type_i_character_density.py').read_bytes()).hexdigest()
p,r,F,L=dt.p,dt.r,dt.gm.F,dt.L
assert j['p']==str(p) and j['r']==str(r)
RL=PolynomialRing(L,'x')
x=RL.gen()
num=RL(1)
den=RL(1)
for i in range(3):
    rho=dt.inf**(p**i)
    C=dt.c0**(p**i)
    D=dt.d0**(p**i)
    roots=[v**(p**(i+1)) for v in dt.quartic_roots]
    num*=prod((rho-v)*(x-C)+D for v in roots)
    den*=(x-C)**4
assert all(c**p==c for c in num.list()+den.list())
R=PolynomialRing(F,'x')
X=R.gen()
to_F=lambda f:R([F(a.polynomial()[0]) for a in f.list()])
numer,denom=to_F(num),to_F(den)
e=X**3+F(dt.gm.a)*X+F(dt.gm.b)
assert e.is_irreducible()
f=R([F(int(a)) for a in j['sextic_coefficients_decimal']])
assert numer==denom.leading_coefficient()*e*f
assert len(list(f.factor()))==2
assert [list(map(str,g.list())) for g,k in f.factor()]==j['sextic_factor_coefficients_decimal']
assert f.gcd(f.derivative())==1 and f.gcd(e)==1
assert [(g.degree(),k) for g,k in denom.factor()]==[(3,4)]
a,b=F(dt.gm.a),F(dt.gm.b)
u=X**4-2*a*X**2-8*b*X+a**2
v=4*e
pulled=R(0)
for i,c in enumerate(f.list()):
    pulled+=c*u**i*v**(6-i)
assert pulled.degree()==24
assert pulled.gcd(pulled.derivative())==1
assert pulled.gcd(f)==1 and pulled.gcd(e)==1
guided=json.loads((HERE/'type_i_geom_character_collector.json').read_text())
rng=random.Random(guided['seed'])
for _ in range(20):
    P=rng.randrange(1,int(r))*dt.G
    x0=F(P[0])
    assert v(x0)!=0
    assert u(x0)/v(x0)==F((2*P)[0])
    assert pulled(x0)==v(x0)**6*f(F((2*P)[0]))

U=isqrt(int(p))+1
t=int(p+1-r)
assert U==int(j['ceil_sqrt_p']) and t==int(j['trace'])
chi=f.leading_coefficient()**((p-1)//2)
delta=1 if chi==1 else 0
assert delta==j['identity_accepted']
err=90*U+3*abs(t)
lower=(int(r)-err+3)//4-delta
assert lower==int(j['accepted_nonzero_scalar_lower_bound'])
assert Fraction(816,lower)==Fraction(j['conditioned_hit_probability_upper_bound'])
assert 816*(1<<180)<lower and j['conditioned_hit_probability_lt_two_pow_minus_180']
assert j['frozen_accepted']==64 and j['frozen_first_rejects']==128 and j['frozen_second_rejects']==67
out={'status':'verified_type_i_character_density_inputs',
     'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'producer_sha256':j['source_sha256'],
     'independent_norm_polynomial_reconstruction':True,
     'branch_polynomials_squarefree_and_disjoint':True,
     'doubling_map_points_replayed':20,
     'accepted_nonzero_lower_bound':str(lower),
     'conditioned_hit_probability_upper_bound':str(Fraction(816,lower)),
     'proof_boundary':'Hasse--Weil and Riemann--Hurwitz applied in the written proof; Lean checks only the arithmetic implication'}
(HERE/'type_i_character_density_verification.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
