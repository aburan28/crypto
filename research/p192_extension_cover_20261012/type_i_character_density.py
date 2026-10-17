"""Exact inputs for a P-192 Type I character-filter density bound.

The finite-field identities, branch polynomial conditions, and frozen query
replay are checked here. The genus and point-count implications use the
written Riemann--Hurwitz and Hasse--Weil argument in BRANCH_CONTINUATION.md.
"""
import contextlib
import hashlib
import io
import json
import random
import time
from fractions import Fraction
from math import isqrt
from pathlib import Path

with contextlib.redirect_stdout(io.StringIO()):
    import type_i_divisor_transfer as dt
from sage.all import PolynomialRing, prod

HERE=Path(__file__).resolve().parent
p,r,F,L=dt.p,dt.r,dt.gm.F,dt.L
start=time.monotonic_ns()
RL=PolynomialRing(L,'x')
x=RL.gen()
rhs=[]
for i in range(3):
    rho,C,D=[a**(p**i) for a in (dt.inf,dt.c0,dt.d0)]
    X=rho+D/(x-C)
    rhs.append(prod(X-z**(p**(i+1)) for z in dt.quartic_roots))
norm=prod(rhs)
assert all(c**p==c for c in norm.numerator().list()+norm.denominator().list())
R=PolynomialRing(F,'x')
X=R.gen()
descend=lambda poly:R([F(c.polynomial()[0]) for c in poly.list()])
num=descend(norm.numerator())
den=descend(norm.denominator())
e=X**3+F(dt.gm.a)*X+F(dt.gm.b)
assert e.is_irreducible()
q,rem=num.quo_rem(e)
assert rem==0 and q.degree()==6
den_factors=list(den.factor())
assert len(den_factors)==1 and den_factors[0][0].degree()==3 and den_factors[0][1]==4
assert den==den.leading_coefficient()*den_factors[0][0]**4
sextic=q/den.leading_coefficient()
assert sextic.degree()==6 and sextic.gcd(e)==1
factors=list(sextic.factor())
assert len(factors)==2 and all(f.degree()==3 and mult==1 for f,mult in factors)
assert sextic.gcd(sextic.derivative())==1
assert all(f.is_irreducible() for f,mult in factors)

# x([2]P) = num2 / den2 on y^2=x^3+a*x+b.
a,b=F(dt.gm.a),F(dt.gm.b)
num2=X**4-2*a*X**2-8*b*X+a**2
den2=4*e
pullback=sum(sextic[i]*num2**i*den2**(6-i) for i in range(7))
assert pullback.degree()==24
assert pullback.gcd(pullback.derivative())==1
assert pullback.gcd(sextic)==1 and pullback.gcd(e)==1
algebra_ns=time.monotonic_ns()-start

# Exact comparison on every scalar trial of the frozen 64-accept query run.
guided=json.loads((HERE/'type_i_geom_character_collector.json').read_text())
assert guided['attempted_scalars']==259 and guided['accepted_queries']==64
rng=random.Random(guided['seed'])
first_reject=second_reject=accepted=0
start=time.monotonic_ns()
def direct_character(P):
    assert P!=dt.E(0) and L(P[0])!=dt.c0
    z=dt.inf+dt.d0/(L(P[0])-dt.c0)
    value=prod(z-v**p for v in dt.quartic_roots)
    assert value!=0 and value**(p**3)==value
    return value**((p**3-1)//2)
def reduced_character(P):
    value=sextic(F(P[0]))
    assert value!=0
    return value**((p-1)//2)
for _ in range(guided['attempted_scalars']):
    scalar=rng.randrange(1,int(r))
    P=scalar*dt.G
    first=direct_character(P)
    assert first==reduced_character(P)
    if first!=1:
        first_reject+=1
        continue
    second=direct_character(2*P)
    assert second==reduced_character(2*P)
    if second!=1:
        second_reject+=1
    else:
        accepted+=1
assert (first_reject,second_reject,accepted)==(128,67,64)
query_replay_ns=time.monotonic_ns()-start

# F(x) branches at 12 points of E; F(x([2]P)) at 48 disjoint points.
# The double covers z^2=F(x) and z^2=F(x)F(x([2]P)) have genera 7,31.
# Hasse--Weil gives |S1| <= 14 sqrt(p)+|t| and
# |S12| <= 62 sqrt(p)+|t|, where S1,S12 are character sums on E(Fp).
g1,g12=7,31
trace=int(p+1-r)
U=isqrt(int(p))+1
assert U*U>int(p) and (U-1)**2<int(p)
chi_lc=sextic.leading_coefficient()**((p-1)//2)
assert chi_lc in (F(1),F(-1))
identity_accepted=1 if chi_lc==1 else 0
error=(4*g1+2*g12)*U+3*abs(trace)
lower_all=(int(r)-error+3)//4
lower_nonzero=lower_all-identity_accepted
assert lower_nonzero>0
hit_bound=Fraction(816,lower_nonzero)
assert 816*(1<<180)<lower_nonzero
out={
    'status':'checked_type_i_character_density_inputs',
    'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'p':str(p),'r':str(r),'trace':str(trace),
    'norm_numerator_degree':int(num.degree()),
    'norm_denominator_degree':int(den.degree()),
    'norm_numerator_has_elliptic_cubic_once':True,
    'norm_denominator_cubic_fourth_power':True,
    'sextic_coefficients_decimal':list(map(str,sextic.list())),
    'sextic_factor_coefficients_decimal':[list(map(str,f.list())) for f,m in factors],
    'sextic_factors_irreducible':True,
    'sextic_squarefree_and_coprime_to_elliptic_cubic':True,
    'doubling_pullback_degree':int(pullback.degree()),
    'doubling_pullback_squarefree':True,
    'doubling_pullback_coprime_to_sextic_and_elliptic_cubic':True,
    'cover_genera_by_riemann_hurwitz':[g1,g12],
    'ceil_sqrt_p':str(U),
    'leading_coefficient_character':1 if chi_lc==1 else -1,
    'identity_accepted':identity_accepted,
    'frozen_attempts_replayed':259,
    'frozen_first_rejects':first_reject,
    'frozen_second_rejects':second_reject,
    'frozen_accepted':accepted,
    'accepted_nonzero_scalar_lower_bound':str(lower_nonzero),
    'fixed_base_three_sum_support_upper_bound':816,
    'conditioned_hit_probability_upper_bound':str(hit_bound),
    'conditioned_hit_probability_lt_two_pow_minus_180':True,
    'algebra_ns':algebra_ns,'query_replay_ns':query_replay_ns,
    'proof_boundary':'exact finite-field checks; genus and character-sum bounds use written Riemann--Hurwitz and Hasse--Weil, not full Lean formalization',
}
print(json.dumps(out,indent=2,sort_keys=True))
