import sys; sys.path.insert(0,'.')
import ghs3, random
from ghs3 import *
src=open('ghs3.py').read()
# instrument: copy build_instance body pieces quickly
rng=random.Random(1); p=67
K=KField(p,rng)
pr=K.rand(rng); r0=pr; r1=K.frob(r0); r2=K.frob(r1)
print('r', r0, r1, r2, 'frob3', K.frob(r2)==r0)
c0=K.sub(r2,r0)
trb=[]
for e in ((1,0,0),(0,1,0),(0,0,1)):
    x=e; s=KZ
    for _ in range(3): s=K.add(s,x); x=K.frob(x)
    trb.append(s)
print('trb',trb)
def Tr(x): return (x[0]*trb[0][0]+x[1]*trb[1][0]+x[2]*trb[2][0])%p
def Qf(s): return Tr(K.mul(c0,K.sqr(s)))
cnt=0; xs=[]
import itertools
pts=[]
for s in itertools.product(range(p),repeat=3):
    if s==(0,0,0): continue
    if Qf(s)==0:
        X0=s; X2=K.frob(K.frob(s))
        n=K.sub(K.mul(r0,K.sqr(X0)),K.mul(r1,K.sqr(X2))); d=K.sub(K.sqr(X0),K.sqr(X2))
        x=K.mul(n,K.inv(d)) if d!=KZ else None
        pts.append((s,x))
print(len(pts), 'conic pts (affine, incl scalars)', len(pts)//(p-1))
print('x values in Fp:', sum(1 for s,x in pts if x is not None and K.in_fp(x)))
print([x for s,x in pts[:5]])
from collections import Counter
cx=Counter(x for s,x in pts if x is not None)
print('x multiplicities', Counter(cx.values()))
a,b = 5, 9
print('pts over a', cx.get((a,0,0)), 'pts over b', cx.get((b,0,0)))
def Bf(s1,s2): return Tr(K.mul(c0,K.mul(s1,s2)))
sstar=pts[0][0]
e1,e2=K.rand(rng),K.rand(rng)
q0,q1,q2=Qf(e1),2*Bf(e1,e2)%p,Qf(e2)
b0,b1=Bf(sstar,e1),Bf(sstar,e2)
c=[K.sub(K.smul(e1,2*b0%p),K.smul(sstar,q0)),
   K.sub(K.add(K.smul(e1,2*b1%p),K.smul(e2,2*b0%p)),K.smul(sstar,q1)),
   K.sub(K.smul(e2,2*b1%p),K.smul(sstar,q2))]
def ev(c,z):
    acc=KZ
    for cc in reversed(c): acc=K.add(K.smul(acc,z),cc)
    return acc
print('Q(s(z)) for z=0..5', [Qf(ev(c,z)) for z in range(6)])
X0=ktr(c[:]); X2=[K.frob(K.frob(cc)) for cc in X0]
X0s,X2s=kpmul(K,X0,X0),kpmul(K,X2,X2)
Delta=kpsub(K,X0s,X2s); Num=kpsub(K,kpscal(K,X0s,r0),kpscal(K,X2s,r1))
kap=K.inv(Delta[-1]); xd=kpscal(K,Delta,kap); xn=kpscal(K,Num,kap)
print('xd',xd); print('xn',xn)
