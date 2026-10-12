"""Same-x-base S4 controls, including the repository's transformed polynomial."""
import re
from pathlib import Path
import decomp
import field
import ir
import nagaodecomp


def terms():
    path=Path(__file__).resolve().parents[3]/'src/cryptanalysis/koblitz_symmetrised.rs'
    s=path.read_text()
    def extract(name):
        block=s.split('pub fn '+name+'(',1)[1].split('3 => vec![',1)[1].split('],\n        _ =>',1)[0]
        return [tuple(map(int,row.split(','))) for row in re.findall(r'vec!\[([0-9, ]+)\]',block)]
    return extract('plain_terms'),extract('symmetrised_terms')


def validateIdentity():
    # Exact polynomial arithmetic over F2, independent of field evaluation.
    def add(*ps):
        out=set()
        for p in ps:out.symmetric_difference_update(p)
        return out
    def mul(p,q):
        out=set()
        for a in p:
            for b in q:
                v=tuple(x+y for x,y in zip(a,b))
                if v in out:out.remove(v)
                else:out.add(v)
        return out
    one={(0,0,0,0)};xs=[{tuple(int(i==j) for i in range(4))} for j in range(4)]
    e1=add(*xs);e3=set()
    for i in range(4):
        t=one
        for j in range(4):
            if i!=j:t=mul(t,xs[j])
        e3=add(e3,t)
    e4=one
    for x in xs:e4=mul(e4,x)
    sq=lambda x:mul(x,x)
    direct=add(sq(sq(e3)),mul(add(e4,one),sq(e3)),mul(mul(e4,add(e4,one)),sq(e1)),sq(sq(e1)))
    plain,sym=terms()
    if direct!=set(plain):raise ArithmeticError('elementary S4 differs from repository polynomial')
    # Clear denominators of transformed polynomial in u_i=1/(x_i+1).
    # Since each u degree is at most four, multiply by product (x_i+1)^4.
    result=set()
    for ex in sym:
        product=one
        for i,power in enumerate(ex[:4]):
            # w_i=u_i²+u_i=x_i/(x_i+1)².
            factor=one
            for _ in range(power):factor=mul(factor,xs[i])
            for _ in range(4-2*power):factor=mul(factor,add(xs[i],one))
            product=mul(product,factor)
        if ex[4]:
            # s=sum u_i: each summand removes one denominator factor.
            product=set()
            for selected in range(4):
                term=one
                for i,power in enumerate(ex[:4]):
                    factor=one
                    for _ in range(power):factor=mul(factor,xs[i])
                    for _ in range(4-2*power-int(i==selected)):factor=mul(factor,add(xs[i],one))
                    term=mul(term,factor)
                product=add(product,term)
        result=add(result,product)
    if result!=set(plain):raise ArithmeticError('transformed S4 denominator identity mismatch')
    return {'elementary_polynomial_exact':True,'transformed_denominator_identity_exact':True,'plain_terms':len(plain),'transformed_terms':len(sym)}


def encode(n,d,target,c,variant):
    f=field.Onb(n);p=ir.Prog()
    inp=lambda name:[p.addInput(name,i) for i in range(n)]
    xs=[inp('p%d'%j) for j in range(3)];r=inp('r');one=inp('one')
    add=lambda *vs:[p.xorList([v[i] for v in vs]) for i in range(n)]
    mul=lambda x,y:decomp.onbMul(p,x,y,n,min(n,12))
    sq=lambda x:decomp.applyPerm(x,decomp.sigmaPerm(n,f.n,1))
    aux=[];roots=[];xr=f.toCoords(target[0]);oneBits=(1<<n)-1
    # x=1 lies outside first-d-coordinate subspaces when d<n. Target x=1
    # uses the exact elementary chart instead of dropping that target.
    if variant=='s4-transformed' and xr!=oneBits:
        aux=[inp('u%d'%j) for j in range(3)]
        ur=inp('ur');us=aux+[ur];ws=[add(sq(u),u) for u in us];ss=add(*us)
        for x,u in zip(xs,aux):roots+=add(mul(add(x,one),u),one)
        variables=ws+[ss];monomials=[]
        for ex in terms()[1]:
            factors=[sq(v) if e==2 else v for v,e in zip(variables,ex) if e]
            t=factors[0]
            for v in factors[1:]:t=mul(t,v)
            monomials.append(t)
        roots+=add(*monomials)
    else:
        x,y,z=xs
        e1=add(x,y,z,r)
        xy=mul(x,y);xz=mul(x,z);yz=mul(y,z);xyz=mul(xy,z)
        e3=add(xyz,mul(r,add(xy,xz,yz)));e4=mul(r,xyz)
        ep=add(e4,one)
        roots=add(sq(sq(e3)),mul(ep,sq(e3)),mul(mul(e4,ep),sq(e1)),sq(sq(e1)))
    lits={};pvars=[[c.newVar() for _ in range(n)] for _ in range(3)]
    for j,vs in enumerate(pvars):
        for i,v in enumerate(vs):lits['p%d'%j,i]=v
    for j in range(len(aux)):
        for i in range(n):lits['u%d'%j,i]=c.newVar()
    fixed={'r':xr,'one':oneBits}
    if aux:fixed['ur']=f.toCoords(f.inv(f.add(target[0],f.fromCoords(oneBits))))
    for name,value in fixed.items():
        for i in range(n):lits[name,i]=c.true if value>>i&1 else c.false
    for v in p.emitCnf(roots,lits,c):c.assertZero(v)
    nagaodecomp.addRestrictedDomain(c,pvars,xr)
    return p,roots,{'pvars':pvars}
