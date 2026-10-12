"""Function-coefficient circuit with conditioned quadratic linear coefficient."""
import decomp
import field
import ir
import nagaoannihilator as scalar


def imageBasis(f,d,uBits):
    u=f.fromCoords(uBits)
    if not 0<uBits<1<<d: raise ValueError('u must be nonzero in V')
    pivots={}
    for i in range(d):
        w=f.fromCoords(1<<i)
        v=f.toCoords(f.add(f.sqr(w),f.mul(u,w)))
        for k,b in sorted(pivots.items(),reverse=True):
            if v>>k&1:v^=b
        if v:pivots[v.bit_length()-1]=v
    assert len(pivots)==d-1
    return pivots


def reduceBits(v,pivots):
    for k,b in sorted(pivots.items(),reverse=True):
        if v>>k&1:v^=b
    return v


def build(n):
    f=field.Onb(n);p=ir.Prog()
    inp=lambda name:[p.addInput(name,i) for i in range(n)]
    a,b,r,s,u=[inp(name) for name in ('a','b','r','s','u')]
    add=lambda *vs:[p.xorList([v[i] for v in vs]) for i in range(n)]
    mul=lambda x,y:decomp.onbMul(p,x,y,n,min(n,12))
    sq=lambda x:decomp.applyPerm(x,decomp.sigmaPerm(n,f.n,1))
    c=add(sq(r),mul(a,r),mul(b,add(r,s)))
    h2=add(b,sq(b),r)
    h1=add(sq(a),mul(a,b),mul(r,h2))
    h0=add(mul(b,c),mul(r,h1))
    x=add(h2,u)
    v=add(h1,mul(x,u))
    hx=add(h0,mul(x,v))
    hr=add(mul(add(mul(add(r,h2),r),h1),r),h0)
    qx=add(sq(x),mul(u,x),v)
    return p,sum([x,v,hx,h0,hr,qx],[])


def encode(p,outputs,f,d,target,uBits,c):
    n=f.m;variables={name:[c.newVar() for _ in range(n)] for name in ('a','b')}
    vals={'r':f.toCoords(target[0]),'s':f.toCoords(target[1]),'u':uBits}
    lits={(name,i):v for name,vs in variables.items() for i,v in enumerate(vs)}
    lits.update({(name,i):c.true if value>>i&1 else c.false for name,value in vals.items() for i in range(n)})
    out=p.emitCnf(outputs,lits,c)
    x,v,hx,h0,hr,qx=[out[k*n:(k+1)*n] for k in range(6)]
    for lit in x[d:]+hx:c.assertZero(lit)
    for row in (h0,hr,qx,variables['b']):c.addClause(row)
    basis=imageBasis(f,d,uBits)
    remainders=[reduceBits(1<<i,basis) for i in range(n)]
    for k in range(n):
        row=[v[i] for i in range(n) if remainders[i]>>k&1]
        if row:c.addXor(row,False)
    return variables
