"""Root-free Riemann--Roch support circuit; uses shared ONB and CNF machinery.

H is monic cubic. L_V mod H is computed by modular Frobenius squaring.
Only a,b are unknown field elements; target and subspace are fixed inputs.
"""
import curves
import decomp
import field
import ir
import nagaoannihilator as scalar


def build(n, d):
    f = field.Onb(n)
    if n < 3 or n % 2 == 0 or not 0 <= d <= n:
        raise ValueError('unsupported field or subspace dimension')
    p = ir.Prog()
    inp = lambda name: [p.addInput(name,i) for i in range(n)]
    a,b,r,s,one = [inp(name) for name in ('a','b','r','s','one')]
    ls = [inp('l%d'%j) for j in range(d+1)]
    zero = [None]*n
    add = lambda *vs: [p.xorList([v[i] for v in vs]) for i in range(n)]
    mul = lambda x,y: decomp.onbMul(p,x,y,n,min(12,n))
    perm = decomp.sigmaPerm(n,f.n,1)
    sq = lambda x: decomp.applyPerm(x,perm)
    c = add(sq(r),mul(a,r),mul(b,add(r,s)))
    h2 = add(b,sq(b),r)
    h1 = add(sq(a),mul(a,b),mul(r,h2))
    h0 = add(mul(b,c),mul(r,h1))
    hs = [h0,h1,h2]

    def squareMod(v):
        tmp = [sq(v[0]),zero,sq(v[1]),zero,sq(v[2])]
        for k in (4,3):
            top = tmp[k]
            for j in range(3):
                tmp[k-3+j] = add(tmp[k-3+j],mul(top,hs[j]))
        return tmp[:3]

    power = [zero,one,zero]
    rem = [zero,zero,zero]
    for j in range(d+1):
        rem = [add(rem[k],mul(ls[j],power[k])) for k in range(3)]
        if j < d:
            power = squareMod(power)
    hr = add(mul(add(mul(add(r,h2),r),h1),r),h0)
    return p, sum(rem,[]), sum([h0,hr],[])


def inputs(f,d,target,a=0,b=0):
    coeff = scalar.subspacePolynomial(f,d)
    vals = {'a':a,'b':b,'r':f.toCoords(target[0]),'s':f.toCoords(target[1]),
            'one':(1<<f.m)-1}
    vals.update({'l%d'%j:f.toCoords(v) for j,v in enumerate(coeff)})
    return vals


def encode(p,roots,nonzero,f,d,target,cnf):
    if target is None or not curves.Curve(f).onCurve(target):
        raise ValueError('valid finite target required')
    vals=inputs(f,d,target)
    lits={}
    variables={name:[cnf.newVar() for _ in range(f.m)] for name in ('a','b')}
    for name,value in vals.items():
        for i in range(f.m):
            lits[name,i] = variables[name][i] if name in variables else (
                cnf.true if value>>i&1 else cnf.false)
    outputs=p.emitCnf(roots+nonzero,lits,cnf)
    for lit in outputs[:3*f.m]:
        cnf.assertZero(lit)
    cnf.addClause(outputs[3*f.m:4*f.m])
    cnf.addClause(outputs[4*f.m:5*f.m])
    cnf.addClause(variables['b'])
    return variables


def recover(f,curve,d,target,aBits,bBits):
    if bBits==0:
        raise ArithmeticError('zero b')
    a,b=f.fromCoords(aBits),f.fromCoords(bBits)
    h,c=scalar.residualNorm(f,target,a,b)
    if not h[0] or not scalar.polyEval(f,h,target[0]):
        raise ArithmeticError('excluded zero or target root')
    xs=[f.fromCoords(x) for x in range(1,1<<d) if scalar.polyEval(f,h,f.fromCoords(x))==0]
    if len(xs)!=3:
        raise ArithmeticError('support model does not have three distinct roots in V')
    inv=f.inv(b)
    pts=[(x,f.mul(f.add(f.add(f.sqr(x),f.mul(a,x)),c),inv)) for x in xs]
    acc=None
    for pt in pts:
        if not curve.onCurve(pt):
            raise ArithmeticError('recovered point off curve')
        acc=curve.add(acc,pt)
    if acc!=target:
        raise ArithmeticError('signed group sum mismatch')
    return sorted(f.toCoords(x) for x in xs)
