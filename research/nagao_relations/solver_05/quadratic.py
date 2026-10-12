"""Hybrid function search: solve b and a by Artin--Schreier equations."""
import curves
import nagaoannihilator as scalar


def asRoots(f,curve,c):
    if f.trace(c):return []
    w=curve.halfTrace(c)
    return [w,f.add(w,f.fromCoords((1<<f.m)-1))]


def candidates(f,curve,d,target):
    """Yield candidate functions (a,b) without scanning either coefficient field."""
    r,s=target
    # h2 equals the sum of the three abscissas, so it must lie in V.
    for hBits in range(1<<d):
        h2=f.fromCoords(hBits)
        for b in asRoots(f,curve,f.add(h2,r)):
            if not b:continue
            for zBits in range(1,1<<d):
                z=f.fromCoords(zBits)
                if z==r or z==h2:continue
                # H(z)=(r+z)*a²+b*z*a+K. Derive K by setting a=0.
                h,_=scalar.residualNorm(f,target,0,b)
                k=scalar.polyEval(f,h,z)
                invRz=f.inv(f.add(r,z))
                t=f.mul(f.mul(b,z),invRz)
                c=f.mul(k,invRz)
                invT=f.inv(t)
                for w in asRoots(f,curve,f.mul(c,f.sqr(invT))):
                    yield f.mul(t,w),b


def recover(f,curve,d,target,a,b):
    # The support condition remains exact H | L_V; evaluate it only after
    # constructing a candidate function, before extracting curve points.
    h,c=scalar.residualNorm(f,target,a,b)
    if not h[0] or not scalar.polyEval(f,h,target[0]):return None
    ls=scalar.subspacePolynomial(f,d)
    if any(scalar.annihilatorRemainder(f,ls,h)):return None
    xs=[f.fromCoords(x) for x in range(1,1<<d) if not scalar.polyEval(f,h,f.fromCoords(x))]
    if len(xs)!=3:raise ArithmeticError('support has wrong root count')
    invB=f.inv(b);total=None
    for x in xs:
        y=f.mul(f.add(f.add(f.sqr(x),f.mul(a,x)),c),invB)
        if not curve.onCurve((x,y)):raise ArithmeticError('off curve')
        total=curve.add(total,(x,y))
    if total!=target:raise ArithmeticError('sum mismatch')
    return tuple(sorted(f.toCoords(x) for x in xs))
