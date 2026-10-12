"""Matched chained S3 and elementary-symmetric S4 for arbitrary a6."""
import decomp
import ir
import nagaodecomp


def encode(f, curve, d, target, c, variant):
    n = f.m
    p = ir.Prog()
    inp = lambda name: [p.addInput(name, i) for i in range(n)]
    xs = [inp('p%d' % j) for j in range(3)]
    r, b = inp('r'), inp('B')
    add = lambda *vs: [p.xorList([v[i] for v in vs]) for i in range(n)]
    mul = lambda x, y: decomp.onbMul(p, x, y, n, min(n, 12))
    sq = lambda x: decomp.applyPerm(x, decomp.sigmaPerm(n, f.n, 1))
    x, y, z = xs
    if variant == 'chained-s3':
        t = inp('t')
        def s3(u, v, w):
            uv = mul(u, v)
            return add(sq(add(uv, mul(u, w), mul(v, w))), mul(uv, w), b)
        outputs = s3(x, y, t) + s3(z, r, t)
    elif variant == 's4-symmetric':
        e1 = add(x, y, z, r)
        xy, xz, yz = mul(x, y), mul(x, z), mul(y, z)
        xyz = mul(xy, z)
        e3 = add(xyz, mul(r, add(xy, xz, yz)))
        e4 = mul(r, xyz)
        ep = add(e4, b)
        outputs = add(sq(sq(e3)), mul(ep, sq(e3)),
                      mul(mul(e4, ep), sq(e1)), mul(sq(b), sq(sq(e1))))
    else:
        raise ValueError(variant)
    pvars = [[c.newVar() for _ in range(n)] for _ in range(3)]
    lits = {('p%d' % j, i): v for j, vs in enumerate(pvars) for i, v in enumerate(vs)}
    if variant == 'chained-s3':
        for i in range(n):
            lits['t', i] = c.newVar()
    for name, value in [('r', target[0]), ('B', curve.a6)]:
        coords = f.toCoords(value)
        for i in range(n):
            lits[name, i] = c.true if coords >> i & 1 else c.false
    for v in p.emitCnf(outputs, lits, c):
        c.assertZero(v)
    nagaodecomp.addRestrictedDomain(c, pvars, f.toCoords(target[0]))
    for vs in pvars:
        for v in vs[d:]:
            c.addClause([-v])
    return pvars


def validateIdentity():
    """Exact F2[x1,x2,x3,x4,B] resultant identity, no field sampling."""
    def add(*ps):
        out = set()
        for p in ps:
            out.symmetric_difference_update(p)
        return out
    def mul(p, q):
        out = set()
        for a in p:
            for b in q:
                v = tuple(x + y for x, y in zip(a, b))
                if v in out:
                    out.remove(v)
                else:
                    out.add(v)
        return out
    sq = lambda p: mul(p, p)
    x, y, z, r, b = [{tuple(int(i == j) for i in range(5))} for j in range(5)]
    # Res(A t²+D t+C, E t²+F t+G) in characteristic two.
    a, d, c = sq(add(x, y)), mul(x, y), add(sq(mul(x, y)), b)
    e, ff, g = sq(add(z, r)), mul(z, r), add(sq(mul(z, r)), b)
    resultant = add(sq(add(mul(a, g), mul(c, e))),
                    mul(add(mul(a, ff), mul(d, e)), add(mul(d, g), mul(c, ff))))
    e1 = add(x, y, z, r)
    e3 = add(mul(mul(x, y), z), mul(r, add(mul(x, y), mul(x, z), mul(y, z))))
    e4 = mul(mul(x, y), mul(z, r))
    symmetric = add(sq(sq(e3)), mul(add(e4, b), sq(e3)),
                    mul(mul(e4, add(e4, b)), sq(e1)), mul(sq(b), sq(sq(e1))))
    if resultant != symmetric:
        raise ArithmeticError('general-a6 symmetric S4/resultant identity failed')
    return {'general_a6_s4_exact_resultant_identity': True, 'terms': len(symmetric)}
