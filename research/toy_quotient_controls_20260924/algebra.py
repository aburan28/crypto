"""Finite symmetric quotient and coordinate controls, restricted to the toy field."""
import sys
from itertools import combinations_with_replacement, product
from pathlib import Path

LEGACY = Path(__file__).resolve().parent.parent/'toy_f5_neighbors_20260924'
sys.path.insert(0, str(LEGACY))
from curves import MUL, SQUARE, O, Curve, decomposition_value  # noqa: E402
from matrix import Ring, step, complete  # noqa: E402


def power(x, n):
    out = 1
    for _ in range(n):
        out = MUL[out][x]
    return out


def elementary(xs):
    coefficients = [1] + [0]*len(xs)
    for i, x in enumerate(xs):
        for j in range(i+1, 0, -1):
            coefficients[j] ^= MUL[x][coefficients[j-1]]
    return tuple(coefficients[1:])


def symmetric_value(e, target, b):
    """S3/S4 written in elementary symmetric coordinates, before root lifting."""
    if len(e) == 2:
        e1, e2 = e
        if target == O:
            return e1
        r = target[0]
        return SQUARE[e2 ^ MUL[r][e1]] ^ MUL[e2][r] ^ b
    e1, e2, e3 = e
    if target == O:
        return SQUARE[e2] ^ e3 ^ b
    r = target[0]
    # Exact symmetric rewrite of the S3/S3 resultant in characteristic two.
    monomials = [(b,b,r,r,r,r), (b,b,e1,e1,e1,e1), (b,r,r,r,e3),
                 (b,r,r,e2,e2), (b,r,e1,e1,e3), (b,e3,e3),
                 (r,r,r,r,e2,e2,e2,e2), (r,r,r,r,e3,e3),
                 (r,r,r,e2,e2,e3), (r,r,e1,e1,e3,e3),
                 (r,e3,e3,e3), (e3,e3,e3,e3)]
    out = 0
    for terms in monomials:
        term = 1
        for factor in terms:
            term = MUL[term][factor]
        out ^= term
    return out


def scaled_value(xs, target, b, u):
    """Evaluate on x'=u^2*x,y'=u^3*y,a1'=u,a6'=u^6*b."""
    u2 = SQUARE[u]
    xs = [MUL[u2][x] for x in xs]
    scaled_b = MUL[power(u,6)][b]

    def f3(x,y,z):
        pairs = MUL[x][y] ^ MUL[x][z] ^ MUL[y][z]
        return SQUARE[pairs] ^ MUL[u2][MUL[MUL[x][y]][z]] ^ MUL[u2][scaled_b]

    if target == O:
        return xs[0] ^ xs[1] if len(xs)==2 else f3(*xs)
    r = MUL[u2][target[0]]
    if len(xs)==2:
        return f3(*xs,r)
    x,y,z = xs
    a,bb,c = SQUARE[x^y], MUL[u2][MUL[x][y]], SQUARE[MUL[x][y]] ^ MUL[u2][scaled_b]
    d,ee,f = SQUARE[z^r], MUL[u2][MUL[z][r]], SQUARE[MUL[z][r]] ^ MUL[u2][scaled_b]
    return SQUARE[MUL[a][f]^MUL[c][d]] ^ MUL[MUL[a][ee]^MUL[bb][d]][MUL[bb][f]^MUL[c][ee]]


def setup(support, m, variant, contract):
    if m not in (2,3) or len(support)!=4 or any(len(p)!=2 for p in support):
        raise ValueError('only fixed four-x, two/three-summand toy supports are supported')
    xs = [p[0][0] for p in support]
    if len(set(xs))!=4:
        raise ValueError('support x-coordinates must be distinct')
    if variant == 'invariant':
        # Only coefficients are retained: reconstruction cannot read preimage tuples.
        coefficients = sorted({elementary([xs[i] for i in slots])
                               for slots in combinations_with_replacement(range(4),m)})
        if len(coefficients)!=(10 if m==2 else 20):
            raise AssertionError('symmetric map must separate multisets, including repetitions')
        return dict(n=(len(coefficients)-1).bit_length(), domain=coefficients)
    labels = contract['label_permutation'] if variant=='relabel' else list(range(4))
    return dict(n=2*m, labels=labels)


def encode(support, target, b, m, variant, prepared, contract):
    ring = Ring(prepared['n'])
    values, invalid = [], []
    for a in range(1<<ring.n):
        if variant=='invariant':
            bad = a>=len(prepared['domain'])
            values.append(0 if bad else symmetric_value(prepared['domain'][a],target,b))
        else:
            labels = tuple((a>>(2*j))&3 for j in range(m))
            slots = tuple(prepared['labels'][j] for j in labels)
            xs = [support[k][0][0] for k in slots]
            value = scaled_value(xs,target,b,contract['coordinate_scale']) if variant=='scale' else decomposition_value(xs,target,b)
            values.append(value)
            bad = variant=='canonical' and labels!=tuple(sorted(labels))
        invalid.append(int(bad))
    generators = [ring.anf([(value>>j)&1 for value in values]) for j in range(8)]
    if variant in ('canonical','invariant'):
        generators.append(ring.anf(invalid))
    if variant=='reverse_equations':
        generators.reverse()
    return ring, generators


def solve(ring, generators):
    """Single F5-criterion path; independent F4/point oracles are outside its cost."""
    generators = [g for g in generators if g]
    traces = []
    for degree in range(2*ring.n+1):
        rows, trace = step(ring,generators,degree,True)
        done, check_xors = complete(ring,rows,generators)
        trace.update(completion_check_xors=check_xors, ideal_complete=done)
        traces.append(trace)
        if done:
            roots = [a for a in range(1<<ring.n) if all(ring.evaluate(f,a)==0 for f in rows)]
            return dict(roots=roots,completion_degree=degree,traces=traces)
    raise AssertionError('bounded Boolean ideal failed to complete')


def divide_by_root(polynomial, root):
    quotient = [polynomial[0]]
    for coefficient in polynomial[1:-1]:
        quotient.append(coefficient ^ MUL[quotient[-1]][root])
    remainder = polynomial[-1] ^ MUL[quotient[-1]][root]
    return quotient, remainder


def recover(e, support):
    polynomial = [1]+list(e)
    slots = []
    for i, pair in enumerate(support):
        while len(polynomial)>1:
            quotient, remainder = divide_by_root(polynomial,pair[0][0])
            if remainder:
                break
            polynomial = quotient
            slots.append(i)
    if len(polynomial)!=1 or len(slots)!=len(e):
        raise ValueError('coefficient point does not split fully on the support')
    return tuple(slots)


def reconstruct(roots, support, m, variant, prepared):
    if variant=='invariant':
        if any(a>=len(prepared['domain']) for a in roots):
            raise AssertionError('padding guard failed')
        return [recover(prepared['domain'][a], support) for a in roots]
    return [tuple(prepared['labels'][(a>>(2*j))&3] for j in range(m)) for a in roots]


def verify(slots_list, support, target):
    answers = set()
    lift_checks = 0
    curve = Curve()  # all five normalized models share the same addition formulas
    # Identical output convention for all variants: unordered SIGNED point tuples.
    for slots in slots_list:
        for signs in product((0,1),repeat=len(slots)):
            points = tuple(support[i][sign] for i,sign in zip(slots,signs))
            lift_checks += 1
            if curve.total(points)==target:
                answers.add(tuple(sorted(points)))
    return sorted(answers), lift_checks
