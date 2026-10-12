"""Exact F2 spaces, including F4-linear spaces stable under fourth powers."""
import algebra


def independent(values):
    pivots = {}
    result = []
    for original in values:
        value = original
        for k, column in sorted(pivots.items(), reverse=True):
            if value >> k & 1:
                value ^= column
        if value:
            pivots[value.bit_length() - 1] = value
            result.append(original)
    return result


def kernel(columns):
    pivots = {}
    result = []
    for i, column in enumerate(columns):
        value, pre = column, 1 << i
        for k, (v, p) in sorted(pivots.items(), reverse=True):
            if value >> k & 1:
                value ^= v
                pre ^= p
        if value:
            pivots[value.bit_length() - 1] = value, pre
        else:
            result.append(pre)
    return result


def smallMul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return out


def remainder(p, q):
    p = list(p)
    while len(p) >= len(q):
        factor = p[-1]  # all divisors are monic
        shift = len(p) - len(q)
        for i, value in enumerate(q):
            p[shift + i] ^= smallMul(factor, value)
        while p and not p[-1]:
            p.pop()
    return p


def divide(p, q):
    p = list(p)
    quotient = [0] * (len(p) - len(q) + 1)
    while len(p) >= len(q):
        factor = p[-1]
        shift = len(p) - len(q)
        quotient[shift] = factor
        for i, value in enumerate(q):
            p[shift + i] ^= smallMul(factor, value)
        while p and not p[-1]:
            p.pop()
    if p:
        raise ArithmeticError('inexact polynomial division')
    return quotient


def factors(ell):
    p = [1] + [0] * (ell - 1) + [1]
    out = []
    for degree in range(1, 4):
        for bits in range(4 ** degree):
            q = [(bits >> (2 * i)) & 3 for i in range(degree)] + [1]
            if not q[0]:
                continue
            if len(p) >= len(q) and not remainder(p, q):
                out.append(q)
                p = divide(p, q)
    if p != [1]:
        raise ValueError('this exploratory field needs higher-degree factors')
    return out


class Space:
    def __init__(self, f, basisCoords):
        if len(independent(basisCoords)) != len(basisCoords):
            raise ValueError('dependent factor-base basis')
        self.f = f
        self.basisCoords = list(basisCoords)
        self.basis = [f.fromCoords(v) for v in basisCoords]
        self.d = len(self.basis)
        self.values = [0]
        for b in self.basis:
            self.values += [f.add(x, b) for x in self.values]
        self.indices = {v: i for i, v in enumerate(self.values)}

    def metadata(self, alpha):
        f = self.f
        return {'dimension': self.d, 'basis': self.basisCoords,
                'f4_linear': all(f.mul(alpha, b) in self.indices for b in self.basis),
                'frobenius4_stable': all(f.frob(b, 2) in self.indices for b in self.basis),
                'proper_subfield_containments': [k for k in range(1, f.m) if f.m % k == 0
                    and all(f.frob(b, k) == b for b in self.basis)]}


def makeBasis(f, d, kind, alpha):
    if not 0 < d <= f.m:
        raise ValueError('invalid dimension')
    if kind == 'prefix':
        return [1 << i for i in range(d)], {'construction': 'first ONB coordinates'}
    if kind != 'f4-stable' or d % 2:
        raise ValueError('F4-linear spaces have even binary dimension')
    ell = f.m // 2
    if ell % 2 == 0:
        raise ValueError('semisimple construction requires odd extension degree')
    ff = factors(ell)
    coeffs = [0, f.one(), alpha, f.add(alpha, f.one())]
    blocks = []
    for poly in ff:
        columns = []
        for i in range(f.m):
            b = f.fromCoords(1 << i)
            v = 0
            for j, a in enumerate(poly):
                if a:
                    v = f.add(v, f.mul(coeffs[a], f.frob(b, 2 * j)))
            columns.append(f.toCoords(v))
        basis = kernel(columns)
        if len(basis) != 2 * (len(poly) - 1):
            raise ArithmeticError('Frobenius component dimension mismatch')
        blocks.append(basis)
    # Deterministic choice: require full ambient-field support, then use the
    # lexicographically first irreducible-component subset of the right size.
    for mask in range(1, 1 << len(blocks)):
        basis = [x for j, block in enumerate(blocks) if mask >> j & 1 for x in block]
        if len(basis) != d:
            continue
        v = Space(f, basis)
        meta = v.metadata(alpha)
        if meta['proper_subfield_containments']:
            continue
        if not meta['f4_linear'] or not meta['frobenius4_stable']:
            raise ArithmeticError('structured basis invariant failed')
        return basis, {'construction': 'kernels of factors of T^(n/2)+1 over F4',
                       'factors': ff, 'selected_factor_indices': [j for j in range(len(blocks)) if mask >> j & 1],
                       **meta}
    raise ValueError('no full-field structured space of this dimension')
