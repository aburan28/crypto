"""solver_16: Weil descent of the two decomposition encodings into WDSat ANF.

Everything here is symbolic over F_2 in a POLYNOMIAL basis of GF(2^n).  The
factor base is V = { x : deg x < d } (the first d polynomial-basis
coordinates), which is exactly Trimoska's layout: products of two elements of
V have degree < 2d - 1 and products of three have degree < 3d - 2, so the
elementary symmetric functions of the summands live in 2d - 1 and 3d - 2
bits *unreduced*.  That is what keeps the S'4 system small.

Two encodings of "P1 + P2 + P3 = R with x(P_i) in V" are emitted, in the ANF
form WDSat reads, with the 3d abscissa bits as variables 1..3d in both:

  s4   Trimoska's S'4: auxiliaries e1 (d bits), e2 (2d-1), e3 (3d-2) defined
       from the x bits, then the 131 bits of S4(E1, E3, E4) = 0 where
       E_k are the elementary symmetric functions of {x1, x2, x3, x(R)}.
  rr   The Nagao / Riemann-Roch norm form: f = x^2 + a x + c + b y in L(4O)
       with f(-R) = 0, and N(X) = A(X)^2 + bX A(X) + b^2 (X^3 + 1) matched to
       (X + r)(X + x1)(X + x2)(X + x3).  The X^3 coefficient equation
       b + b^2 = r + e1 is F_2-linear in b, so b is solved symbolically
       (b = b0(e1) + beta, one free bit); a stays a block of 131 free
       variables.  See contract.json for why b is eliminated and a is not.

No type hints, camelCase identifiers, no itertools (project convention).
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CODE = ROOT / 'ecc2k130/codegen'
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
import field   # noqa: E402
import curves  # noqa: E402

CHALLENGE_R = 680564733841876926932320129493409985129


# ---------------------------------------------------------------- monomials

class Monomials:
    """Global monomial table.  Column 0 is the constant 1; columns 1..nvars
    are the variables; higher columns are non-unary monomials, allocated on
    first use.  A Boolean polynomial is an int bitmask over columns."""

    def __init__(self, nvars):
        self.nvars = nvars
        self.byKey = {}
        self.keys = [None] * (nvars + 1)
        for v in range(1, nvars + 1):
            self.byKey[(v,)] = v
            self.keys[v] = (v,)

    def column(self, key):
        c = self.byKey.get(key)
        if c is None:
            c = len(self.keys)
            self.byKey[key] = c
            self.keys.append(key)
        return c

    def mulMono(self, i, j):
        """Column of the product of two monomial columns (x^2 = x)."""
        if i == 0:
            return j
        if j == 0:
            return i
        key = tuple(sorted(set(self.keys[i]) | set(self.keys[j])))
        return self.column(key)

    def mulPoly(self, p, q):
        """Product of two Boolean polynomials (bitmask form)."""
        if p == 0 or q == 0:
            return 0
        r = 0
        pi = p
        while pi:
            i = (pi & -pi).bit_length() - 1
            pi &= pi - 1
            qj = q
            while qj:
                j = (qj & -qj).bit_length() - 1
                qj &= qj - 1
                r ^= 1 << self.mulMono(i, j)
        return r

    def degree(self, column):
        return 0 if column == 0 else len(self.keys[column])

    def polyTerms(self, p):
        """(constant, [columns]) of a polynomial."""
        const = p & 1
        cols = []
        p >>= 1
        c = 1
        while p:
            if p & 1:
                cols.append(c)
            p >>= 1
            c += 1
        return const, cols


CONST_ONE = 1  # bitmask of the constant polynomial 1


def varPoly(v):
    return 1 << v


# ------------------------------------------------------- symbolic GF(2^n)

class Descent:
    """Field elements whose n coordinate bits are Boolean polynomials."""

    def __init__(self, pb, mono):
        self.pb = pb
        self.m = pb.m
        self.mono = mono
        # red[k] = t^k mod F as a bitmask over the m coordinates, k < 3m
        self.red = []
        t = 1
        for k in range(3 * self.m):
            self.red.append(t)
            t <<= 1
            if (t >> self.m) & 1:
                t ^= pb.poly
        self.sqrPattern = [self.red[2 * i] for i in range(self.m)]

    def zero(self):
        return [0] * self.m

    def constant(self, c):
        return [CONST_ONE if (c >> k) & 1 else 0 for k in range(self.m)]

    def fromVars(self, varIds):
        """Element sum_k var_k t^k for the given coordinate variables."""
        e = self.zero()
        for k, v in enumerate(varIds):
            e[k] = varPoly(v)
        return e

    def fromCoefficients(self, polys):
        """Element sum_k p_k t^k with symbolic coefficients p_k (k may reach
        beyond m; reduction is applied)."""
        e = self.zero()
        for k, p in enumerate(polys):
            if p:
                pat = self.red[k]
                b = 0
                while pat:
                    if pat & 1:
                        e[b] ^= p
                    pat >>= 1
                    b += 1
        return e

    def add(self, a, b):
        return [x ^ y for x, y in zip(a, b)]

    def addConst(self, a, c):
        return [a[k] ^ (CONST_ONE if (c >> k) & 1 else 0) for k in range(self.m)]

    def mulConst(self, a, c):
        """a * c for a constant c: an F_2-linear map on the coordinates."""
        out = self.zero()
        for i in range(self.m):
            if not a[i]:
                continue
            pat = self.pb.mul(c, 1 << i)
            b = 0
            while pat:
                if pat & 1:
                    out[b] ^= a[i]
                pat >>= 1
                b += 1
        return out

    def square(self, a):
        out = self.zero()
        for i in range(self.m):
            if not a[i]:
                continue
            pat = self.sqrPattern[i]
            b = 0
            while pat:
                if pat & 1:
                    out[b] ^= a[i]
                pat >>= 1
                b += 1
        return out

    def mulUnreduced(self, a, b, degA, degB):
        """Coefficients c_0..c_{degA+degB-2} of the product polynomial, before
        reduction, where a and b are supported on coordinates < degA, < degB."""
        out = [0] * (degA + degB - 1)
        for i in range(degA):
            if not a[i]:
                continue
            for j in range(degB):
                if b[j]:
                    out[i + j] ^= self.mono.mulPoly(a[i], b[j])
        return out

    def mul(self, a, b):
        """Full symbolic product, reduced."""
        prods = {}
        for i in range(self.m):
            if not a[i]:
                continue
            for j in range(self.m):
                if b[j]:
                    prods[i + j] = prods.get(i + j, 0) ^ self.mono.mulPoly(a[i], b[j])
        out = self.zero()
        for k, p in prods.items():
            if p:
                pat = self.red[k]
                bit = 0
                while pat:
                    if pat & 1:
                        out[bit] ^= p
                    pat >>= 1
                    bit += 1
        return out

    def evaluate(self, e, assignment):
        """Numeric value of a symbolic element under a {var: bit} dict."""
        v = 0
        for k in range(self.m):
            if evalPoly(self.mono, e[k], assignment):
                v |= 1 << k
        return v


def evalPoly(mono, p, assignment):
    const, cols = mono.polyTerms(p)
    r = const
    for c in cols:
        t = 1
        for v in mono.keys[c]:
            if not assignment.get(v, 0):
                t = 0
                break
        r ^= t
    return r


# ------------------------------------------------------------ the systems

def s4Symmetric(D, E1, E3, E4):
    """S4 for y^2 + xy = x^3 + 1 in elementary-symmetric form (solver_06/15):
    e3^4 + (e4 + 1) e3^2 + e4 (e4 + 1) e1^2 + e1^4, on the four abscissas."""
    e3sq = D.square(E3)
    e1sq = D.square(E1)
    t1 = D.square(e3sq)
    t2 = D.mul(D.addConst(E4, 1), e3sq)
    t3 = D.mul(D.add(D.square(E4), E4), e1sq)   # e4^2 + e4 = e4 (e4 + 1)
    t4 = D.square(e1sq)
    return D.add(D.add(t1, t2), D.add(t3, t4))


class System:
    def __init__(self, name, nvars, equations, mono, layout):
        self.name = name
        self.nvars = nvars
        self.equations = equations   # list of Boolean polynomials (== 0)
        self.mono = mono
        self.layout = layout

    def stats(self):
        cols = set()
        maxTerms = 0
        maxDeg = 1
        totalTerms = 0
        for p in self.equations:
            const, cs = self.mono.polyTerms(p)
            maxTerms = max(maxTerms, len(cs) + 1)
            totalTerms += len(cs) + 1
            for c in cs:
                cols.add(c)
                maxDeg = max(maxDeg, self.mono.degree(c))
        nonUnary = sum(1 for c in cols if c > self.nvars)
        return {'vars': self.nvars, 'equations': len(self.equations),
                'distinct_monomials': len(cols), 'non_unary_monomials': nonUnary,
                'max_terms_per_equation': maxTerms, 'total_terms': totalTerms, 'max_degree': maxDeg,
                'wdsat_max_id': self.nvars + nonUnary}

    def cmsClauses(self):
        """CNF-XOR for CryptoMiniSat: one Tseitin AND per non-unary monomial
        (monomial column c <-> AND of its variables) and one native XOR clause
        per equation.  Returns (nvars_total, clauses, xors) with xors as
        (list_of_vars, rhs)."""
        colVar = {}
        nxt = self.nvars
        clauses = []
        xors = []
        for p in self.equations:
            const, cs = self.mono.polyTerms(p)
            lits = []
            for c in cs:
                key = self.mono.keys[c]
                if len(key) == 1:
                    lits.append(key[0])
                    continue
                v = colVar.get(c)
                if v is None:
                    nxt += 1
                    v = nxt
                    colVar[c] = v
                    for u in key:
                        clauses.append([-v, u])
                    clauses.append([v] + [-u for u in key])
                lits.append(v)
            # p = const + sum(lits) == 0  <=>  xor(lits) == const
            xors.append((lits, bool(const)))
        return nxt, clauses, xors

    def anfText(self):
        lines = ['p cnf %d %d' % (self.nvars, len(self.equations))]
        for p in self.equations:
            const, cols = self.mono.polyTerms(p)
            parts = ['x']
            # XOR clause is asserted TRUE; p == 0 <=> (terms) xor (const xor 1) == 1
            if const == 0:
                parts.append('T')
            for c in cols:
                key = self.mono.keys[c]
                if len(key) == 1:
                    parts.append(str(key[0]))
                else:
                    parts.append('.%d' % len(key))
                    parts.extend(str(v) for v in key)
            parts.append('0')
            line = ' '.join(parts)
            if len(line) >= 29000:
                raise ValueError('ANF line exceeds WDSat static clause buffer (30000)')
            lines.append(line)
        return '\n'.join(lines) + '\n'


def buildS4(pb, d, xR):
    """Trimoska's S'4 layout: x blocks, then e1, e2, e3 auxiliaries."""
    m = pb.m
    if 3 * d - 2 > m:
        raise ValueError('3d - 2 must not exceed n for the unreduced layout')
    nx = 3 * d
    nAux = d + (2 * d - 1) + (3 * d - 2)
    nvars = nx + nAux
    mono = Monomials(nvars)
    D = Descent(pb, mono)
    x = [D.fromVars(range(1 + blk * d, 1 + (blk + 1) * d)) for blk in range(3)]
    w = list(range(nx + 1, nx + d + 1))
    u = list(range(nx + d + 1, nx + d + 2 * d))
    v = list(range(nx + d + 2 * d, nx + nAux + 1))
    equations = []
    # e1 = x1 + x2 + x3, coordinatewise (d equations)
    for k in range(d):
        equations.append(varPoly(w[k]) ^ x[0][k] ^ x[1][k] ^ x[2][k])
    # e2 = x1x2 + x1x3 + x2x3, unreduced (2d - 1 equations)
    p12 = D.mulUnreduced(x[0], x[1], d, d)
    p13 = D.mulUnreduced(x[0], x[2], d, d)
    p23 = D.mulUnreduced(x[1], x[2], d, d)
    for k in range(2 * d - 1):
        equations.append(varPoly(u[k]) ^ p12[k] ^ p13[k] ^ p23[k])
    # e3 = x1x2x3, unreduced (3d - 2 equations)
    p123 = D.mulUnreduced(p12, x[2], 2 * d - 1, d)
    for k in range(3 * d - 2):
        equations.append(varPoly(v[k]) ^ p123[k])
    # S4 on {x1, x2, x3, xR} through the auxiliaries
    e1 = D.fromCoefficients([varPoly(t) for t in w])
    e2 = D.fromCoefficients([varPoly(t) for t in u])
    e3 = D.fromCoefficients([varPoly(t) for t in v])
    E1 = D.addConst(e1, xR)
    E3 = D.add(e3, D.mulConst(e2, xR))
    E4 = D.mulConst(e3, xR)
    s4 = s4Symmetric(D, E1, E3, E4)
    equations.extend(s4)
    layout = {'x_blocks': [[1 + blk * d, (blk + 1) * d] for blk in range(3)],
              'e1': [w[0], w[-1]], 'e2': [u[0], u[-1]], 'e3': [v[0], v[-1]]}
    return System('s4', nvars, equations, mono, layout)


def solveTraceEquation(pb, target):
    """Return b with b + b^2 = target, or None if Tr(target) = 1."""
    curve = curves.CurvePb(pb)
    if curve.trace(target):
        return None
    b = curve.halfTrace(target)
    assert pb.mul(b, b) ^ b == target
    return b


def buildRR(pb, d, r, s):
    """Nagao norm form with b eliminated linearly and a kept as variables.

    Variables: x blocks (3d), beta (1 bit, the kernel of b -> b + b^2),
    a (n bits).  Equations, each n Boolean bits:
      E2: a^2 + a b + r e1 + e2 = 0
      E3: r^2 b + r a b + (r + s) b^2 + r e2 + e3 = 0
      E4: r^4 + r^2 a^2 + ((r+s)^2 + 1) b^2 + r e3 = 0
    where b = b0(e1) + beta is the symbolic solution of E1: b + b^2 = r + e1.
    E1 requires Tr(r + e1) = 0; since e1 in V and the trace is linear, this
    is Tr(r) + Tr(e1) = 0, one linear equation on the x bits, emitted as E1.
    """
    m = pb.m
    nx = 3 * d
    beta = nx + 1
    aVars = list(range(nx + 2, nx + 2 + m))
    nvars = nx + 1 + m
    mono = Monomials(nvars)
    D = Descent(pb, mono)
    curve = curves.CurvePb(pb)
    x = [D.fromVars(range(1 + blk * d, 1 + (blk + 1) * d)) for blk in range(3)]
    e1 = D.add(D.add(x[0], x[1]), x[2])
    e2 = D.add(D.add(D.mul(x[0], x[1]), D.mul(x[0], x[2])), D.mul(x[1], x[2]))
    e3 = D.mul(D.mul(x[0], x[1]), x[2])
    # E1: b + b^2 = T with T = r + e1.  For odd n the half-trace H is
    # F_2-linear with H(y)^2 + H(y) = y + Tr(y), so b0 = H(T) solves E1
    # exactly when Tr(T) = 0, which is the one equation emitted as E1; the
    # kernel {0, 1} of b -> b + b^2 is the free bit beta.
    if m % 2 == 0:
        raise ValueError('half-trace elimination needs odd n')
    target = D.addConst(e1, r)
    traceEq = 0
    b0 = D.zero()
    for k in range(m):
        if not target[k]:
            continue
        if curve.trace(1 << k):
            traceEq ^= target[k]
        pat = curve.halfTrace(1 << k)
        bit = 0
        while pat:
            if pat & 1:
                b0[bit] ^= target[k]
            pat >>= 1
            bit += 1
    b = list(b0)
    b[0] ^= varPoly(beta)   # the kernel {0, 1} of b -> b + b^2
    a = D.fromVars(aVars)
    aSq = D.square(a)
    bSq = D.square(b)
    ab = D.mul(a, b)
    equations = [traceEq]
    E2 = D.add(D.add(aSq, ab), D.add(D.mulConst(e1, r), e2))
    rs = r ^ s
    E3 = D.add(D.add(D.mulConst(b, pb.mul(r, r)), D.mulConst(ab, r)),
               D.add(D.mulConst(bSq, rs), D.add(D.mulConst(e2, r), e3)))
    r2 = pb.mul(r, r)
    r4 = pb.mul(r2, r2)
    coef = pb.mul(rs, rs) ^ 1
    E4 = D.addConst(D.add(D.add(D.mulConst(aSq, r2), D.mulConst(bSq, coef)), D.mulConst(e3, r)), r4)
    equations.extend(E2)
    equations.extend(E3)
    equations.extend(E4)
    layout = {'x_blocks': [[1 + blk * d, (blk + 1) * d] for blk in range(3)],
              'beta': beta, 'a': [aVars[0], aVars[-1]], 'b': 'eliminated: b = b0(e1) + beta'}
    sysm = System('rr', nvars, equations, mono, layout)
    sysm.bSymbolic = b
    sysm.descent = D
    return sysm


# --------------------------------------------------------------- curve side

def admissibleAbscissae(pb, d):
    """Nonzero x of degree < d that lift to the curve."""
    curve = curves.CurvePb(pb)
    out = []
    for x in range(1, 1 << d):
        if curve.pointFromX(x) is not None:
            out.append(x)
    return out


def sampleAdmissible(pb, d, rng):
    curve = curves.CurvePb(pb)
    while True:
        x = rng.getrandbits(d)
        if x and curve.pointFromX(x) is not None:
            return x


def plantedTarget(pb, d, rng, maxAttempts=100000):
    """R = P1 + P2 + P3 with distinct admissible abscissae in V, by sampling.
    Bounded: a base with fewer than three admissible abscissae (tiny d) raises
    instead of looping, the solver_14 defect."""
    curve = curves.CurvePb(pb)
    for _ in range(maxAttempts):
        xs = [sampleAdmissible(pb, d, rng) for _ in range(3)]
        if len(set(xs)) < 3:
            continue
        pts = []
        for x in xs:
            p = curve.pointFromX(x)
            if rng.getrandbits(1):
                p = curve.neg(p)
            pts.append(p)
        R = curve.add(curve.add(pts[0], pts[1]), pts[2])
        if R is None or R[0] == 0 or R[0] in xs:
            continue
        return R, xs, pts
    raise ValueError('no planted target after %d attempts at d = %d' % (maxAttempts, d))


def uniformTarget(pb, rng):
    curve = curves.CurvePb(pb)
    while True:
        x = rng.getrandbits(pb.m)
        p = curve.pointFromX(x)
        if p is not None and x:
            if rng.getrandbits(1):
                p = curve.neg(p)
            return p


def decodeBlocks(bits, d):
    """x1, x2, x3 as ints from the first 3d characters of a WDSat assignment."""
    xs = []
    for blk in range(3):
        v = 0
        for k in range(d):
            if bits[blk * d + k] == '1':
                v |= 1 << k
        xs.append(v)
    return xs


def verifyProjected(pb, R, xs):
    """Independent check of a projected witness: lift each abscissa and look
    for signs with P1 + P2 + P3 = +-R.  Returns a dict with the verdict."""
    curve = curves.CurvePb(pb)
    # x = 0 is the 2-torsion point (0, 1); CurvePb.pointFromX declines it.
    pts = [(0, 1) if x == 0 else curve.pointFromX(x) for x in xs]
    degenerate = len(set(xs)) < 3 or 0 in xs or R[0] in xs
    if any(p is None for p in pts):
        return {'lifts': False, 'relation': False, 'degenerate': degenerate, 'signs': None}
    for signs in range(8):
        acc = None
        for i in range(3):
            p = pts[i] if not (signs >> i) & 1 else curve.neg(pts[i])
            acc = curve.add(acc, p)
        if acc is not None and acc[0] == R[0]:
            return {'lifts': True, 'relation': True, 'degenerate': degenerate, 'signs': signs,
                    'sum_is_R': acc == R}
    return {'lifts': True, 'relation': False, 'degenerate': degenerate, 'signs': None}


def pairEnumerationCount(nAdmissible):
    """Null object: unordered abscissa pairs, C(A, 2)."""
    return nAdmissible * (nAdmissible - 1) // 2


def wdsatConfig(stats, findAll):
    """config.h text sized to a System's stats (plus margin)."""
    maxId = stats['wdsat_max_id']
    anf = stats['vars'] + 1
    deg = stats['max_degree'] + 1
    xeq = stats['equations'] + 2
    xeqSize = maxId + 2
    eq = (stats['max_degree'] + 1) * stats['non_unary_monomials'] + 200
    # dimacs.c asserts the total number of XOR atoms stays below this
    buf = max(5000, 2 * stats['total_terms'] + 1000)
    lines = ['#define __XG_ENHANCED__']
    if findAll:
        lines.append('#define __FIND_ALL_SOLUTIONS__')
    lines += ['#ifdef __XG_ENHANCED__',
              '#define __MAX_ANF_ID__ %d' % anf,
              '#define __MAX_DEGREE__ %d' % deg,
              '#endif',
              '#define __MAX_ID__ %d' % maxId,
              '#define __MAX_BUFFER_SIZE__ %d' % buf,
              '#define __MAX_EQ__ %d' % eq,
              '#define __MAX_EQ_SIZE__ %d' % (deg + 1),
              '#define __MAX_XEQ__ %d' % xeq,
              '#define __MAX_XEQ_SIZE__ %d' % xeqSize]
    return '\n'.join(lines) + '\n', {'MAX_ANF_ID': anf, 'MAX_DEGREE': deg, 'MAX_ID': maxId,
                                     'MAX_BUFFER_SIZE': buf, 'MAX_EQ': eq, 'MAX_EQ_SIZE': deg + 1,
                                     'MAX_XEQ': xeq, 'MAX_XEQ_SIZE': xeqSize, 'FIND_ALL_SOLUTIONS': findAll}


def cmsSolve(sysm, d, timeLimit, enumerate=False, maxSolutions=64):
    """CryptoMiniSat (pycryptosat) on the CNF-XOR form.  Returns a dict with
    status ('SAT'/'UNSAT'/'TIMEOUT'), the projected witnesses found, and the
    wall time.  In enumerate mode each witness is blocked on its 3d abscissa
    bits and solving continues until UNSAT or the time limit."""
    import time
    import pycryptosat
    nv, clauses, xors = sysm.cmsClauses()
    solver = pycryptosat.Solver(threads=1, time_limit=float(timeLimit))
    solver.add_clauses(clauses)
    for lits, rhs in xors:
        if lits:
            solver.add_xor_clause(lits, rhs)
        elif rhs:
            solver.add_clause([])   # 0 == 1: trivially unsatisfiable
    t0 = time.time()
    witnesses = []
    status = None
    while True:
        remaining = timeLimit - (time.time() - t0)
        if remaining <= 0:
            status = 'TIMEOUT'
            break
        sat, model = solver.solve()
        if sat is None:
            status = 'TIMEOUT'
            break
        if not sat:
            status = 'UNSAT' if witnesses == [] or enumerate else status
            break
        bits = ''.join('1' if model[v] else '0' for v in range(1, 3 * d + 1))
        witnesses.append(bits)
        status = 'SAT'
        if not enumerate or len(witnesses) >= maxSolutions:
            break
        solver.add_clause([-(v) if model[v] else v for v in range(1, 3 * d + 1)])
    return {'status': status, 'witnesses': witnesses, 'wall': time.time() - t0,
            'cnf_vars': nv, 'clauses': len(clauses), 'xors': len(xors)}


def cmsConflictBracket(sysm, d, timeLimit, maxConflicts=1 << 40):
    """Smallest power-of-two conflict budget under which CryptoMiniSat decides
    the instance (SAT or UNSAT), by doubling then bisection.  Each probe is a
    fresh solver on identical input, so the result is deterministic for this
    pycryptosat build.  Returns (lo, hi, status): the instance is undecided at
    lo conflicts and decided at hi; None if the time limit is exhausted."""
    import time
    import pycryptosat
    nv, clauses, xors = sysm.cmsClauses()

    def probe(limit, budget):
        solver = pycryptosat.Solver(threads=1, confl_limit=int(limit), time_limit=float(budget))
        solver.add_clauses(clauses)
        for lits, rhs in xors:
            if lits:
                solver.add_xor_clause(lits, rhs)
            elif rhs:
                solver.add_clause([])
        sat, _ = solver.solve()
        return sat

    t0 = time.time()
    lo, hi = 0, 1
    status = None
    while hi <= maxConflicts:
        left = timeLimit - (time.time() - t0)
        if left <= 0:
            return None
        sat = probe(hi, left)
        if sat is not None:
            status = 'SAT' if sat else 'UNSAT'
            break
        lo, hi = hi, hi * 2
    if status is None:
        return None
    while hi - lo > 1:
        left = timeLimit - (time.time() - t0)
        if left <= 0:
            break
        mid = (lo + hi) // 2
        sat = probe(mid, left)
        if sat is None:
            lo = mid
        else:
            hi = mid
    return {'undecided_at': lo, 'decided_at': hi, 'status': status, 'wall': time.time() - t0}


def historyBytes(cfg):
    """Approximate static footprint of WDSat's largest arrays for a config."""
    idSize = cfg['MAX_ID'] + 1
    szGauss = idSize // 64 + 1
    hist = cfg['MAX_ANF_ID'] * idSize * szGauss * 8
    maskList = idSize * (idSize + 1) * 8
    mono = cfg['MAX_ANF_ID'] * idSize * (cfg['MAX_DEGREE'] - 1) * 8
    xeq = cfg['MAX_XEQ'] * cfg['MAX_XEQ_SIZE'] * 8
    return {'xorgauss_equivalency_history': hist, 'xorgauss_mask_list': maskList,
            'monomials_to_column': mono, 'xor_equations': xeq,
            'total': hist + maskList + mono + xeq}
