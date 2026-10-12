"""Target-only four-point decomposition on E: y^2 + xy = x^3 + 1 over F_2^n,
encoded as an explicit finite-field circuit for CryptoMiniSat.

Unknowns: the factor-space coordinates c[i] (s bits per summand), the
x-coordinate t of P1 + P2 (up to signs), and either u (x of P3 + P4,
"specialized") or the slope l of the line through P1+P2, P3+P4 and -R
("line").  Every field product z = x*y of linear forms becomes n^2 AND gates
plus n XOR constraints; squaring, sqrt, constant multiplication and trace are
linear and stay inside XOR constraints.  With --cnf every XOR is expanded into
ordinary clauses instead (Tseitin chunks of 4).

The model is given only the target R, the field, and the factor-space basis.
It is never given valid payloads, orbit keys, or a decomposition.
"""
import itertools, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..",
                                "frobenius_quotient_followup_20260926"))
from verify_claims import F, Curve  # noqa: E402


class Lin:
    """Linear form over GF(2): XOR of variables plus a constant."""
    __slots__ = ("v", "c")

    def __init__(self, v=frozenset(), c=0):
        self.v, self.c = frozenset(v), c

    def __add__(self, o):
        return Lin(self.v ^ o.v, self.c ^ o.c)


def lsum(forms):
    v, c = set(), 0
    for f in forms:
        v ^= f.v
        c ^= f.c
    return Lin(v, c)


class Model:
    def __init__(self, n, native_xor=True):
        self.f = F(n)
        self.n = n
        self.native_xor = native_xor
        self.nv = 0
        self.clauses, self.xors = [], []
        self.stats = {"and_gates": 0, "products": 0}
        # linear maps as matrices of column images
        self.sq_cols = [self.f.sq(1 << i) for i in range(n)]
        self.sqrt_cols = [self.f.sqrt(1 << i) for i in range(n)]
        self.tr_row = [self.f.tr(1 << i) for i in range(n)]
        self.mul_tab = [self.f.mul(1 << i, 1 << j)
                        for i in range(n) for j in range(n)]

    # -- variables -------------------------------------------------------
    def var(self):
        self.nv += 1
        return self.nv

    def elem(self):
        return [Lin({self.var()}) for _ in range(self.n)]

    def const(self, a):
        return [Lin((), (a >> k) & 1) for k in range(self.n)]

    # -- linear maps -----------------------------------------------------
    def lin_map(self, x, cols):
        return [lsum(x[i] for i in range(self.n) if (cols[i] >> k) & 1)
                for k in range(self.n)]

    def add(self, x, y):
        return [a + b for a, b in zip(x, y)]

    def sq(self, x):
        return self.lin_map(x, self.sq_cols)

    def sqrt(self, x):
        return self.lin_map(x, self.sqrt_cols)

    def cmul(self, a, x):
        return self.lin_map(x, [self.f.mul(a, 1 << i) for i in range(self.n)])

    def trace(self, x):
        return lsum(x[i] for i in range(self.n) if self.tr_row[i])

    # -- constraints -----------------------------------------------------
    def xor_eq(self, form, rhs=0):
        """form == rhs."""
        vs, r = sorted(form.v), form.c ^ rhs
        if not vs:
            assert r == 0, "inconsistent constant constraint"
            return
        self.xors.append((vs, bool(r)))

    def zero(self, x):
        for form in x:
            self.xor_eq(form)

    def as_var(self, form):
        if len(form.v) == 1 and form.c == 0:
            return next(iter(form.v))
        v = self.var()
        self.xor_eq(form + Lin({v}))
        return v

    def mul(self, x, y):
        """Field product of two elements given as linear forms."""
        n = self.n
        self.stats["products"] += 1
        xv = [self.as_var(a) for a in x]
        yv = [self.as_var(b) for b in y]
        terms = [[] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                g = self.var()
                a, b = xv[i], yv[j]
                self.clauses += [[-g, a], [-g, b], [g, -a, -b]]
                self.stats["and_gates"] += 1
                m = self.mul_tab[i * n + j]
                for k in range(n):
                    if (m >> k) & 1:
                        terms[k].append(g)
        return [Lin(t) for t in terms]

    def differ(self, a, b):
        """Bit vectors a, b (lists of var ids) are not equal."""
        ds = []
        for p, q in zip(a, b):
            d = self.var()
            self.xor_eq(Lin({p, q, d}))
            ds.append(d)
        self.clauses.append(ds)

    # -- output ----------------------------------------------------------
    def load(self, solver):
        for c in self.clauses:
            solver.add_clause(c)
        if self.native_xor:
            for vs, r in self.xors:
                solver.add_xor_clause(vs, r)
        else:
            for c in self._xor_cnf():
                solver.add_clause(c)

    def _xor_cnf(self):
        out = []
        nv = self.nv
        for vs, r in self.xors:
            vs = list(vs)
            while len(vs) > 4:  # chop: head XOR = aux
                nv += 1
                head, vs = vs[:3], vs[3:] + [nv]
                out += _xor_clauses(head + [nv], False)
            out += _xor_clauses(vs, r)
        self.nv = nv
        return out

    def size(self):
        return {"vars": self.nv, "clauses": len(self.clauses),
                "xors": len(self.xors), **self.stats}


def _xor_clauses(vs, r):
    out = []
    for signs in itertools.product((0, 1), repeat=len(vs)):
        # forbid assignments whose parity != r: clause negates that assignment
        if sum(signs) % 2 != r:
            out.append([-v if s else v for v, s in zip(vs, signs)])
    return out


def f3(M, x1, x2, x3_is_var):
    """Constrain S3(x1, x2, t) = 0; returns nothing.  Semaev S3 for
    y^2+xy=x^3+1:  (x1x2 + x1x3 + x2x3)^2 + x1x2x3 + 1."""
    t = x3_is_var
    m = M.mul(x1, x2)
    w = M.mul(t, M.add(x1, x2))
    v = M.mul(m, t)
    M.zero(M.add(M.add(M.sq(M.add(m, w)), v), M.const(1)))


def build(n, basis, target, final="line", subgroup=False, native_xor=True,
          phases=False):
    """Return (model, c_vars, e_vars) for target R = (a, b).

    phases=False: each x_i is the factor-space representative v_i itself.
    phases=True:  x_i = v_i^(2^k_i) with the Frobenius phase k_i unknown,
    selected by a one-hot vector e_i (n AND gates per bit of x_i)."""
    M = Model(n, native_xor)
    a, b = target
    s = len(basis)
    cvars = [[M.var() for _ in range(s)] for _ in range(4)]
    vs = []
    for ci in cvars:
        # v = sum_j c_j * basis_j  (linear)
        vs.append([lsum(Lin({ci[j]}) for j in range(s) if (basis[j] >> k) & 1)
                   for k in range(n)])
    for v in vs:
        M.xor_eq(M.trace(v), 0)
        if subgroup:  # q^2 + sqrt(v) q + 1 = 0, Tr(q) = 0 (Frobenius-invariant)
            q = M.elem()
            g = M.mul(M.sqrt(v), q)
            M.zero(M.add(M.add(M.sq(q), g), M.const(1)))
            M.xor_eq(M.trace(q), 0)
    if not subgroup:  # certificate already excludes v = 0
        for ci in cvars:
            M.clauses.append(list(ci))
    evars = None
    if phases:
        evars, xs = [], []
        for v in vs:
            e = [M.var() for _ in range(n)]
            M.clauses.append(list(e))
            for i, j in itertools.combinations(e, 2):
                M.clauses.append([-i, -j])
            evars.append(e)
            frob = [[Lin({M.as_var(f)}) for f in v]]
            for _ in range(n - 1):  # Frob^k(v), linear in v
                frob.append(M.sq(frob[-1]))
            x = []
            for bit in range(n):
                terms = set()
                for k in range(n):
                    form = frob[k][bit]
                    g = M.var()  # g = e_k AND form
                    fv = M.as_var(form) if form.v else None
                    if fv is None:
                        if form.c:
                            M.xor_eq(Lin({g, e[k]}))
                        else:
                            M.clauses.append([-g])
                    else:
                        M.clauses += [[-g, e[k]], [-g, fv], [g, -e[k], -fv]]
                        M.stats["and_gates"] += 1
                    terms ^= {g}
                x.append(Lin(terms))
            xs.append(x)
        xvars = [[M.as_var(f) for f in x] for x in xs]
        for i, j in itertools.combinations(range(4), 2):
            M.differ(xvars[i], xvars[j])
    else:
        xs = vs
        for i, j in itertools.combinations(range(4), 2):
            M.differ(cvars[i], cvars[j])
    t = M.elem()
    if final == "line":
        ell = M.elem()
        # u = t + l^2 + l + a ;  t u = a l^2 + a^2 + b + a
        u = M.add(M.add(M.add(t, M.sq(ell)), ell), M.const(a))
        p = M.mul(t, u)
        rhs = M.add(M.cmul(a, M.sq(ell)), M.const(M.f.sq(a) ^ b ^ a))
        M.zero(M.add(p, rhs))
    elif final == "specialized":
        u = M.elem()
        # (tu + a(t+u))^2 + a tu + 1 = 0
        p = M.mul(t, u)
        inner = M.add(p, M.cmul(a, M.add(t, u)))
        M.zero(M.add(M.add(M.sq(inner), M.cmul(a, p)), M.const(1)))
    else:
        raise ValueError(final)
    f3(M, xs[0], xs[1], t)
    f3(M, xs[2], xs[3], u)
    return M, cvars, evars


def decode(sol, cvars, basis, evars=None, f=None):
    xs = []
    for i, ci in enumerate(cvars):
        x = 0
        for j, v in enumerate(ci):
            if sol[v]:
                x ^= basis[j]
        if evars is not None:
            k = next(k for k, e in enumerate(evars[i]) if sol[e])
            for _ in range(k):
                x = f.sq(x)
        xs.append(x)
    return xs


def verify(E, xs, target):
    """Independent check: distinct x in factor space are rational points in
    the odd subgroup, and some signed sum equals the target."""
    if len(set(xs)) != 4:
        return None
    pts = []
    for x in xs:
        P = E.lift(x)
        if P is None or x == 0 or E.mul(E.h, P) is not None:
            return None
        pts.append(P)
    for signs in itertools.product((0, 1), repeat=4):
        S = None
        for P, sg in zip(pts, signs):
            S = E.add(S, E.neg(P) if sg else P)
        if S == target:
            return [P if not sg else E.neg(P) for P, sg in zip(pts, signs)]
    return None
