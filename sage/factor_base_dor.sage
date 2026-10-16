# Factor base -> degree of regularity, measured.
#
# Index calculus decomposes a target R over a factor base {P : x(P) in V}
# by solving a Weil-descended Semaev system and running a Groebner basis.
# The cost is set by the solving degree D* (the degree of regularity): the
# highest degree the computation reaches. For a non-decomposable target the
# system is unsatisfiable and the computation derives 1; the degree at which
# 1 first appears in the Macaulay row space is D*, equal to the
# Polynomial-Calculus refutation degree and, up to a small additive term,
# the Huang-Kosters-Yeo last fall degree. See RESEARCH_FFD_PROOF_COMPLEXITY.md.
#
# This engine measures D* as a function of the factor base alone, so you can
# search for a base that lowers it. It reproduces the repo's live result:
# a subfield (multiplicatively closed) base has a lower D* than a random
# subspace of the same dimension, an algebraic effect that no graph invariant
# of the descent sees (RESEARCH_FFD_WORKFLOW.md, P2/P2'' killed, P3-alg live).
#
# Each summand may use its own subspace, which measures the symmetry
# breaking of Galbraith-Granger-Merz-Petit (SAC 2020, section 3.1): the
# nb_break structure takes P1 from a normal-basis window V1 and P2 from
# V2 = pi^dim(V1). The report gives the decomposable fraction next to D*,
# so the m! probability gain and its degree cost are compared directly.
#
# Scope. Over a binary field F_{2^n} the factor base is an F_2-subspace V of
# x-coordinates, and the descent is over F_2 -- the ECC2K-130 setting at toy
# size. Over a prime field there is no subfield to descend into, so DoR is a
# different object (Diem's degree-k extension, or symmetric-function solving
# degree over F_p); this engine does not cover it. See the note printed with
# --curve p256 and the ideas in the README.
#
# The curve is the Koblitz shape y^2 + x*y = x^3 + b with b = 1, so S_3 is
#   S_3(X1,X2,X3) = (X1+X2)^2 X3^2 + X1 X2 X3 + (X1 X2)^2 + b
# matching src/cryptanalysis/binary_semaev.rs::binary_semaev_s3.
#
#   sage factor_base_dor.sage --field 7 --dim 3 --targets 24
#   sage factor_base_dor.sage --field 9 --dim 3 --structures subfield,random,standard
#   sage factor_base_dor.sage --field 7 --dim 3 --dmax 6 --out dor-n7.json
#   sage factor_base_dor.sage --field 11 --dim 5 --targets 160 --structures nb_window,nb_break
#
# This file is preparsed: bitwise xor is ^^, and integer power stays **.

import argparse
import itertools
import json
import random
import sys


def popcount(mask):
    return int(mask).bit_count()


def rng_for(seed):
    """Stdlib RNG; Sage preparses integer defaults into Integer, which
    random.seed rejects, so coerce once here."""
    return random.Random(int(seed))


def struct_rng_for(seed, tag):
    """RNG for factor-base structure draws (random generators, geometric
    theta/g, normal elements). Must NOT reuse the target stream: with
    rng_for(seed) on both sides the first structure draws coincide with the
    first targets, so a target can land inside the base it is tested
    against and inflate the decomposable fraction (measured 2026-10-14:
    target[0] == theta at seed 1). A distinct string tag breaks the
    overlap while staying deterministic per cell."""
    return random.Random("%s|%s" % (int(seed), tag))


# ── GF(2^n) as F_2-coordinate vectors over the polynomial basis ────────

class BinaryField:
    """F_{2^n} = F_2[z]/(f) with the basis 1, z, ..., z^{n-1}.

    Field elements are length-n GF(2) coordinate vectors packed into an int.
    The F_2-bilinear multiplication tensor is precomputed so the same code
    multiplies constants and symbolic (polynomial-coefficient) elements.
    """

    def __init__(self, n):
        self.n = int(n)
        self.K = GF(2**self.n, name="z")
        self.z = self.K.gen()
        powers = [self.z**i for i in range(self.n)]
        self.basis_ints = [self._to_int(p) for p in powers]
        # tensor[i][j] = coordinate int of z^i * z^j
        self.tensor = [[self._to_int(powers[i] * powers[j])
                        for j in range(self.n)] for i in range(self.n)]
        self.one_int = self._to_int(self.K(1))

    def _to_int(self, elt):
        value = 0
        for i, coeff in enumerate(self.K(elt).polynomial().list()):
            if coeff:
                value |= 1 << i
        return value

    def from_int(self, value):
        value = int(value)
        acc = self.K(0)
        term = self.K(1)
        i = 0
        while value:
            if value & 1:
                acc += self.z**i
            value >>= 1
            i += 1
        return acc

    def to_int(self, elt):
        return self._to_int(elt)


# ── symbolic field elements: each coordinate is a GF(2) polynomial ─────
# A polynomial is a Python set of monomial bitmasks over the descent
# variables; the empty mask 0 is the constant 1. Addition is symmetric
# difference; multiplication unions masks, which reduces x_i^2 -> x_i, so
# everything stays multilinear (the field equations are built in).

def poly_add(a, b):
    return a ^^ b


def poly_mul(a, b):
    out = set()
    for ma in a:
        for mb in b:
            m = ma | mb
            if m in out:
                out.discard(m)
            else:
                out.add(m)
    return out


def fe_zero(n):
    return [set() for _ in range(n)]


def fe_add(u, v):
    return [poly_add(u[k], v[k]) for k in range(len(u))]


def fe_mul(field, u, v):
    n = field.n
    w = [set() for _ in range(n)]
    tensor = field.tensor
    for i in range(n):
        if not u[i]:
            continue
        for j in range(n):
            if not v[j]:
                continue
            prod = poly_mul(u[i], v[j])
            if not prod:
                continue
            bits = tensor[i][j]
            k = 0
            while bits:
                if bits & 1:
                    w[k] ^^= prod
                bits >>= 1
                k += 1
    return w


def fe_square(field, u):
    return fe_mul(field, u, u)


def fe_const(field, value_int):
    """Field constant from its coordinate int, as a symbolic element."""
    out = fe_zero(field.n)
    k = 0
    v = int(value_int)
    while v:
        if v & 1:
            out[k] = {0}
        v >>= 1
        k += 1
    return out


def fe_linear(field, basis_ints, var_offset):
    """X = sum_i a_i v_i with a_i a fresh variable (bit var_offset + i)."""
    out = fe_zero(field.n)
    for i, vi in enumerate(basis_ints):
        var_mask = 1 << (var_offset + i)
        k = 0
        bits = int(vi)
        while bits:
            if bits & 1:
                out[k] ^^= {var_mask}
            bits >>= 1
            k += 1
    return out


# ── the descended S_3 system for a factor base V ───────────────────────

def descent_equations(field, basis1, basis2, target_int, b_int):
    """The n GF(2) polynomials whose common zero is a decomposition of the
    target x-coordinate as x(P1), x(P2) with x(P1) in span(basis1) and
    x(P2) in span(basis2). Equal bases give the symmetric decomposition;
    distinct shifted bases give GGMP-style symmetry breaking."""
    X1 = fe_linear(field, basis1, 0)
    X2 = fe_linear(field, basis2, len(basis1))
    x3 = fe_const(field, target_int)
    bvec = fe_const(field, b_int)
    sum12 = fe_add(X1, X2)
    prod12 = fe_mul(field, X1, X2)
    term1 = fe_mul(field, fe_square(field, sum12), fe_square(field, x3))
    term2 = fe_mul(field, prod12, x3)
    term3 = fe_add(fe_square(field, prod12), bvec)
    s3 = fe_add(fe_add(term1, term2), term3)
    return [eq for eq in s3 if eq]


def descent_equations_m3(field, basis1, basis2, basis3, target_int, b_int):
    """Chained Semaev for a 3-point decomposition: x(P1)+x(P2) meet at an
    intermediate abscissa U, which with x(P3) sums to the target. Variables
    are the V-coordinates of P1, P2, P3 and the n full-field coordinates of
    U. Common zero <=> some sign choice of P1, P2, P3 sums to the target.
    This is the repo's Weil-chained form; it is how m >= 3 is handled."""
    d1, d2, d3 = len(basis1), len(basis2), len(basis3)
    X1 = fe_linear(field, basis1, 0)
    X2 = fe_linear(field, basis2, d1)
    X3 = fe_linear(field, basis3, d1 + d2)
    U = fe_linear(field, [1 << i for i in range(field.n)], d1 + d2 + d3)
    xr = fe_const(field, target_int)
    bvec = fe_const(field, b_int)

    def s3(a, c, e):
        s = fe_add(a, c)
        p = fe_mul(field, a, c)
        t1 = fe_mul(field, fe_square(field, s), fe_square(field, e))
        t2 = fe_mul(field, p, e)
        t3 = fe_add(fe_square(field, p), bvec)
        return fe_add(fe_add(t1, t2), t3)

    eqs = [eq for eq in s3(X1, X2, U) if eq]
    eqs += [eq for eq in s3(U, X3, xr) if eq]
    return eqs, d1 + d2 + d3 + field.n


def poly_degree(poly):
    return max((popcount(m) for m in poly), default=0)


# ── Macaulay refutation degree over GF(2) ──────────────────────────────

def refutation_degree(equations, nvars, dmax, shapes=None):
    """Smallest D in 2..dmax at which the constant 1 lies in the GF(2) row
    space of the degree-D Macaulay matrix (the solving degree D*), or None.
    Returns (D*, rank_profile) with the row rank at each degree reached.
    If `shapes` is a dict it receives {D: (rows, cols)} for each degree."""
    monomials = [0]
    for k in range(1, nvars + 1):
        for comb in itertools.combinations(range(nvars), k):
            m = 0
            for i in comb:
                m |= 1 << i
            monomials.append(m)
    profile = {}
    for D in range(2, dmax + 1):
        index = {m: i for i, m in enumerate(mm for mm in monomials if popcount(mm) <= D)}
        multipliers = [m for m in monomials if popcount(m) <= D]
        rows = []
        for eq in equations:
            edeg = poly_degree(eq)
            for mu in multipliers:
                if popcount(mu) + edeg > D:
                    continue
                row = 0
                ok = True
                for mono in eq:
                    prod = mu | mono
                    if popcount(prod) > D:
                        ok = False
                        break
                    row ^^= 1 << index[prod]
                if ok and row:
                    rows.append(row)
        rank, const_in_span = gf2_rank(rows, index[0])
        profile[D] = rank
        if shapes is not None:
            shapes[D] = (len(rows), len(index))
        if const_in_span:
            return D, profile
    return None, profile


def gf2_reduce(vec, basis):
    """Reduce vec by a basis keyed on each vector's lowest set bit. Every
    basis vector has no bits below its pivot, so clearing the lowest bit
    only touches higher bits and the loop terminates fully reduced."""
    while vec:
        piv = (vec & -vec).bit_length() - 1
        if piv not in basis:
            return vec, piv
        vec ^^= basis[piv]
    return 0, None


def gf2_rank(rows, const_index):
    """Row rank over GF(2) and whether the unit vector e_{const_index} is in
    the row space (i.e. the constant polynomial 1 is derivable)."""
    basis = {}
    for row in rows:
        r, piv = gf2_reduce(row, basis)
        if r:
            basis[piv] = r
    target, _ = gf2_reduce(1 << const_index, basis)
    return len(basis), target == 0


# ── factor base structures: F_2-subspaces of x-coordinates ─────────────

def subspace_from_generators(field, generators):
    """F_2-span of the given field elements, returned as an independent
    basis of coordinate ints (deduplicated by Gaussian elimination)."""
    basis = {}
    for g in generators:
        r, piv = gf2_reduce(int(field.to_int(g)), basis)
        if r:
            basis[piv] = r
    return [basis[p] for p in sorted(basis)]


def normal_conjugates(field, seed, pool=64):
    """theta, theta^2, theta^4, ... for a normal element theta chosen by
    seed among the first `pool` normal elements. A normal basis exists for
    every n, so this works at n = 131 too."""
    n = field.n
    found = []
    for v in range(1, 2**n):
        theta = field.from_int(v)
        conj, cur = [], theta
        for _ in range(n):
            conj.append(cur)
            cur = cur**2
        if len(subspace_from_generators(field, conj)) == n:
            found.append(conj)
            if len(found) >= pool:
                break
    if not found:
        raise RuntimeError("no normal element found")
    return rng_for(seed).choice(found)


def structure_bases(field, name, dim, seed):
    """(basis1, basis2) for the two summands, or (None, reason).

    nb_window puts both summands in V1 = <theta^(2^j) : j < dim>.
    nb_break is GGMP section 3.1: P1 in F1 from V1, P2 in F2 = pi^dim(F1)
    from V2 = <theta^(2^j) : dim <= j < 2 dim>. A relation P1 + P2 is
    rewritten over F1 as P1 + [lambda^dim] P2' with P2 = pi^dim(P2'), which
    recovers the 2! symmetry factor without enlarging F1."""
    if name in ("nb_window", "nb_break"):
        if 2 * dim > field.n:
            return None, "nb windows need 2*dim <= n (%d > %d)" % (2 * dim, field.n)
        conj = normal_conjugates(field, seed)
        v1 = subspace_from_generators(field, conj[:dim])
        if name == "nb_window":
            return (v1, v1), None
        v2 = subspace_from_generators(field, conj[dim:2 * dim])
        return (v1, v2), None
    basis, err = structure_basis(field, name, dim, seed)
    if basis is None:
        return None, err
    return (basis, basis), None


def structure_basis(field, name, dim, seed):
    K, z, n = field.K, field.z, field.n
    if name == "standard":
        gens = [z**i for i in range(dim)]
        return subspace_from_generators(field, gens), None
    if name == "subfield":
        if n % dim != 0:
            return None, "subfield needs dim | n (%d does not divide %d)" % (dim, n)
        # z is primitive (Conway modulus), so w has order exactly 2^dim - 1
        # and generates F_{2^dim}; its powers 0..dim-1 are an F_2-basis.
        w = z**((2**n - 1) // (2**dim - 1))
        gens = [w**i for i in range(dim)]
        return subspace_from_generators(field, gens), None
    if name == "frobenius":
        # A Frobenius-stable subspace is the squaring-closure of a seed.
        # Truncating a larger closure would break stability, so accept only
        # seeds whose closure has dimension exactly dim. These exist iff
        # x^n - 1 has a degree-dim divisor over F_2.
        candidates = []
        for v in range(1, 2**n):
            gens, cur = [], field.from_int(v)
            for _ in range(n):
                gens.append(cur)
                cur = cur**2
            basis = subspace_from_generators(field, gens)
            if len(basis) == dim:
                candidates.append(tuple(basis))
        if not candidates:
            return None, ("no Frobenius-stable subspace of dim %d: x^%d - 1 has "
                          "no degree-%d divisor over F_2" % (dim, n, dim))
        distinct = sorted(set(candidates))
        return list(rng_for(seed).choice(distinct)), None
    if name == "random":
        rng = struct_rng_for(seed, ("bases", field.n, name, dim))
        for _ in range(256):
            gens = [field.from_int(rng.randrange(1, 2**n)) for _ in range(dim)]
            basis = subspace_from_generators(field, gens)
            if len(basis) == dim:
                return basis, None
        return None, "could not draw %d independent random generators" % dim
    return None, "unknown structure %r" % name


# ── curve-side labelling and coverage (the factor base is real) ────────

def koblitz_points_with_x_in(field, basis_ints):
    """Affine points (x, y) on y^2 + x*y = x^3 + 1 with x in span(basis)."""
    K = field.K
    E = EllipticCurve(K, [1, 0, 0, 0, 1])
    xs = set()
    dim = len(basis_ints)
    for bits in range(1 << dim):
        v = 0
        for i in range(dim):
            if (bits >> i) & 1:
                v ^^= basis_ints[i]
        xs.add(v)
    points = []
    for v in xs:
        if v == 0:
            continue
        x = field.from_int(v)
        for P in E.lift_x(x, all=True):
            points.append(P)
    return E, points


def span_members(basis_ints):
    members = []
    dim = len(basis_ints)
    for bits in range(1 << dim):
        v = 0
        for i in range(dim):
            if (bits >> i) & 1:
                v ^^= basis_ints[i]
        members.append(v)
    return members


def reachable_targets_m2(field, basis1, basis2, b_int):
    """Exact reachable set of the m=2 label: every t with a decomposition
    over span(basis1) x span(basis2). For a target t, decomposability is
    the condition s = (X1+X2)^2 t^2 + X1 X2 t + (X1 X2)^2 + b == 0 for some
    pair; solving the same equation for t per pair enumerates the set, so
    the yield is measured exactly where sampling cannot see it."""
    K = field.K
    T = PolynomialRing(K, "T").gen()
    members1 = [field.from_int(v) for v in span_members(basis1)]
    members2 = [field.from_int(v) for v in span_members(basis2)]
    b = field.from_int(b_int)
    S = set()
    for X1 in members1:
        for X2 in members2:
            c2 = (X1 + X2)**2
            c1 = X1 * X2
            c0 = (X1 * X2)**2 + b
            for r, _ in (c2 * T**2 + c1 * T + c0).roots():
                S.add(field.to_int(r))
    return S


def reachable_targets_m3(field, basis1, basis2, basis3, b_int):
    """Exact reachable set of the chained m=3 label (label_targets_m3):
    every target t with a 3-decomposition over the three spans."""
    K = field.K
    T = PolynomialRing(K, "T").gen()
    m1 = [field.from_int(v) for v in span_members(basis1)]
    m2 = [field.from_int(v) for v in span_members(basis2)]
    m3 = [field.from_int(v) for v in span_members(basis3)]
    b = field.from_int(b_int)
    roots_u = []
    for x1 in m1:
        for x2 in m2:
            poly = (x1 + x2)**2 * T**2 + x1 * x2 * T + (x1 * x2)**2 + b
            roots_u.extend(r for r, _ in poly.roots())
    roots_u = list(set(roots_u))
    S = set()
    for u in roots_u:
        for x3 in m3:
            p = (u + x3)**2 * T**2 + u * x3 * T + (u * x3)**2 + b
            for r, _ in p.roots():
                S.add(field.to_int(r))
    return S


def label_targets(field, basis1, basis2, targets_int, b_int):
    """For each target x-coordinate, True if it has no decomposition with
    x(P1) in span(basis1) and x(P2) in span(basis2) (unsatisfiable system),
    which is the case D* measures."""
    members1 = [field.from_int(v) for v in span_members(basis1)]
    members2 = [field.from_int(v) for v in span_members(basis2)]
    b = field.from_int(b_int)
    labels = {}
    # The descent ranges over all of V, including x = 0 (the 2-torsion point
    # (0, 1) is a valid abscissa here), so the label covers the same domain.
    for t_int in targets_int:
        x3 = field.from_int(t_int)
        x3sq = x3 * x3
        unsat = True
        for X1 in members1:
            for X2 in members2:
                s = (X1 + X2)**2 * x3sq + X1 * X2 * x3 + (X1 * X2)**2 + b
                if s == 0:
                    unsat = False
                    break
            if not unsat:
                break
        labels[t_int] = unsat
    return labels


def label_targets_m3(field, basis1, basis2, basis3, targets_int, b_int):
    """True per target if the chained 3-point system has no common zero.
    For fixed x1, x2 the first S_3 is a quadratic in U (never identically
    zero since b = 1), so its at most two roots stand in for enumerating U
    over the whole field; this labels the descended system exactly."""
    K = field.K
    T = PolynomialRing(K, "T").gen()
    m1 = [field.from_int(v) for v in span_members(basis1)]
    m2 = [field.from_int(v) for v in span_members(basis2)]
    m3 = [field.from_int(v) for v in span_members(basis3)]
    b = field.from_int(b_int)
    roots_u = []
    for x1 in m1:
        for x2 in m2:
            poly = (x1 + x2)**2 * T**2 + x1 * x2 * T + (x1 * x2)**2 + b
            roots_u.extend(r for r, _ in poly.roots())
    roots_u = list(set(roots_u))
    labels = {}
    for t_int in targets_int:
        xr = field.from_int(t_int)
        xr2 = xr * xr
        unsat = True
        for u in roots_u:
            for x3 in m3:
                if (u + x3)**2 * xr2 + u * x3 * xr + (u * x3)**2 + b == 0:
                    unsat = False
                    break
            if not unsat:
                break
        labels[t_int] = unsat
    return labels


# ── self-test: the measurement is only as good as the descent ──────────

def eval_poly(poly, assignment_mask):
    """Value in GF(2) of a multilinear poly at a Boolean assignment: a
    monomial is 1 iff all its variables are set."""
    value = 0
    for m in poly:
        if (assignment_mask & m) == m:
            value ^^= 1
    return value


def selftest(field, dim, dmax, seed, samples):
    """Two independent checks per structure and target:
      1. the descended system has a Boolean zero iff brute force finds a
         2-decomposition over V (validates descent_equations);
      2. a satisfiable target never derives 1 (no false refutation)."""
    b_int = field.one_int
    failures = 0
    checked = 0
    rng = rng_for(seed)
    for name in ("standard", "random", "frobenius", "subfield", "nb_window", "nb_break"):
        pair, err = structure_bases(field, name, dim, seed)
        if pair is None:
            continue
        basis1, basis2 = pair
        nvars = len(basis1) + len(basis2)
        targets = [rng.randrange(1, 2**field.n) for _ in range(samples)]
        labels = label_targets(field, basis1, basis2, targets, b_int)
        for t in targets:
            eqs = descent_equations(field, basis1, basis2, t, b_int)
            has_zero = any(all(eval_poly(eq, a) == 0 for eq in eqs)
                           for a in range(1 << nvars))
            if has_zero == labels[t]:
                failures += 1
                sys.stderr.write("DESCENT MISMATCH %s target 0x%x: zero=%s unsat=%s\n"
                                 % (name, t, has_zero, labels[t]))
            if not labels[t]:
                dstar, _ = refutation_degree(eqs, nvars, dmax)
                if dstar is not None:
                    failures += 1
                    sys.stderr.write("FALSE REFUTATION %s target 0x%x at D=%d\n"
                                     % (name, t, dstar))
            checked += 1
    sys.stderr.write("selftest: %d cases, %d failures\n" % (checked, failures))
    return failures == 0


def verify_frobenius_windows(n, dim, m, seed, tries=400):
    """Check the GGMP section 3.1 mechanism: for a normal element theta, the
    windows V_i = <theta^(2^j)> over consecutive blocks of dim satisfy
    pi^dim(V_i) = V_{i+1} and are disjoint, so a summand from F_i folds into
    F_1 by a power of lambda. This gives the m! saving, not the n-fold one,
    which needs each base to be Frobenius-invariant. Uses a random normal
    element so it runs at n = 131."""
    field = BinaryField(n)
    rng = struct_rng_for(seed, ("normal", n))
    conj = None
    for _ in range(tries):
        theta = field.from_int(rng.randrange(1, 2**n))
        cand, cur = [], theta
        for _ in range(n):
            cand.append(cur)
            cur = cur**2
        if len(subspace_from_generators(field, cand)) == n:
            conj = cand
            break
    if conj is None:
        sys.stderr.write("no normal element in %d random draws\n" % tries)
        sys.exit(int(1))
    if m * dim > n:
        sys.stderr.write("need m*dim <= n\n")
        sys.exit(int(1))
    windows = [conj[i * dim:(i + 1) * dim] for i in range(m)]

    def rank(gens):
        return len(subspace_from_generators(field, gens))

    shift = 2**dim
    permuted = all(
        rank(windows[i + 1]) == dim
        and rank([g**shift for g in windows[i]] + windows[i + 1]) == dim
        for i in range(m - 1))
    disjoint = rank([g for w in windows for g in w]) == m * dim
    sys.stderr.write("n=%d dim=%d m=%d: normal element found; pi^dim(V_i) = V_(i+1): %s; "
                     "windows disjoint: %s\n" % (n, dim, m, permuted, disjoint))
    sys.exit(int(0 if permuted and disjoint else 1))


def stacking_report(n, dim, weights, seed):
    """Can symmetry breaking and the n-fold invariance saving stack where no
    invariant subspace exists? Compare, for m = 2, four bases by the repo's
    relation-collection metric trials = (unknowns + 1) / coverage:

      win_sym   both summands in a normal-basis window V1 (no orbit saving)
      win_break P1 in V1, P2 in pi^dim(V1), folded by lambda (GGMP 3.1)
      wt_sym    both summands in W_a u W_b, unknowns = Frobenius orbits
      wt_break  P1 in W_a, P2 in W_b with a != b: classes are invariant and
                disjoint, so the orbit saving and the broken symmetry coexist

    W_w = {x : normal-basis Hamming weight of x is w}. Frobenius rotates the
    normal-basis coordinates, so each W_w is invariant and orbits are
    rotation classes. Coverage counts target abscissae over the whole field.
    pairs_per_target is an enumeration proxy for per-trial cost; a SAT or
    Groebner oracle costs differently."""
    field = BinaryField(n)
    K = field.K
    T = PolynomialRing(K, "T").gen()
    b = K(1)
    conj = normal_conjugates(field, seed)
    full = (1 << n) - 1

    def element(mask):
        acc = K(0)
        for j in range(n):
            if (mask >> j) & 1:
                acc += conj[j]
        return acc

    def valid(x):
        return x != 0 and (x + x**(-2)).trace() == 0

    def rot_class(mask):
        best, cur = mask, mask
        for _ in range(n - 1):
            cur = ((cur << 1) | (cur >> (n - 1))) & full
            best = min(best, cur)
        return best

    def masks_of_weight(w):
        out = []
        for comb in itertools.combinations(range(n), w):
            m = 0
            for j in comb:
                m |= 1 << j
            out.append(m)
        return out

    def window(i):
        masks = []
        for bits in range(1, 1 << dim):
            m = 0
            for j in range(dim):
                if (bits >> j) & 1:
                    m |= 1 << (i * dim + j)
            masks.append(m)
        return masks

    def coverage(pairs):
        reach = set()
        for x1, x2 in pairs:
            poly = (x1 + x2)**2 * T**2 + x1 * x2 * T + (x1 * x2)**2 + b
            for r, _ in poly.roots():
                if r != 0:
                    reach.add(r)
        return len(reach) / (2**n - 1)

    def scheme(name, side1, side2, unknown_masks, use_orbits, symmetric):
        xs1 = [element(m) for m in side1 if valid(element(m))]
        xs2 = [element(m) for m in side2 if valid(element(m))]
        if symmetric:
            pairs = [(xs1[i], xs1[j]) for i in range(len(xs1)) for j in range(i, len(xs1))]
        else:
            pairs = [(a, c) for a in xs1 for c in xs2]
        valid_masks = [m for m in unknown_masks if valid(element(m))]
        unknowns = len(set(rot_class(m) for m in valid_masks)) if use_orbits else len(valid_masks)
        cov = coverage(pairs)
        trials = (unknowns + 1) / cov if cov else None
        return {"scheme": name, "base_abscissae": len(valid_masks), "unknowns": unknowns,
                "pairs_per_target": len(pairs), "coverage": float(cov),
                "trials": float(trials) if trials else None,
                "enumeration_work": float(trials * len(pairs)) if trials else None}

    wa, wb = weights
    v1, v2 = window(0), window(1)
    a_masks, b_masks = masks_of_weight(wa), masks_of_weight(wb)
    rows = [
        scheme("win_sym", v1, v1, v1, False, True),
        scheme("win_break", v1, v2, v1, False, False),
        scheme("wt_sym", a_masks + b_masks, a_masks + b_masks, a_masks + b_masks, True, True),
        scheme("wt_break", a_masks, b_masks, a_masks + b_masks, True, False),
    ]
    sys.stderr.write("n=%d windows dim %d, weight classes %d and %d\n" % (n, dim, wa, wb))
    sys.stderr.write("  %-10s %9s %9s %8s %9s %11s %14s\n" % (
        "scheme", "abscissae", "unknowns", "pairs", "coverage", "trials", "enum work"))
    for r in rows:
        sys.stderr.write("  %-10s %9d %9d %8d %9.4f %11s %14s\n" % (
            r["scheme"], r["base_abscissae"], r["unknowns"], r["pairs_per_target"], r["coverage"],
            "%.1f" % r["trials"] if r["trials"] else "-",
            "%.0f" % r["enumeration_work"] if r["enumeration_work"] else "-"))
    return rows


# ── three summands ─────────────────────────────────────────────────────

def structure_bases3(field, name, dim, seed):
    """(basis1, basis2, basis3), or (None, reason). nb_break uses three
    consecutive normal-basis windows V1, pi^dim(V1), pi^(2 dim)(V1), the
    m = 3 symmetry breaking whose expected probability gain is 3! = 6."""
    if name in ("nb_window", "nb_break"):
        if 3 * dim > field.n:
            return None, "3 nb windows need 3*dim <= n (%d > %d)" % (3 * dim, field.n)
        conj = normal_conjugates(field, seed)
        v1 = subspace_from_generators(field, conj[:dim])
        if name == "nb_window":
            return (v1, v1, v1), None
        v2 = subspace_from_generators(field, conj[dim:2 * dim])
        v3 = subspace_from_generators(field, conj[2 * dim:3 * dim])
        return (v1, v2, v3), None
    basis, err = structure_basis(field, name, dim, seed)
    if basis is None:
        return None, err
    return (basis, basis, basis), None


def selftest_m3(field, dim, dmax, seed, samples):
    """Same two checks as selftest, for the chained 3-point system. Brute
    force over every Boolean assignment, so keep 3*dim + n small."""
    b_int = field.one_int
    failures = 0
    checked = 0
    rng = rng_for(seed)
    for name in ("standard", "random", "subfield", "nb_window", "nb_break"):
        triple, err = structure_bases3(field, name, dim, seed)
        if triple is None:
            continue
        targets = [rng.randrange(1, 2**field.n) for _ in range(samples)]
        labels = label_targets_m3(field, *triple, targets, b_int)
        for t in targets:
            eqs, nvars = descent_equations_m3(field, *triple, t, b_int)
            if nvars > 18:
                sys.stderr.write("selftest_m3: %d vars too many to brute force\n" % nvars)
                return False
            has_zero = any(all(eval_poly(eq, a) == 0 for eq in eqs)
                           for a in range(1 << nvars))
            if has_zero == labels[t]:
                failures += 1
                sys.stderr.write("M3 DESCENT MISMATCH %s target 0x%x: zero=%s unsat=%s\n"
                                 % (name, t, has_zero, labels[t]))
            if not labels[t]:
                dstar, _ = refutation_degree(eqs, nvars, dmax)
                if dstar is not None:
                    failures += 1
                    sys.stderr.write("M3 FALSE REFUTATION %s target 0x%x at D=%d\n"
                                     % (name, t, dstar))
            checked += 1
    sys.stderr.write("selftest_m3: %d cases, %d failures\n" % (checked, failures))
    return failures == 0


def measure_structure_m3(field, name, dim, dmax, targets_int, seed, dstar_sample=None):
    triple, err = structure_bases3(field, name, dim, seed)
    if triple is None:
        return {"structure": name, "skipped": err}
    b_int = field.one_int
    labels = label_targets_m3(field, *triple, targets_int, b_int)
    decomposable = sum(1 for t in targets_int if not labels[t])
    results = []
    nvars = None
    for t_int in targets_int:
        if not labels[t_int]:
            continue
        if dstar_sample is not None and len(results) >= dstar_sample:
            break
        eqs, nvars = descent_equations_m3(field, *triple, t_int, b_int)
        dstar, profile = refutation_degree(eqs, nvars, dmax)
        results.append({"target": "0x%x" % t_int, "d_star": dstar,
                        "rank_profile": {str(k): v for k, v in profile.items()}})
    solved = [r["d_star"] for r in results if r["d_star"] is not None]
    summary = {
        "structure": name, "dim": dim, "summands": 3, "nvars": nvars,
        "decomposable_fraction": decomposable / len(targets_int),
        "unsat_targets": len(results),
        "censored_above_dmax": sum(1 for r in results if r["d_star"] is None),
    }
    if solved:
        summary["d_star_min"] = min(solved)
        summary["d_star_max"] = max(solved)
        summary["d_star_mean"] = sum(solved) / len(solved)
    summary["targets"] = results
    return summary


# ── direct symmetric S_4 (m = 3): Galbraith-Gebregiyorgis special subspaces ──
#
# The chained system above carries an n-bit intermediate U and is cubic.
# Galbraith-Gebregiyorgis (ePrint 2014/806) instead write S_4(x1,x2,x3,x_R)
# in the elementary symmetric functions e1, e2, e3 of the summands.  For
# y^2 + xy = x^3 + b that is (src/cryptanalysis/koblitz_sym_s4.rs, checked
# against the resultant):
#
#   S_4 = x^4 b^2 + e3 x^3 b + e3^2 (x^4 + b) + e3^3 x + e3^4 + e2^2 x^2 b
#       + e2^2 e3 x^3 + e2^4 x^4 + e1^2 e3 x b + e1^2 e3^2 x^2 + e1^4 b^2
#
# which is Boolean-quadratic after descent.  The unknowns are the
# coordinates of e_k in a space W_k that contains it.  For a generic V the
# products fill span(V.V) (dim up to l(l+1)/2) and span(V.V.V); for the
# "special subspace" V = theta<1, g, ..., g^(l-1)> they stay inside
# theta^k <1, ..., g^(k(l-1))>, dimension k(l-1)+1.  Fewer unknowns at the
# same base size is the lever; this measures what it does to D*.

def product_span(field, basis_a, basis_b):
    gens = []
    for a in basis_a:
        fa = field.from_int(a)
        for b in basis_b:
            gens.append(fa * field.from_int(b))
    return subspace_from_generators(field, gens)


def symmetric_spaces(field, name, dim, seed):
    """Bases (W1, W2, W3) holding e1, e2, e3 for x_i in V = span(W1).

    geometric: V = theta<1, g, ..., g^(dim-1)>, W_k = theta^k <1..g^(k(dim-1))>.
    standard:  theta = 1, g = z (the polynomial-basis window).
    subfield:  V = F_{2^dim}, closed under products, so W_k = V.
    random:    V random, W2 = span(V.V), W3 = span(V.V.V)."""
    K, z, n = field.K, field.z, field.n
    rng = struct_rng_for(seed, ("sym", field.n, name, dim))
    if name in ("geometric", "standard"):
        if name == "standard":
            theta, g = K(1), z
        else:
            theta = field.from_int(rng.randrange(1, 2**n))
            g = field.from_int(rng.randrange(2, 2**n))
        spaces = []
        for k in (1, 2, 3):
            gens = [theta**k * g**j for j in range(k * (dim - 1) + 1)]
            spaces.append(subspace_from_generators(field, gens))
        if len(spaces[0]) != dim:
            return None, "geometric generators dependent at dim %d" % dim
        return tuple(spaces), None
    if name == "subfield":
        v, err = structure_basis(field, "subfield", dim, seed)
        if v is None:
            return None, err
        return (v, v, v), None
    if name == "random":
        v, err = structure_basis(field, "random", dim, seed)
        if v is None:
            return None, err
        w2 = product_span(field, v, v)
        w3 = product_span(field, w2, v)
        return (v, w2, w3), None
    return None, "unknown symmetric structure %r" % name


def descent_equations_sym_s4(field, w1, w2, w3, target_int, b_int):
    """Weil descent of the symmetrised S_4 with e_k = sum of W_k coordinates
    times fresh variables. Returns (equations, nvars)."""
    d1, d2, d3 = len(w1), len(w2), len(w3)
    E1 = fe_linear(field, w1, 0)
    E2 = fe_linear(field, w2, d1)
    E3 = fe_linear(field, w3, d1 + d2)
    x = field.from_int(target_int)
    b = field.from_int(b_int)

    def c(v):
        return fe_const(field, field.to_int(v))

    def mul(u, v):
        return fe_mul(field, u, v)

    def sq(u):
        return fe_square(field, u)

    E1_2, E2_2, E3_2 = sq(E1), sq(E2), sq(E3)
    terms = [
        c(x**4 * b**2),
        mul(E3, c(x**3 * b)),
        mul(E3_2, c(x**4 + b)),
        mul(mul(E3_2, E3), c(x)),
        sq(E3_2),
        mul(E2_2, c(x**2 * b)),
        mul(mul(E2_2, E3), c(x**3)),
        mul(sq(E2_2), c(x**4)),
        mul(mul(E1_2, E3), c(x * b)),
        mul(mul(E1_2, E3_2), c(x**2)),
        mul(sq(E1_2), c(b**2)),
    ]
    acc = fe_zero(field.n)
    for t in terms:
        acc = fe_add(acc, t)
    return [eq for eq in acc if eq], d1 + d2 + d3


def boolean_unsat(equations, nvars):
    """Exact satisfiability of a multilinear GF(2) system through PolyBoRi:
    the reduced Groebner basis is [1] iff there is no Boolean solution."""
    B = BooleanPolynomialRing(max(int(nvars), 1), "v")
    gens = B.gens()
    polys = []
    for eq in equations:
        p = B(0)
        for m in eq:
            term = B(1)
            k, bits = 0, int(m)
            while bits:
                if bits & 1:
                    term *= gens[k]
                bits >>= 1
                k += 1
            p += term
        polys.append(p)
    gb = B.ideal(polys).groebner_basis()
    return list(gb) == [B(1)]


def selftest_sym_s4(field, dim, seed, samples):
    """The symmetric system must be satisfiable whenever the target is a
    true 3-sum over V (completeness); spurious Boolean solutions (e's not
    coming from any x_i in V) are allowed and counted, not failures."""
    b_int = field.one_int
    rng = rng_for(seed)
    failures = checked = 0
    for name in ("standard", "geometric", "subfield", "random"):
        spaces, err = symmetric_spaces(field, name, dim, seed)
        if spaces is None:
            continue
        v = spaces[0]
        targets = [rng.randrange(1, 2**field.n) for _ in range(samples)]
        point_unsat = label_targets_m3(field, v, v, v, targets, b_int)
        for t in targets:
            eqs, nvars = descent_equations_sym_s4(field, *spaces, t, b_int)
            unsat = boolean_unsat(eqs, nvars)
            if unsat and not point_unsat[t]:
                failures += 1
                sys.stderr.write("SYM S4 INCOMPLETE %s target 0x%x: decomposable but refuted\n"
                                 % (name, t))
            checked += 1
    sys.stderr.write("selftest_sym_s4: %d cases, %d failures\n" % (checked, failures))
    return failures == 0


def measure_symmetric_s4(field, name, dim, dmax, targets_int, seed, dstar_sample=None):
    spaces, err = symmetric_spaces(field, name, dim, seed)
    if spaces is None:
        return {"structure": name, "skipped": err}
    b_int = field.one_int
    v = spaces[0]
    point_unsat = label_targets_m3(field, v, v, v, targets_int, b_int)
    results = []
    bool_unsat_count = spurious = 0
    nvars = None
    for t_int in targets_int:
        eqs, nvars = descent_equations_sym_s4(field, *spaces, t_int, b_int)
        unsat = boolean_unsat(eqs, nvars)
        if not unsat:
            if point_unsat[t_int]:
                spurious += 1
            continue
        bool_unsat_count += 1
        if dstar_sample is not None and len(results) >= dstar_sample:
            continue
        shapes = {}
        dstar, profile = refutation_degree(eqs, nvars, dmax, shapes)
        results.append({"target": "0x%x" % t_int, "d_star": dstar,
                        "rank_profile": {str(k): val for k, val in profile.items()},
                        "macaulay_shape": {str(k): list(val) for k, val in shapes.items()}})
    solved = [r["d_star"] for r in results if r["d_star"] is not None]
    summary = {
        "structure": name, "dim": dim, "summands": 3, "encoding": "symmetric_s4",
        "space_dims": [len(w) for w in spaces], "nvars": nvars,
        "chain_nvars_for_comparison": 3 * dim + field.n,
        "equation_degree": max((poly_degree(eq) for eq in descent_equations_sym_s4(
            field, *spaces, targets_int[0], b_int)[0]), default=0),
        "point_decomposable_fraction": sum(1 for t in targets_int if not point_unsat[t]) / len(targets_int),
        "decomposable_fraction": sum(1 for t in targets_int if not point_unsat[t]) / len(targets_int),
        "reachable_set_size": len(reachable_targets_m3(field, v, v, v, b_int)),
        "exact_yield": len(reachable_targets_m3(field, v, v, v, b_int)) / float(2**field.n),
        "boolean_unsat_fraction": bool_unsat_count / len(targets_int),
        "spurious_satisfiable": spurious,
        "unsat_targets": len(results),
        "censored_above_dmax": sum(1 for r in results if r["d_star"] is None),
    }
    if solved:
        summary["d_star_min"] = min(solved)
        summary["d_star_max"] = max(solved)
        summary["d_star_mean"] = sum(solved) / len(solved)
    summary["targets"] = results
    return summary


# ── Shantz-Teske delta method: m*dim = n - delta ───────────────────────
#
# Shantz & Teske (2013) pick the base dimension so the descended system is
# over-determined: m*l = n - delta with delta > 0 lowers the decomposition
# probability by about 2^-delta but each system is smaller and (they
# report) solves at a lower degree.  The sweep reports both sides and the
# product that matters: Macaulay work at D* divided by the decomposable
# fraction, i.e. a cost-per-relation proxy.

def delta_sweep(field, name, summands, deltas, dmax, targets_int, seed, dstar_sample=None):
    n = field.n
    rows = []
    for delta in deltas:
        if (n - delta) % summands != 0:
            rows.append({"delta": delta, "skipped": "m*dim = n - delta has no integer dim"})
            continue
        dim = (n - delta) // summands
        if dim < 1:
            rows.append({"delta": delta, "skipped": "dim < 1"})
            continue
        if summands == 2:
            pair, err = structure_bases(field, name, dim, seed)
            if pair is None:
                rows.append({"delta": delta, "dim": dim, "skipped": err})
                continue
            basis1, basis2 = pair
            nvars = len(basis1) + len(basis2)
            labels = label_targets(field, basis1, basis2, targets_int, field.one_int)
            build = lambda t: (descent_equations(field, basis1, basis2, t, field.one_int), nvars)
        else:
            triple, err = structure_bases3(field, name, dim, seed)
            if triple is None:
                rows.append({"delta": delta, "dim": dim, "skipped": err})
                continue
            labels = label_targets_m3(field, *triple, targets_int, field.one_int)
            build = lambda t: descent_equations_m3(field, *triple, t, field.one_int)
        decomposable = sum(1 for t in targets_int if not labels[t])
        dstars, shapes_at_dstar = [], []
        nv = None
        censored = 0
        for t in targets_int:
            if not labels[t]:
                continue
            if dstar_sample is not None and len(dstars) + censored >= dstar_sample:
                break
            eqs, nv = build(t)
            shapes = {}
            dstar, _ = refutation_degree(eqs, nv, dmax, shapes)
            if dstar is not None:
                dstars.append(dstar)
                shapes_at_dstar.append(shapes[dstar])
            else:
                censored += 1
        row = {"delta": delta, "dim": dim, "nvars": nv,
               "decomposable_fraction": decomposable / len(targets_int),
               "unsat_measured": len(dstars), "censored_above_dmax": censored}
        if dstars:
            row["d_star_mean"] = sum(dstars) / len(dstars)
            row["d_star_max"] = max(dstars)
            mean_rows = sum(s[0] for s in shapes_at_dstar) / len(shapes_at_dstar)
            mean_cols = sum(s[1] for s in shapes_at_dstar) / len(shapes_at_dstar)
            row["macaulay_rows_mean"] = mean_rows
            row["macaulay_cols_mean"] = mean_cols
            row["work_proxy"] = mean_rows * mean_cols
            if decomposable:
                row["work_per_relation_proxy"] = row["work_proxy"] / (decomposable / len(targets_int))
        rows.append(row)
    base = next((r for r in rows if r.get("delta") == 0 and "decomposable_fraction" in r), None)
    for r in rows:
        if base and "decomposable_fraction" in r and base["decomposable_fraction"]:
            r["fraction_ratio_to_delta0"] = r["decomposable_fraction"] / base["decomposable_fraction"]
            r["predicted_ratio_2^-delta"] = 2.0 ** (-r["delta"])
    return rows


# ── driver ─────────────────────────────────────────────────────────────

def measure_structure(field, name, dim, dmax, targets_int, seed):
    pair, err = structure_bases(field, name, dim, seed)
    if pair is None:
        return {"structure": name, "skipped": err}
    basis1, basis2 = pair
    nvars = len(basis1) + len(basis2)
    b_int = field.one_int
    labels = label_targets(field, basis1, basis2, targets_int, b_int)
    E, points = koblitz_points_with_x_in(field, basis1)
    decomposable = sum(1 for t in targets_int if not labels[t])
    results = []
    for t_int in targets_int:
        if not labels[t_int]:
            continue  # satisfiable targets do not refute; D* is for UNSAT
        eqs = descent_equations(field, basis1, basis2, t_int, b_int)
        shapes = {}
        dstar, profile = refutation_degree(eqs, nvars, dmax, shapes)
        results.append({
            "target": "0x%x" % t_int,
            "d_star": dstar,
            "rank_profile": {str(k): v for k, v in profile.items()},
            "macaulay_shape": {str(k): list(val) for k, val in shapes.items()},
        })
    solved = [r["d_star"] for r in results if r["d_star"] is not None]
    censored = sum(1 for r in results if r["d_star"] is None)
    summary = {
        "structure": name,
        "dim": dim,
        "basis": ["0x%x" % v for v in basis1],
        "basis2": None if basis2 == basis1 else ["0x%x" % v for v in basis2],
        "factor_base_points": len(points),
        "decomposable_fraction": decomposable / len(targets_int),
        "unsat_targets": len(results),
        "censored_above_dmax": censored,
        "reachable_set_size": len(reachable_targets_m2(field, basis1, basis2, b_int)),
    }
    if solved:
        solved_sorted = sorted(solved)
        summary["d_star_min"] = solved_sorted[0]
        summary["d_star_max"] = solved_sorted[-1]
        summary["d_star_mean"] = sum(solved) / len(solved)
    summary["targets"] = results
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Measure factor base -> degree of regularity (solving degree)."
    )
    parser.add_argument("--field", type=int, default=7, help="binary field degree n (odd, prime 2n+1 not required)")
    parser.add_argument("--dim", type=int, default=3, help="factor-base subspace dimension n'")
    parser.add_argument("--dmax", type=int, default=6, help="maximum Macaulay degree to search for D*")
    parser.add_argument("--targets", type=int, default=24, help="number of random target x-coordinates")
    parser.add_argument("--seed", type=int, default=1, help="seed for targets and random structures")
    parser.add_argument("--structures", default="subfield,standard,frobenius,random,nb_window,nb_break",
                        help="comma-separated factor-base structures to compare")
    parser.add_argument("--curve", help="p192/p224/p256 print the prime-field DoR note and exit")
    parser.add_argument("--stacking", metavar="A,B",
                        help="compare window vs invariant weight-class bases W_A, W_B for m=2 by "
                             "trials = (unknowns+1)/coverage, then exit")
    parser.add_argument("--frobenius-windows", action="store_true",
                        help="check that Frobenius maps each symmetry-breaking window to the next "
                             "(the GGMP section 3.1 fold), then exit; works at n=131")
    parser.add_argument("--out", help="write the JSON report here; default stdout")
    parser.add_argument("--selftest", type=int, default=0, metavar="N",
                        help="validate descent and refutation on N targets per structure, then exit")
    parser.add_argument("--summands", type=int, default=2, choices=(2, 3),
                        help="points per decomposition; 3 uses the chained Semaev system")
    parser.add_argument("--dstar-sample", type=int, default=None, metavar="K",
                        help="with --summands 3, compute D* on at most K non-decomposable targets; "
                             "the decomposable fraction still uses every target")
    parser.add_argument("--symmetric-s4", action="store_true",
                        help="m = 3 through the symmetrised S_4 in e1, e2, e3 (Galbraith-"
                             "Gebregiyorgis special subspaces); structures geometric,standard,"
                             "subfield,random")
    parser.add_argument("--delta-sweep", metavar="D0,D1,...",
                        help="Shantz-Teske: sweep dim = (n - delta)/m for the listed deltas on the "
                             "first --structures entry, report D*, decomposable fraction and a "
                             "work-per-relation proxy, then exit")
    args = parser.parse_args()

    if args.selftest:
        if args.symmetric_s4:
            ok = selftest_sym_s4(BinaryField(args.field), args.dim, args.seed, args.selftest)
        else:
            check = selftest_m3 if args.summands == 3 else selftest
            ok = check(BinaryField(args.field), args.dim, args.dmax, args.seed, args.selftest)
        sys.exit(int(0 if ok else 1))

    if args.frobenius_windows:
        verify_frobenius_windows(args.field, args.dim, args.summands, args.seed)
        return

    if args.stacking:
        wa, wb = (int(w) for w in args.stacking.split(","))
        rows = stacking_report(args.field, args.dim, (wa, wb), args.seed)
        if args.out:
            with open(args.out, "w") as handle:
                handle.write(json.dumps({"schema": "factor-base-stacking-v1", "field_degree": args.field,
                                         "window_dim": args.dim, "weights": [wa, wb], "seed": int(args.seed),
                                         "rows": rows}, indent=2, default=int) + "\n")
        return

    if args.curve in ("p192", "p224", "p256"):
        sys.stderr.write(
            "%s is a prime field. Semaev index calculus has no subfield to Weil-descend\n"
            "into, so the degree of regularity is a different object: either Diem's\n"
            "descent over a degree-k extension F_{p^k}, or the solving degree of the\n"
            "symmetric Semaev system in m variables over F_p. This engine measures the\n"
            "binary-field DoR (the ECC2K-130 analog). See the README for the prime-field\n"
            "ideas. Use --field/--dim to run the binary measurement.\n" % args.curve)
        return

    field = BinaryField(args.field)
    rng = rng_for(args.seed)
    seen = set()
    targets_int = []
    while len(targets_int) < args.targets and len(seen) < 2**args.field - 1:
        t = rng.randrange(1, 2**args.field)
        if t in seen:
            continue
        seen.add(t)
        targets_int.append(t)

    structures = [s.strip() for s in args.structures.split(",") if s.strip()]

    if args.delta_sweep:
        deltas = [int(d) for d in args.delta_sweep.split(",")]
        name = structures[0]
        rows = delta_sweep(field, name, args.summands, deltas, args.dmax, targets_int,
                           args.seed, args.dstar_sample)
        report = {"schema": "factor-base-delta-sweep-v1", "field_degree": args.field,
                  "structure": name, "summands": args.summands, "dmax": args.dmax,
                  "targets": len(targets_int), "rows": rows}
        text = json.dumps(report, indent=2, default=lambda o: int(o) if isinstance(o, Integer) else float(o))
        if args.out:
            with open(args.out, "w") as handle:
                handle.write(text + "\n")
        else:
            sys.stdout.write(text + "\n")
        for r in rows:
            if "skipped" in r:
                sys.stderr.write("  delta=%+d skipped: %s\n" % (r["delta"], r["skipped"]))
                continue
            sys.stderr.write("  delta=%+d dim=%d vars=%s D*=%s frac=%.3f (x%.2f vs 2^-d %.2f) work/rel=%s\n" % (
                r["delta"], r["dim"], r["nvars"], r.get("d_star_mean"), r["decomposable_fraction"],
                r.get("fraction_ratio_to_delta0", float("nan")), r.get("predicted_ratio_2^-delta", float("nan")),
                ("%.3g" % r["work_per_relation_proxy"]) if "work_per_relation_proxy" in r else None))
        return

    if args.symmetric_s4:
        if args.structures == parser.get_default("structures"):
            structures = ["geometric", "standard", "subfield", "random"]

        def measure(field, s, dim, dmax, targets_int, seed):
            return measure_symmetric_s4(field, s, dim, dmax, targets_int, seed, args.dstar_sample)
    elif args.summands == 3:
        def measure(field, s, dim, dmax, targets_int, seed):
            return measure_structure_m3(field, s, dim, dmax, targets_int, seed, args.dstar_sample)
    else:
        measure = measure_structure
    report = {
        "schema": "factor-base-dor-v1",
        "curve": "Koblitz y^2+xy=x^3+1 over F_2^%d" % args.field,
        "field_degree": args.field,
        "dim": args.dim,
        "summands": 3 if args.symmetric_s4 else args.summands,
        "encoding": "symmetric_s4" if args.symmetric_s4 else "chained_s3",
        "dmax": args.dmax,
        "targets": len(targets_int),
        "observable": "D* = solving/refutation degree on non-decomposable targets",
        "results": [measure(field, s, args.dim, args.dmax, targets_int, args.seed)
                    for s in structures],
    }
    ranked = [r for r in report["results"] if "d_star_mean" in r]
    ranked.sort(key=lambda r: r["d_star_mean"])
    report["ranking_by_mean_d_star"] = [
        {"structure": r["structure"], "d_star_mean": r["d_star_mean"],
         "d_star_min": r["d_star_min"], "d_star_max": r["d_star_max"],
         "decomposable_fraction": r["decomposable_fraction"],
         "unsat_targets": r["unsat_targets"]}
        for r in ranked
    ]

    text = json.dumps(report, indent=2, default=lambda o: int(o) if isinstance(o, Integer) else float(o))
    if args.out:
        with open(args.out, "w") as handle:
            handle.write(text + "\n")
        sys.stderr.write("wrote %s\n" % args.out)
    else:
        sys.stdout.write(text + "\n")

    for row in report["ranking_by_mean_d_star"]:
        sys.stderr.write("  %-10s mean D*=%.3f  [%d..%d]  decomposable %.3f  over %d unsat targets\n" % (
            row["structure"], row["d_star_mean"], row["d_star_min"],
            row["d_star_max"], row["decomposable_fraction"], row["unsat_targets"]))
    if args.symmetric_s4:
        for r in report["results"]:
            if "skipped" in r:
                sys.stderr.write("  %-10s skipped: %s\n" % (r["structure"], r["skipped"]))
            else:
                sys.stderr.write("  %-10s W dims %s vars=%d (chain %d) point-decomp %.3f bool-unsat %.3f spurious %d\n" % (
                    r["structure"], r["space_dims"], r["nvars"], r["chain_nvars_for_comparison"],
                    r["point_decomposable_fraction"], r["boolean_unsat_fraction"], r["spurious_satisfiable"]))


if __name__ == "__main__":
    main()
