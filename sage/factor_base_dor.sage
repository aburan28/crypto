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

def descent_equations(field, basis_ints, target_int, b_int):
    """The n GF(2) polynomials whose common zero is a 2-decomposition of
    the target x-coordinate over V = span(basis_ints)."""
    dim = len(basis_ints)
    X1 = fe_linear(field, basis_ints, 0)
    X2 = fe_linear(field, basis_ints, dim)
    x3 = fe_const(field, target_int)
    bvec = fe_const(field, b_int)
    sum12 = fe_add(X1, X2)
    prod12 = fe_mul(field, X1, X2)
    term1 = fe_mul(field, fe_square(field, sum12), fe_square(field, x3))
    term2 = fe_mul(field, prod12, x3)
    term3 = fe_add(fe_square(field, prod12), bvec)
    s3 = fe_add(fe_add(term1, term2), term3)
    return [eq for eq in s3 if eq]


def poly_degree(poly):
    return max((popcount(m) for m in poly), default=0)


# ── Macaulay refutation degree over GF(2) ──────────────────────────────

def refutation_degree(equations, nvars, dmax):
    """Smallest D in 2..dmax at which the constant 1 lies in the GF(2) row
    space of the degree-D Macaulay matrix (the solving degree D*), or None.
    Returns (D*, rank_profile) with the row rank at each degree reached."""
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
        rng = rng_for(seed)
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


def label_targets(field, basis_ints, targets_int, b_int):
    """For each target x-coordinate, True if it has no 2-decomposition over
    V (unsatisfiable system), which is the case D* measures."""
    K = field.K
    members = []
    dim = len(basis_ints)
    for bits in range(1 << dim):
        v = 0
        for i in range(dim):
            if (bits >> i) & 1:
                v ^^= basis_ints[i]
        members.append(v)
    b = field.from_int(b_int)
    labels = {}
    # The descent ranges over all of V, including x = 0 (the 2-torsion point
    # (0, 1) is a valid abscissa here), so the label covers the same domain.
    for t_int in targets_int:
        x3 = field.from_int(t_int)
        x3sq = x3 * x3
        unsat = True
        for v1 in members:
            X1 = field.from_int(v1)
            for v2 in members:
                X2 = field.from_int(v2)
                s = (X1 + X2)**2 * x3sq + X1 * X2 * x3 + (X1 * X2)**2 + b
                if s == 0:
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
    nvars = 2 * dim
    failures = 0
    checked = 0
    rng = rng_for(seed)
    for name in ("standard", "random", "frobenius", "subfield"):
        basis, err = structure_basis(field, name, dim, seed)
        if basis is None:
            continue
        targets = [rng.randrange(1, 2**field.n) for _ in range(samples)]
        labels = label_targets(field, basis, targets, b_int)
        for t in targets:
            eqs = descent_equations(field, basis, t, b_int)
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


# ── driver ─────────────────────────────────────────────────────────────

def measure_structure(field, name, dim, dmax, targets_int, seed):
    basis, err = structure_basis(field, name, dim, seed)
    if basis is None:
        return {"structure": name, "skipped": err}
    b_int = field.one_int
    labels = label_targets(field, basis, targets_int, b_int)
    E, points = koblitz_points_with_x_in(field, basis)
    results = []
    for t_int in targets_int:
        if not labels[t_int]:
            continue  # satisfiable targets do not refute; D* is for UNSAT
        eqs = descent_equations(field, basis, t_int, b_int)
        dstar, profile = refutation_degree(eqs, 2 * dim, dmax)
        results.append({
            "target": "0x%x" % t_int,
            "d_star": dstar,
            "rank_profile": {str(k): v for k, v in profile.items()},
        })
    solved = [r["d_star"] for r in results if r["d_star"] is not None]
    censored = sum(1 for r in results if r["d_star"] is None)
    summary = {
        "structure": name,
        "dim": dim,
        "basis": ["0x%x" % v for v in basis],
        "factor_base_points": len(points),
        "unsat_targets": len(results),
        "censored_above_dmax": censored,
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
    parser.add_argument("--structures", default="subfield,standard,frobenius,random",
                        help="comma-separated factor-base structures to compare")
    parser.add_argument("--curve", help="p192/p224/p256 print the prime-field DoR note and exit")
    parser.add_argument("--out", help="write the JSON report here; default stdout")
    parser.add_argument("--selftest", type=int, default=0, metavar="N",
                        help="validate descent and refutation on N targets per structure, then exit")
    args = parser.parse_args()

    if args.selftest:
        ok = selftest(BinaryField(args.field), args.dim, args.dmax, args.seed, args.selftest)
        sys.exit(int(0 if ok else 1))

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
    report = {
        "schema": "factor-base-dor-v1",
        "curve": "Koblitz y^2+xy=x^3+1 over F_2^%d" % args.field,
        "field_degree": args.field,
        "dim": args.dim,
        "summands": 2,
        "dmax": args.dmax,
        "targets": len(targets_int),
        "observable": "D* = solving/refutation degree on non-decomposable targets",
        "results": [measure_structure(field, s, args.dim, args.dmax, targets_int, args.seed)
                    for s in structures],
    }
    ranked = [r for r in report["results"] if "d_star_mean" in r]
    ranked.sort(key=lambda r: r["d_star_mean"])
    report["ranking_by_mean_d_star"] = [
        {"structure": r["structure"], "d_star_mean": r["d_star_mean"],
         "d_star_min": r["d_star_min"], "d_star_max": r["d_star_max"],
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
        sys.stderr.write("  %-10s mean D*=%.3f  [%d..%d] over %d unsat targets\n" % (
            row["structure"], row["d_star_mean"], row["d_star_min"],
            row["d_star_max"], row["unsat_targets"]))


if __name__ == "__main__":
    main()
