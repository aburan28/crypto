#!/usr/bin/env python3
"""Isogeny-loop endomorphism sweeper for GLV-style scalar multiplication.

Implements, as a tool, the construction of Koshelev–Sanso, "Endomorphisms for
Faster Cryptography on Elliptic Curves of Moderate CM Discriminants"
(ePrint 2024/1985; De Cifris Koine, CIFRIS25 ACTA, doi:10.69091/koine/vol-7-W09):

  An ordinary curve E/F_q with End(E) = O (an order of discriminant D < 0 in
  K = Q(sqrt(D))) sits on the crater of its isogeny volcano, on which the ideal
  class group Cl(O) acts freely and transitively by horizontal isogenies.  A
  closed non-backtracking walk of horizontal isogenies of small prime degrees
  l_1, ..., l_m from E back to E is an endomorphism phi of E of degree
  l_1 ... l_m that costs about sum_i c(l_i) field multiplications to evaluate
  (each step is a precomputed Vélu/Kohel isogeny), not c(l_1 ... l_m).  Such
  loops are exactly the vectors e in the *relation lattice*
      Lambda = { e in Z^k : prod_i [l_i]^{e_i} = 1 in Cl(O) },
  and the cheapest loop is the shortest nonzero vector of Lambda in the
  weighted 1-norm  ||e||_c = sum_i c(l_i) |e_i|  (paper, Section 4).  The
  endomorphism is a generator alpha of the principal ideal prod_i l_i^{e_i};
  it acts on the prime-order subgroup E(F_q)[r] as the scalar
  lambda = alpha mod (r, pi - 1), which is what GLV needs.

What this script does
  * `curve` mode: for a curve given by its Frobenius trace (registry slug or
    alias in docs/curves/registry.json, a built-in Certicom challenge curve, or
    explicit --p/--m/--trace/--order), compute Delta = t^2 - 4q = f^2 D_K, the
    class group of the chosen endomorphism order O (default: the maximal order
    O_K, which gives the cheapest possible loops; --conductor g pins O = O_g),
    the ideal classes of all primes l <= L that split in O, the relation
    lattice, its cheapest vectors under the weighted 1-norm (LLL on the weighted
    Gram matrix, then Fincke–Pohst enumeration of the ball that is guaranteed to
    contain every cheaper vector), the explicit endomorphism alpha = a + b*sqrt(D_K),
    its eigenvalue on E(F_q)[r], the GLV basis quality, the per-prime Frobenius
    eigenvalue that selects each kernel, and the cost against the doublings the
    GLV split saves.
  * `disc` mode: the same for a CM discriminant alone (curve-independent part).
  * `scan` mode: rank a range of fundamental discriminants by cheapest loop,
    for choosing a CM discriminant before running the CM method.
  * `selftest` mode: build a crater curve for D = -71 over a ~30-bit prime,
    walk the loops found by the lattice search with PARI's Vélu isogenies,
    select each kernel by its Frobenius eigenvalue, and check that the
    composite acts on a point of prime order r exactly as [lambda].

Requires cypari2 (PARI/GP).  Nothing here is constant-time or a production
implementation; it is a search and verification tool.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction

try:
    import cypari2
except ImportError:  # pragma: no cover
    sys.stderr.write("cypari2 is required: pip install cypari2\n")
    raise

pari = cypari2.Pari()
pari.allocatemem(1 << 31)
pari.default("parisizemax", 1 << 33)
PariError = cypari2.PariError

HERE = os.path.dirname(os.path.abspath(__file__))
REGISTRY = os.path.join(HERE, "..", "docs", "curves", "registry.json")

# --------------------------------------------------------------------------
# Built-in curves: the Certicom ECC challenge curves that are not in the
# registry (parameters from the Certicom challenge document, hex).  #E = h*n.
# Every entry was checked with PARI: a random point P on the stated curve
# satisfies [h*n]P = O, |t| <= 2 sqrt(q), and n is prime.
# --------------------------------------------------------------------------
BUILTIN = {
    "ECCp-131": dict(family="prime", p=0x048E1D43F293469E33194C43186B3ABC0B,
                     n=0x048E1D43F293469E317F7ED728F6B8E6F1, h=1,
                     a=0x041CB121CE2B31F608A76FC8F23D73CB66, b=0x02F74F717E8DEC90991E5EA9B2FF03DA58),
    "ECCp-109": dict(family="prime", p=0x1BD579792B380B5B521E6D9FB599,
                     n=0x1BD579792B380B049C4D13A75AE5, h=1,
                     a=0x0FD4C926FD178E9805E663021744, b=0x153D3CBB508FFE3A7F31FF4FAFFD),
    "ECCp-97": dict(family="prime", p=0x016EA1595ED21AE4D8D8420E35,
                    n=0x016EA1595ED21AE98FB6CCA20D, h=1,
                    a=0x0047370916A603B07657C305C4, b=0x01124DF86D04064F503D9925AF),
    "ECCp-89": dict(family="prime", p=0x0158685C903F1643908BA955,
                    n=0x0158685C903EF906D7F58D47, h=1,
                    a=0x0C8AE4F7DE8918AA9FAB2260, b=0x00647E7EA1062AE69A7D1037),
    "ECCp-79": dict(family="prime", p=0x62CE5177412ACA899CF5,
                    n=0x62CE5177407B7258DC31, h=1,
                    a=0x39C95E6DDDB1BC45733C, b=0x1F16D880E89D5A1C0ED1),
    "ECC2-131": dict(family="binary", m=131, modulus=(1 << 131) | (1 << 13) | (1 << 2) | (1 << 1) | 1,
                     n=0x0400000000000000026ABB991FE311FE83, h=2,
                     a=0x07EBCB7EECC296A1C4A1A14F2C9E44352E, b=0x00610B0A57C73649AD0093BDD622A61D81),
    "ECC2-109": dict(family="binary", m=109, modulus=(1 << 109) | (1 << 9) | (1 << 2) | (1 << 1) | 1,
                     n=0x10000000000000053701AB26100B, h=2,
                     a=0x14BAA8C4131E992C7E35FCF70CE3, b=0x1333BE219E61625E4C2B6B1032D9),
    "ECC2-97": dict(family="binary", m=97, modulus=(1 << 97) | (1 << 6) | 1,
                    n=0x10000000000007383E2DE1E81, h=2,
                    a=0x01EA5CE2B7F0A58E01B4389418, b=0x009687742B6329E70680231988),
    "ECC2-89": dict(family="binary", m=89, modulus=(1 << 89) | (1 << 38) | 1,
                    n=0x0100000000000B41C8C9D8FD, h=2,
                    a=0x0095AA3E660B75E77315E94E, b=0x01AC2701C6C54021D1BF0A72),
    "ECC2-79": dict(family="binary", m=79, modulus=(1 << 79) | (1 << 9) | 1,
                    n=0x40000000004531A2562B, h=2,
                    a=0x4A2E38A8F66D7F4C385F, b=0x2C0BB31C6BECC03D68A7),
    "ECC2K-108": dict(family="koblitz", m=109, modulus=(1 << 109) | (1 << 9) | (1 << 2) | (1 << 1) | 1,
                      n=0x0FFFFFFFFFFFFFA621B02C383E9B, h=2, a=1, b=1),
}


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def ffelt_to_int(z):
    """Integer value of a prime-field t_FFELT (or t_INTMOD / t_INT)."""
    if str(pari.type(z)) == "t_FFELT":
        return int(pari("(a)->lift(polcoef(a.pol,0))")(z))
    return int(pari.lift(z))


def mk_matrix(rows):
    """PARI t_MAT from a list of rows (lists of ints/strings)."""
    m, n = len(rows), len(rows[0])
    return pari.matrix(m, n, [pari(str(x)) for row in rows for x in row])


def primes_upto(n):
    return [int(q) for q in pari.primes(pari.primepi(n))]


def koblitz_trace(a, n):
    """Trace of y^2 + xy = x^3 + a x^2 + 1 over F_{2^n}: t_0 = 2, t_1 = mu, t_k = mu t_{k-1} - 2 t_{k-2}."""
    mu = (-1) ** (1 - a)
    t0, t1 = 2, mu
    for _ in range(n - 1):
        t0, t1 = t1, mu * t1 - 2 * t0
    return t1


def factor_with_budget(N, trial_limit=1 << 20, seconds=60):
    """Return (list of (p, e), unfactored cofactor).  The cofactor is 1 when the
    factorization is complete; otherwise it has no prime factor <= trial_limit
    and factorint did not finish within the budget."""
    N = abs(int(N))
    f = pari.factor(N, trial_limit)
    fac, cof = [], 1
    nrows = int(pari.matsize(f)[0])
    for i in range(nrows):
        q, e = int(f[i, 0]), int(f[i, 1])
        if q > trial_limit and not pari.isprime(q):
            cof *= q ** e
        else:
            fac.append((q, e))
    if cof > 1:
        try:
            g = pari(f"alarm({seconds}, factorint({cof}))")
            for i in range(int(pari.matsize(g)[0])):
                fac.append((int(g[i, 0]), int(g[i, 1])))
            cof = 1
        except PariError:
            pass
    fac.sort()
    return fac, cof


def fundamental_and_conductor(delta, fac, cof):
    """Delta = f^2 D_K with D_K fundamental, from a (possibly partial) factorization
    of |Delta|.  An unfactored cofactor is assumed squarefree (stated in output)."""
    sqf = 1
    f = 1
    for q, e in fac:
        if e % 2:
            sqf *= q
        f *= q ** (e // 2)
    sqf *= cof
    D = -sqf
    if D % 4 != 1:
        D *= 4
        assert f % 2 == 0, "parity bookkeeping"
        f //= 2
    assert f * f * D == delta, "Delta != f^2 D_K"
    return D, f


def is_fundamental(D):
    return int(pari.isfundamental(D)) == 1


# --------------------------------------------------------------------------
# Cost model (field multiplications per isogeny step)
# --------------------------------------------------------------------------
def step_cost(ell, char2, model):
    if char2 and ell == 2:
        # Frobenius (inseparable 2-isogeny) or Verschiebung: a couple of
        # squarings; charged 2 M so the weight stays positive.
        return 2.0
    if model == "projective":
        return 7.5 * ell          # paper, Section 2.2: Horner on homogenised Vélu maps
    if model == "affine":
        return 2.0 * ell + 2.0    # Kohel-affine, inversion not charged (batched)
    raise ValueError(model)


# --------------------------------------------------------------------------
# Class group / relation lattice machinery (maximal order via bnf)
# --------------------------------------------------------------------------
class ClassGroupCtx:
    def __init__(self, DK, g=1, verbose=False):
        assert is_fundamental(DK), "DK must be a fundamental discriminant"
        self.DK, self.g = DK, g
        self.D = DK * g * g
        t0 = time.time()
        self.bnf = pari.bnfinit(pari(f"x^2-({DK})"), 0)
        self.t_bnf = time.time() - t0
        self.cyc = [int(c) for c in pari("(b)->b.cyc")(self.bnf)]
        self.hK = int(pari("(b)->b.no")(self.bnf))
        # class number of the order of conductor g (exact sequence)
        hg = self.hK
        if g > 1:
            w = {-3: 3, -4: 2}.get(DK, 1)  # [O_K^* : O_g^*] = index of unit groups
            hg = self.hK * g
            for q, _ in factor_with_budget(g)[0]:
                hg = hg * (q - int(pari.kronecker(DK, q))) // q
            hg //= w
        self.h = hg
        self.bid = None
        if g > 1:
            self.bid = pari.idealstar(self.bnf, g, 1)
            self.bid_cyc = [int(c) for c in pari("(b)->b.cyc")(self.bid)]
            # subgroup of (O_K/g)^* generated by (Z/g)^* and O_K^*
            gens = []
            zn = pari.znstar(g, 1)
            for u in pari("(z)->z.gen")(zn):
                gens.append(self.ideallog(int(u)))
            gens.append(self.ideallog(-1))
            tu = pari("(b)->b.tu")(self.bnf)[1]
            gens.append(self.ideallog(tu))
            self.H_gens = gens

    def ideallog(self, x):
        return [int(c) for c in pari.ideallog(self.bnf, x, self.bid)]

    def frobenius_element(self, t, f):
        """pi = (t + f sqrt(DK))/2 as an nf element."""
        return pari(f"({t} + {f}*x)/2")

    def prime_data(self, ell, t=None, f=None):
        """One prime ideal above a split prime ell, with its class-group dlog,
        the order of its class, and the Frobenius eigenvalue pi mod l (the
        eigenvalue of Frobenius on the kernel E[l] of the isogeny it defines)."""
        pr = pari.idealprimedec(self.bnf, ell)
        if len(pr) != 2:
            return None
        P = pr[0]
        v = [int(c) for c in pari.bnfisprincipal(self.bnf, P, 0)]
        order = 1
        for vi, ci in zip(v, self.cyc):
            order = math.lcm(order, ci // math.gcd(ci, vi))
        d = dict(ell=ell, dlog=v, order_in_ClK=order, pr=P, pr_conj=pr[1])
        modpr = pari.nfmodprinit(self.bnf, P)
        gen2 = pari.lift(pari.nfbasistoalg(self.bnf, P[1]))
        d["generator"] = str(gen2)
        if ell % 2 == 1:
            # P = (l, x - c): c = sqrt(DK) mod l selecting this prime
            c0, c1 = Fraction(str(pari.polcoef(gen2, 0))), Fraction(str(pari.polcoef(gen2, 1)))
            inv2 = pow(2, -1, ell)
            n0 = int(c0.numerator) * pow(int(c0.denominator), -1, ell) % ell
            n1 = int(c1.numerator) * pow(int(c1.denominator), -1, ell) % ell
            d["sqrtD_mod_l"] = (-n0 * pow(n1, -1, ell)) % ell if n1 % ell else None
        if t is not None:
            pi = self.frobenius_element(t, f)
            lam = pari.nfmodpr(self.bnf, pi, modpr)
            d["frobenius_eigenvalue_mod_l"] = ffelt_to_int(lam)
        return d

    def relation_lattice(self, primes):
        """Basis (columns) of Lambda = {e : prod [l_i]^{e_i} = 1 in Cl(O)}."""
        k, m = len(primes), len(self.cyc)
        if m == 0:
            B = pari.matid(k)  # class number 1: every prime ideal is principal
        else:
            rows = []
            for i in range(m):
                rows.append([pd["dlog"][i] for pd in primes] + [self.cyc[i] if j == i else 0 for j in range(m)])
            A = mk_matrix(rows)
            Kmat = pari.matkerint(A)
            ncols = int(pari.matsize(Kmat)[1])
            proj = mk_matrix([[Kmat[i, j] for j in range(ncols)] for i in range(k)])
            B = pari.mathnf(proj)
        assert int(pari.matsize(B)[1]) == k, "relation lattice not full rank"
        if self.g > 1:
            B = self._restrict_to_order(B, primes)
        return B

    def generator(self, primes, e):
        """Generator alpha of prod l_i^{e_i} (principal by construction) as
        (a0, a1) with alpha = a0 + a1*sqrt(DK), a0, a1 in (1/2)Z."""
        prs = [pd["pr"] if int(x) >= 0 else pd["pr_conj"] for pd, x in zip(primes, e)]
        I = pari.idealfactorback(self.bnf, prs, [abs(int(x)) for x in e])
        res = pari.bnfisprincipal(self.bnf, I, 1)
        assert all(int(c) == 0 for c in res[0]), "ideal not principal"
        al = pari.lift(pari.nfbasistoalg(self.bnf, res[1]))
        a0 = Fraction(str(pari.polcoef(al, 0)))
        a1 = Fraction(str(pari.polcoef(al, 1)))
        return a0, a1, res[1]

    def _restrict_to_order(self, B, primes):
        """Sublattice of relations whose generator lies in O_g (up to O_K^*)."""
        k = int(pari.matsize(B)[1])
        logs = []
        for j in range(k):
            e = [int(B[i, j]) for i in range(k)]
            _, _, alg = self.generator(primes, e)
            logs.append(self.ideallog(alg))
        m = len(self.bid_cyc)
        cols = logs + self.H_gens + [[self.bid_cyc[i] if r == i else 0 for i in range(m)] for r in range(m)]
        A = mk_matrix([[c[i] for c in cols] for i in range(m)])
        Kmat = pari.matkerint(A)
        ncols = int(pari.matsize(Kmat)[1])
        proj = mk_matrix([[Kmat[i, j] for j in range(ncols)] for i in range(k)])
        Mg = pari.mathnf(proj)  # coefficients w.r.t. the columns of B
        Bg = pari.mathnf(B * Mg)
        return Bg


def weighted_shortest(B, weights, enum_cap=200000, dlogs=None, cyc=None, radius_cap=None, node_cap=20_000_000,
                      minimum_only=False):
    """Cheapest nonzero vectors of the relation lattice in the weighted 1-norm
    sum_i w_i |e_i|.

    Two stages: (1) LLL on the weighted Gram matrix B^T W B gives an upper
    bound `best` (the paper's Q_w approximation); (2) an exact depth-first
    enumeration of the weighted-1-norm ball {v in Z^k : sum w_i |v_i| <= R},
    R = min(best, radius_cap), testing membership v in Lambda through the
    class-group dlogs (sum v_i dlog_i = 0 mod cyc).  The 1-norm ball is far
    smaller than the Q_w ellipsoid that Lemma 1 of the paper would require, so
    the search is exhaustive within R whenever it finishes under `node_cap`.
    Returns (ranked list of (cost, e), exhaustive, R)."""
    k = int(pari.matsize(B)[1])
    wint = [max(1, int(round(2 * w))) for w in weights]  # integer weights (half-M units)
    W = pari.matdiagonal(wint)
    G = pari.mattranspose(B) * W * B
    T = pari.qflllgram(G)
    Bred = B * T

    def l1cost(e):
        return sum(w * abs(x) for w, x in zip(weights, e))

    cands = {}
    for j in range(k):
        e = tuple(int(Bred[i, j]) for i in range(k))
        cands[e] = l1cost(e)
    best = min(cands.values())
    lll_best = best
    R = best if radius_cap is None else radius_cap
    if minimum_only:
        # enumerating up to the LLL bound is enough to certify the minimum
        R = min(R, best)
    exhaustive = False
    if dlogs is not None:
        m = len(cyc)
        order = sorted(range(k), key=lambda i: -weights[i])  # heavy coordinates first
        nodes = 0
        hit_cap = False
        v = [0] * k
        acc = [0] * m

        def dfs(pos, budget):
            nonlocal nodes, hit_cap
            if hit_cap:
                return
            if pos == k:
                nodes += 1
                if any(v) and all(a % c == 0 for a, c in zip(acc, cyc)):
                    e = tuple(v)
                    cands[e] = l1cost(e)
                return
            i = order[pos]
            w = weights[i]
            lim = int(budget // w)
            dl = dlogs[i]
            for x in range(-lim, lim + 1):
                nodes += 1
                if nodes > node_cap:
                    hit_cap = True
                    return
                v[i] = x
                for t_ in range(m):
                    acc[t_] += x * dl[t_]
                dfs(pos + 1, budget - w * abs(x))
                for t_ in range(m):
                    acc[t_] -= x * dl[t_]
            v[i] = 0

        dfs(0, R)
        exhaustive = not hit_cap
    else:
        # fallback: Fincke-Pohst ball on the weighted quadratic form (not exhaustive in the 1-norm)
        res = pari.qfminim(pari.mattranspose(T) * G * T, int((2 * best) ** 2), enum_cap, 2)
        M = res[2]
        for j in range(int(pari.matsize(M)[1])):
            vv = [int(M[i, j]) for i in range(k)]
            e = tuple(int(sum(int(Bred[i, c]) * vv[c] for c in range(k))) for i in range(k))
            cands[e] = l1cost(e)
    # e and -e are the conjugate loops (alpha and alpha-bar); both are kept because
    # their eigenvalues lambda and N(alpha)/lambda are not interchangeable when
    # N(alpha) is large (e.g. pi acts as 1 but pi-bar acts as q).
    for e in list(cands):
        ne = tuple(-x for x in e)
        if ne not in cands:
            cands[ne] = cands[e]
    ranked = sorted(cands.items(), key=lambda kv: (kv[1], sum(abs(x) for x in kv[0]), [-x for x in kv[0]]))
    out = [(c, e) for e, c in ranked]
    return out, exhaustive, R, lll_best


# --------------------------------------------------------------------------
# Curve-level analysis
# --------------------------------------------------------------------------
def load_registry():
    with open(REGISTRY) as fh:
        d = json.load(fh)
    cs = d["curves"]
    return cs if isinstance(cs, list) else list(cs.values())


def resolve_curve(name, args):
    """Return a dict: name, family, q, char2, n (ext degree), trace, order, r, cofactor, a, j."""
    if name in BUILTIN:
        c = BUILTIN[name]
        if c["family"] == "prime":
            q = c["p"]
            order = c["h"] * c["n"]
            return dict(name=name, family="prime", q=q, char2=False, n=1, trace=q + 1 - order,
                        order=order, r=c["n"], cofactor=c["h"], a=c.get("a"), b=c.get("b"))
        q = 1 << c["m"]
        order = c["h"] * c["n"]
        return dict(name=name, family=c["family"], q=q, char2=True, n=c["m"], trace=q + 1 - order,
                    order=order, r=c["n"], cofactor=c["h"], a=c.get("a"), b=c.get("b"))
    if name == "explicit":
        if args.p:
            q, char2, n = int(args.p, 0), False, 1
        else:
            n = int(args.m)
            q, char2 = 1 << n, True
        if args.trace is not None:
            t = int(args.trace, 0)
            order = q + 1 - t
        else:
            order = int(args.order, 0)
            t = q + 1 - order
        r = int(args.r, 0) if args.r else largest_prime_factor(order)
        return dict(name="explicit", family="binary" if char2 else "prime", q=q, char2=char2, n=n,
                    trace=t, order=order, r=r, cofactor=order // r, a=None, b=None)
    for c in load_registry():
        names = set(c.get("aliases", [])) | set(c.get("standard_names", [])) | {c.get("slug")}
        if name in names or name.lower() in {s.lower() for s in names if s}:
            fam = c["family"]
            t, order = int(c["trace"]), int(c["order"])
            if fam in ("binary", "koblitz"):
                n = int(c["params"].get("n") or c["params"].get("m"))
                q, char2 = 1 << n, True
            else:
                q, char2, n = int(c["params"]["p"]), False, 1
            r = largest_prime_factor(order)
            a = c["params"].get("a")
            return dict(name=name, family=fam, q=q, char2=char2, n=n, trace=t, order=order, r=r,
                        cofactor=order // r, a=int(str(a), 0) if a is not None else None,
                        b=c["params"].get("b"), slug=c.get("slug"), end=c.get("end"))
    raise SystemExit(f"unknown curve {name!r}: use a registry alias, one of {sorted(BUILTIN)}, or 'explicit'")


def largest_prime_factor(N):
    fac, cof = factor_with_budget(N, seconds=120)
    cands = [q for q, _ in fac] + ([cof] if cof > 1 else [])
    return max(cands)


def glv_basis(r, lam):
    """Gauss-reduced basis of {(k1,k2): k1 + k2*lam = 0 mod r} and its quality."""
    M = pari(f"[{r}, {-lam % r}; 0, 1]")
    T = pari.qflll(M)
    R = M * T
    b1 = (int(R[0, 0]), int(R[1, 0]))
    b2 = (int(R[0, 1]), int(R[1, 1]))
    inf = max(abs(x) for x in b1 + b2)
    return dict(b1=b1, b2=b2, max_abs_entry_bits=inf.bit_length(),
                half_bits=(r.bit_length() + 1) // 2, ratio_to_sqrt_r=inf / math.sqrt(r))


def analyse_curve(cv, args):
    out = dict(curve={k: (str(v) if isinstance(v, int) and v.bit_length() > 60 else v)
                      for k, v in cv.items() if k not in ("b",)})
    q, t, r = cv["q"], cv["trace"], cv["r"]
    char2 = cv["char2"]
    delta = t * t - 4 * q
    if delta >= 0:
        out["error"] = "not ordinary / not an elliptic curve trace"
        return out
    t0 = time.time()
    fac, cof = factor_with_budget(delta, seconds=args.factor_seconds)
    DK, fpi = fundamental_and_conductor(delta, fac, cof)
    out["frobenius"] = dict(delta=str(delta), delta_bits=abs(delta).bit_length(),
                            factorization=[(str(p), e) for p, e in fac],
                            unfactored_cofactor=str(cof) if cof > 1 else None,
                            DK=str(DK), DK_bits=abs(DK).bit_length(), conductor_f=str(fpi),
                            conductor_factorization=[(str(p), e) for p, e in factor_with_budget(fpi)[0]] if fpi > 1 else [],
                            factor_seconds=round(time.time() - t0, 3))
    if cof > 1:
        out["frobenius"]["assumption"] = "unfactored cofactor treated as squarefree (a prime > 2^20 dividing the conductor would be missed)"
    # which orders to analyse
    conductors = [1]
    if args.conductor:
        conductors = [int(x) for x in args.conductor.split(",") if x]
        for g in conductors:
            if fpi % g:
                raise SystemExit(f"--conductor {g} does not divide f_pi = {fpi}")
    out["orders"] = []
    # eigenvalue bookkeeping: sqrt(DK) mod r with pi = 1 on E[r]
    s = None
    if r and fpi % r:
        s = ((2 - t) * pow(fpi, -1, r)) % r
        if (s * s - DK) % r:
            out["eigenvalue_error"] = "sqrt(DK) mod r inconsistent with pi = 1 (r does not divide #E?)"
            s = None
    out["sqrtDK_mod_r"] = str(s) if s is not None else None
    bits_r = r.bit_length()
    lprime = (bits_r + 1) // 2
    c_dbl = args.doubling_cost
    out["budget"] = dict(log2_r=bits_r, half_bits=lprime,
                         doublings_saved_cost_M=c_dbl * lprime,
                         cost_model=args.cost_model, doubling_cost_M=c_dbl)
    if abs(DK).bit_length() > args.max_disc_bits:
        out["skipped"] = (f"|DK| has {abs(DK).bit_length()} bits > --max-disc-bits {args.max_disc_bits}; "
                          "class group not computed")
        out["heuristic"] = heuristic_estimate(DK, args, char2, lprime * c_dbl)
        return out
    for g in conductors:
        out["orders"].append(analyse_order(DK, g, t, fpi, r, s, char2, args, lprime * c_dbl, cv))
    return out


def heuristic_estimate(DK, args, char2, saving):
    """Gaussian-heuristic estimate of the cheapest relation when the class group
    is too large to compute: det(Lambda) = h ~ sqrt(|D|) L(1,chi)/pi, shortest
    vector ~ sqrt(k/(2 pi e)) h^{1/k} in L2, cost >= c_min * L1 >= c_min * L2."""
    h_est = math.sqrt(abs(DK)) / math.pi  # L(1,chi) ~ 1
    rows = []
    for L in (50, 100, 200, 500):
        ps = [p for p in primes_upto(L) if int(pari.kronecker(DK, p)) == 1]
        k = len(ps)
        if k == 0:
            continue
        l2 = math.sqrt(k / (2 * math.pi * math.e)) * h_est ** (1.0 / k)
        l1 = l2 * math.sqrt(2 * k / math.pi)  # 1-norm of a random direction of that 2-norm
        cmin = min(step_cost(p, char2, args.cost_model) for p in ps)
        cmean = sum(step_cost(p, char2, args.cost_model) for p in ps) / k
        rows.append(dict(prime_bound=L, split_primes=k, log2_h_est=round(math.log2(h_est), 1),
                         shortest_L2_est=round(l2, 1), shortest_L1_est=round(l1, 1),
                         cost_lower_bound_M=round(cmin * l2),
                         cost_typical_M=round(cmean * l1), saving_M=saving,
                         ratio_lower_bound=round(cmin * l2 / saving, 2), ratio_typical=round(cmean * l1 / saving, 1)))
    return rows


def analyse_order(DK, g, t, fpi, r, s, char2, args, saving, cv):
    res = dict(conductor=g, D=str(DK * g * g))
    t0 = time.time()
    try:
        ctx = ClassGroupCtx(DK, g)
    except PariError as e:
        res["error"] = f"bnfinit failed: {e}"
        return res
    res["class_group"] = dict(h=str(ctx.h), hK=str(ctx.hK), cyc_K=[str(c) for c in ctx.cyc],
                              bnf_seconds=round(ctx.t_bnf, 2))
    primes = []
    for ell in primes_upto(args.prime_bound):
        if g % ell == 0 or int(pari.kronecker(DK, ell)) != 1:
            continue
        pd = ctx.prime_data(ell, t, fpi)
        if pd is None:
            continue
        pd["cost"] = step_cost(ell, char2, args.cost_model)
        primes.append(pd)
    if not primes:
        res["error"] = "no split primes below the bound"
        return res
    B = ctx.relation_lattice(primes)
    weights = [pd["cost"] for pd in primes]
    ranked, exhaustive, R, lll_best = weighted_shortest(B, weights, args.enum_cap,
                                              dlogs=[pd["dlog"] for pd in primes], cyc=ctx.cyc,
                                              radius_cap=args.radius_factor * saving)
    res["primes"] = [dict(ell=pd["ell"], cost_M=pd["cost"], dlog=pd["dlog"], order_in_ClK=pd["order_in_ClK"],
                          frobenius_eigenvalue_mod_l=pd.get("frobenius_eigenvalue_mod_l"),
                          sqrtDK_mod_l=pd.get("sqrtD_mod_l"), generator=pd["generator"]) for pd in primes]
    # single-prime loops (paper, Section 3): order m of [l] in Cl(O_g)
    single = []
    for idx, pd in enumerate(primes):
        m = pd["order_in_ClK"]
        if g > 1:
            # order of [l] in Cl(O_g): the least m with m*u_idx in the lattice spanned by B
            unit = pari.matrix(len(primes), 1, [1 if i == idx else 0 for i in range(len(primes))])
            coords = pari.matsolve(B, unit)
            m = 1
            for c in [Fraction(str(x)) for x in coords]:
                m = math.lcm(m, c.denominator)
        single.append(dict(ell=pd["ell"], loop_length=m, degree_bits=round(m * math.log2(pd["ell"]), 1),
                           cost_M=m * pd["cost"], ratio_to_saving=round(m * pd["cost"] / saving, 3)))
    single.sort(key=lambda d: d["cost_M"])
    res["single_prime_loops"] = single[: args.top]
    res["relation_search"] = dict(rank=len(primes), exhaustive_within_radius=exhaustive, radius_M=R,
                                  lll_upper_bound_M=lll_best, loops_within_radius=sum(1 for c, _ in ranked if c <= R),
                                  candidates=len(ranked))
    loops = []
    for cost, e in ranked[: args.top]:
        a0, a1, _ = ctx.generator(primes, e)
        norm = a0 * a0 - DK * a1 * a1
        trace = 2 * a0
        assert norm.denominator == 1 and trace.denominator == 1
        norm, trace = int(norm), int(trace)
        deg_check = 1
        for pd, x in zip(primes, e):
            deg_check *= pd["ell"] ** abs(x)
        assert norm == deg_check, (norm, deg_check)
        A2, B2 = int(2 * a0), int(2 * a1)
        # alpha = c + d*pi lies in Z[pi] iff f_pi | B and A = d*t (mod 2); then alpha acts
        # on all of E(F_q) as the integer c + d (pi = 1): a known scalar, i.e. the
        # "folklore" split of the paper (lambda = q^j type), not a genuine GLV endomorphism.
        in_zpi = (B2 % fpi == 0) and ((A2 - (B2 // fpi) * t) % 2 == 0)
        scalar_on_EFq = None
        if in_zpi:
            d_ = B2 // fpi
            scalar_on_EFq = (A2 - d_ * t) // 2 + d_
        # in characteristic 2 the loop {2: m*n} is pi^m (m > 0) or pi-bar^|m| (m < 0):
        # the Frobenius / Verschiebung chain, which acts on E(F_q) as 1 or q^|m|.
        nz = [(pd["ell"], int(x)) for pd, x in zip(primes, e) if x]
        frob_power = None
        if char2 and len(nz) == 1 and nz[0][0] == 2 and nz[0][1] % cv["n"] == 0:
            frob_power = nz[0][1] // cv["n"]
        entry = dict(cost_M=cost, ratio_to_saving=round(cost / saving, 3),
                     exponents={pd["ell"]: int(x) for pd, x in zip(primes, e) if x},
                     chain_length=sum(abs(x) for x in e), degree=str(norm), degree_bits=round(math.log2(norm), 1),
                     alpha=f"({A2} + {B2}*sqrt({DK}))/2", trace=str(trace),
                     in_order_Og=(2 * a1).numerator % (2 * g) == 0 if g > 1 else True,
                     in_Z_pi=in_zpi, Z_pi_is_OK=(fpi == 1), acts_as_integer=str(scalar_on_EFq) if in_zpi else None,
                     frobenius_power=frob_power)
        if s is not None:
            inv2 = pow(2, -1, r)
            lam = ((2 * a0).numerator * inv2 + (2 * a1).numerator * inv2 * s) % r
            lamc = ((2 * a0).numerator * inv2 - (2 * a1).numerator * inv2 * s) % r
            assert (lam * lam - trace * lam + norm) % r == 0
            entry["lambda"] = str(lam)
            entry["lambda_conj"] = str(lamc)
            entry["lambda_is_pm1"] = lam in (1, r - 1)
            if in_zpi:
                assert (scalar_on_EFq - lam) % r == 0, "Z[pi] scalar disagrees with eigenvalue"
            entry["glv"] = glv_basis(r, lam)
        loops.append(entry)
    res["loops"] = loops
    # cheapest loop that is not a Frobenius/Verschiebung power (binary curves):
    # from the exhaustive list under the radius, else from the LLL basis vectors
    genuine = None
    for cost, e in ranked:
        nz = [(pd["ell"], int(x)) for pd, x in zip(primes, e) if x]
        if char2 and len(nz) == 1 and nz[0][0] == 2 and nz[0][1] % cv["n"] == 0:
            continue
        genuine = dict(cost_M=cost, ratio_to_saving=round(cost / saving, 3),
                       exponents={pd["ell"]: int(x) for pd, x in zip(primes, e) if x},
                       chain_length=sum(abs(x) for x in e), within_radius=cost <= R)
        break
    res["cheapest_non_frobenius_loop"] = genuine
    res["seconds"] = round(time.time() - t0, 2)
    # Koblitz / Frobenius sanity: is the cheapest loop the Frobenius itself?
    if cv.get("family") == "koblitz" and loops and s is not None:
        a_k = cv.get("a") or 0
        mu = (-1) ** (1 - a_k)
        inv2 = pow(2, -1, r)
        lam_tau = ((mu + s) * inv2) % r
        o = 1
        x = lam_tau
        while x != 1 and o <= 4 * cv["n"]:
            x = x * lam_tau % r
            o += 1
        res["tau"] = dict(mu=mu, lambda_tau=str(lam_tau), order_mod_r=o, extension_degree=cv["n"],
                          note="tau^n = pi acts as 1 on E(F_q)[r]; <tau> x <-1> gives the 2n-element orbits used by rho")
    if cv.get("family") == "binary" and loops and any(lp.get("in_Z_pi") for lp in loops):
        res["note"] = ("binary ordinary curve: the only 2-isogenies are the Frobenius F (inseparable) and the "
                       "Verschiebung V, so the 2-loops are {2: +n} = pi (acts as 1) and {2: -n} = pi-bar = V^n "
                       "(acts as q = t - 1 on E(F_q)); both lie in Z[pi]. A GLV split with lambda = t - 1 is the "
                       "paper's folklore split: it trades the n/2 doublings of [t-1]P for n Verschiebung steps, "
                       "and the gain is model-dependent (x-only V step ~1M+2S, but y must be recovered).")
    if cv.get("family") == "koblitz" and loops:
        res["note"] = ("class number 1: every prime ideal is principal, every loop has length 1 and is an "
                       "explicit element a + b*tau of Z[tau]; the cheapest is tau itself (the Frobenius), "
                       "which the tau-adic (tau-NAF) expansion already exploits. Nothing beyond tau-NAF.")
    return res


# --------------------------------------------------------------------------
# Disc / scan modes
# --------------------------------------------------------------------------
def analyse_disc(DK, args, saving_bits=None):
    char2 = args.char2
    lprime = (saving_bits + 1) // 2 if saving_bits else args.bits // 2
    saving = args.doubling_cost * lprime
    ctx = ClassGroupCtx(DK, 1)
    primes = []
    for ell in primes_upto(args.prime_bound):
        if int(pari.kronecker(DK, ell)) != 1:
            continue
        pd = ctx.prime_data(ell)
        if pd is None:
            continue
        pd["cost"] = step_cost(ell, char2, args.cost_model)
        primes.append(pd)
    res = dict(DK=DK, h=ctx.hK, cyc=ctx.cyc, dmin=(-DK if DK % 4 == 0 else (1 - DK)) // 4,
               saving_M=saving, half_bits=lprime)
    if not primes:
        res["error"] = "no split primes below the bound"
        return res
    B = ctx.relation_lattice(primes)
    ranked, exhaustive, R, lll_best = weighted_shortest(B, [pd["cost"] for pd in primes], args.enum_cap,
                                              dlogs=[pd["dlog"] for pd in primes], cyc=ctx.cyc,
                                              radius_cap=args.radius_factor * saving, node_cap=2_000_000,
                                              minimum_only=True)
    best = []
    for cost, e in ranked[: args.top]:
        a0, a1, _ = ctx.generator(primes, e)
        norm = int(a0 * a0 - DK * a1 * a1)
        best.append(dict(cost_M=cost, ratio_to_saving=round(cost / saving, 3),
                         exponents={pd["ell"]: int(x) for pd, x in zip(primes, e) if x},
                         chain_length=sum(abs(x) for x in e), degree=norm,
                         alpha=f"({2*a0} + {2*a1}*sqrt({DK}))/2"))
    single = sorted((dict(ell=pd["ell"], loop_length=pd["order_in_ClK"], cost_M=pd["order_in_ClK"] * pd["cost"])
                     for pd in primes), key=lambda d: d["cost_M"])
    res.update(primes=[pd["ell"] for pd in primes], exhaustive_within_radius=exhaustive, radius_M=R, lll_upper_bound_M=lll_best, loops=best, single_prime_loops=single[:3])
    return res


def scan(args):
    rows = []
    lo, hi = args.scan_from, args.scan_to
    t0 = time.time()
    for D in range(-lo, -hi - 1, -1):
        if not is_fundamental(D):
            continue
        try:
            r = analyse_disc(D, args)
        except (PariError, AssertionError) as e:
            rows.append(dict(DK=D, error=str(e)[:80]))
            continue
        if "error" in r:
            continue
        b = r["loops"][0]
        rows.append(dict(DK=D, h=r["h"], cyc=r["cyc"], dmin=r["dmin"], best_cost_M=b["cost_M"],
                         ratio=b["ratio_to_saving"], exponents=b["exponents"], degree=b["degree"],
                         single_best=r["single_prime_loops"][0] if r["single_prime_loops"] else None))
    rows.sort(key=lambda d: d.get("best_cost_M", 1e18))
    return dict(mode="scan", from_=lo, to=hi, prime_bound=args.prime_bound, cost_model=args.cost_model,
                bits=args.bits, seconds=round(time.time() - t0, 1), rows=rows)


# --------------------------------------------------------------------------
# Self-test: explicit Vélu loops on a D = -71 crater curve over F_p
# --------------------------------------------------------------------------
def selftest(args):
    DK = args.selftest_disc
    log = []
    ctx = ClassGroupCtx(DK, 1)
    log.append(f"D = {DK}, h = {ctx.hK}, Cl = {ctx.cyc}")
    # a prime p = a'^2 - DK b'^2 with t = 2a' even (so 2 | f_pi and the three
    # 2-isogenies are rational: two horizontal, one descending), p ~ 2^30
    p = None
    for bp in range(2, 2000):
        for ap in range(10000, 60000):
            cand = ap * ap - DK * bp * bp
            if cand.bit_length() >= 31 and pari.isprime(cand):
                p, a_, b_ = cand, ap, bp
                break
        if p:
            break
    t = 2 * a_
    fpi = 2 * b_
    assert t * t - 4 * p == fpi * fpi * DK
    H = pari.polclass(DK)
    roots = [int(x) for x in pari.polrootsmod(H, p)]
    assert len(roots) == ctx.hK, "crater not fully F_p-rational for this p"
    log.append(f"p = {p} ({p.bit_length()} bits), t = +-{t}, f_pi = {fpi}; {len(roots)} crater j-invariants")

    def curve_from_j(j):
        c = pari.Mod(j, p) / (1728 - j)
        return pari.ellinit([3 * c, 2 * c])

    E0 = curve_from_j(roots[0])
    N0 = int(pari.ellcard(E0))
    tE = p + 1 - N0
    assert tE in (t, -t)
    t = tE  # the trace of the model we actually use
    r = largest_prime_factor(N0)
    s = ((2 - t) * pow(fpi, -1, r)) % r
    assert (s * s - DK) % r == 0
    log.append(f"E0: j = {roots[0]}, #E0 = {N0} = {N0 // r} * {r}, trace {t}")
    cof = N0 // r
    P = pari.ellmul(E0, pari("(E)->random(E)")(E0), cof)
    while int(pari.ellorder(E0, P)) != r:
        P = pari.ellmul(E0, pari("(E)->random(E)")(E0), cof)
    jset = set(roots)

    def curve_AB(E):
        v = pari.ellinit(E)
        return pari.lift(v[3]), pari.lift(v[4])

    def frob_eigenvalue_on_kernel(E, g, ell):
        """Eigenvalue (mod l) of Frobenius on the cyclic subgroup with kernel
        polynomial g, by lifting one kernel point to F_{p^{2d}} and comparing
        (x^p, y^p) with [lam](x, y) for lam = 0..l-1 (distinguishes lam from -lam)."""
        A, Bc = curve_AB(E)
        fa = pari.factormod(g, p)
        f1 = fa[0, 0]
        d = int(pari.poldegree(f1))
        T = pari.ffinit(p, 2 * d, "y")
        w = pari.ffgen(T, "w")
        xq = pari("(f,w)->polrootsmod(f, w)")(f1, w)[0]
        yq = pari.sqrt(xq ** 3 + A * xq + Bc)
        Eext = pari.ellinit([A * w ** 0, Bc * w ** 0])
        Q = [xq, yq]
        assert int(pari.ellisoncurve(Eext, Q)) == 1
        fr = [xq ** p, yq ** p]
        for lam in range(ell):
            if pari.ellmul(Eext, Q, lam) == fr:
                return lam
        raise AssertionError("no eigenvalue found")

    def eigenspace_kernels(E, ell, lam):
        """Kernel polynomials (degree (l-1)/2) of the cyclic l-subgroups on which
        Frobenius acts as [lam]: g = gcd(psi_l, x^p * D_lam - N_lam) where
        x([lam]Q) = N_lam/D_lam (x-only, so lam and -lam are merged; when the
        two eigenvalues are opposite, the factors are separated with a point lift)."""
        psi = pari.elldivpol(E, ell)
        xp = pari.lift(pari.Mod(pari("x"), psi) ** p)
        N, Dn = pari.ellxn(E, lam % ell)
        g = pari.gcd(psi, (xp * Dn - N) % psi)
        g = g / pari.pollead(g)
        if int(pari.poldegree(g)) == (ell - 1) // 2:
            return [g]
        # both eigenspaces share this x-polynomial (lam = -lam_conj): split by factoring
        from itertools import combinations
        fa = pari.factormod(g, p)
        irr = [fa[i, 0] for i in range(int(pari.matsize(fa)[0])) for _ in range(int(fa[i, 1]))]
        out = []
        for kk in range(1, len(irr) + 1):
            for comb in combinations(range(len(irr)), kk):
                h = pari("1")
                for i in comb:
                    h = h * irr[i]
                if int(pari.poldegree(h)) != (ell - 1) // 2:
                    continue
                try:
                    if frob_eigenvalue_on_kernel(E, h, ell) == lam % ell:
                        out.append(h)
                except (PariError, AssertionError):
                    continue
        return out

    def horizontal_neighbours(E, ell, pd=None):
        """[(kernel polynomial, codomain E', maps, j', lam)] for the horizontal
        l-isogenies of E; lam is the Frobenius eigenvalue on the kernel (None
        for l = 2 with rational 2-torsion, where it is 1 on every kernel)."""
        out = []
        if ell == 2:
            A, Bc = curve_AB(E)
            rts = pari.polrootsmod(pari("x^3 + %s*x + %s" % (A, Bc)), p)
            cands = [(pari("x - %s" % pari.lift(x)), None) for x in rts]
        else:
            lam_l = pd["frobenius_eigenvalue_mod_l"]
            lam_c = (t - lam_l) % ell
            cands = [(g, lam_l) for g in eigenspace_kernels(E, ell, lam_l)]
            if lam_c != lam_l:
                cands += [(g, lam_c) for g in eigenspace_kernels(E, ell, lam_c)]
        for g, lam in cands:
            try:
                iso = pari.ellisogeny(E, g)
            except PariError:
                continue
            E2 = pari.ellinit(iso[0])
            j2 = int(pari.lift(pari("(E)->E.j")(E2)))
            if j2 in jset:
                out.append((g, E2, iso[1], j2, lam))
        return out

    primes = []
    for ell in (2, 3, 5, 7, 11, 13):
        if int(pari.kronecker(DK, ell)) != 1:
            continue
        pd = ctx.prime_data(ell, t, fpi)
        pd["cost"] = step_cost(ell, False, "projective")
        primes.append(pd)
    B = ctx.relation_lattice(primes)
    ranked, _, _, _ = weighted_shortest(B, [pd["cost"] for pd in primes], dlogs=[pd["dlog"] for pd in primes], cyc=ctx.cyc, radius_cap=400)
    cyc = ctx.cyc

    def addlab(a, b, sign=1):
        return tuple((x + sign * y) % c for x, y, c in zip(a, b, cyc))

    # Label the crater as the Cayley graph of Cl(O): E0 <-> identity, and the
    # l-neighbour of a vertex with label c in the direction of the prime ideal
    # l (the one whose Frobenius eigenvalue pi mod l was recorded) has label
    # c + dlog(l).  Primes with l | f_pi have both eigenvalues equal (E[l] is
    # rational) and cannot be oriented by Frobenius; they are oriented by the
    # labels once the unambiguous primes have labelled the crater.
    labels = {roots[0]: tuple([0] * len(cyc))}
    models = {roots[0]: E0}
    unamb = [pd for pd in primes if fpi % pd["ell"] != 0]
    amb = [pd for pd in primes if fpi % pd["ell"] == 0]
    frontier = [roots[0]]
    while frontier and len(labels) < len(roots):
        nxt = []
        for j in frontier:
            E = models[j]
            for pd in unamb:
                ell = pd["ell"]
                lam_l = pd["frobenius_eigenvalue_mod_l"]
                lam_c = (t - lam_l) % ell
                for g, E2, maps, j2, lam in horizontal_neighbours(E, ell, pd):
                    if lam == lam_l:
                        lab = addlab(labels[j], pd["dlog"], +1)
                    elif lam == lam_c:
                        lab = addlab(labels[j], pd["dlog"], -1)
                    else:
                        raise AssertionError("Frobenius eigenvalue on kernel is neither root")
                    if j2 in labels:
                        assert labels[j2] == lab, "inconsistent crater labelling"
                    else:
                        labels[j2] = lab
                        models[j2] = E2
                        nxt.append(j2)
        frontier = nxt
    assert len(labels) == len(roots), f"unambiguous primes {[pd['ell'] for pd in unamb]} do not generate Cl(O); labelled {len(labels)}/{len(roots)}"
    assert len(set(labels.values())) == len(roots), "labels not injective (class group action not free?)"
    log.append(f"crater labelled by Cl(O) via primes {[pd['ell'] for pd in unamb]}; primes needing label orientation (l | f_pi): {[pd['ell'] for pd in amb]}")

    def step(E, j, pd, sign):
        ell, dlog = pd["ell"], pd["dlog"]
        want = addlab(labels[j], dlog, sign)
        for g, E2, maps, j2, lam in horizontal_neighbours(E, ell, pd):
            if labels[j2] == want:
                if lam is not None:
                    # cross-check: the label direction agrees with the Frobenius eigenvalue
                    assert lam == (pd["frobenius_eigenvalue_mod_l"] if sign > 0 else (t - pd["frobenius_eigenvalue_mod_l"]) % ell), "label/eigenvalue mismatch"
                return E2, maps, j2
        raise AssertionError(f"no horizontal {ell}-isogeny from j={j} to label {want}")

    # test: the single-prime 2-loop and the cheapest multi-prime relations
    tests = []
    idx2 = next((i for i, pd in enumerate(primes) if pd["ell"] == 2), None)
    if idx2 is not None:
        e = [0] * len(primes)
        e[idx2] = primes[idx2]["order_in_ClK"]
        tests.append(tuple(e))
    for c, e in ranked[:3]:
        if tuple(e) not in tests:
            tests.append(tuple(e))
    results = []
    for e in tests:
        a0, a1, _ = ctx.generator(primes, e)
        norm = int(a0 * a0 - DK * a1 * a1)
        tr = int(2 * a0)
        inv2 = pow(2, -1, r)
        lam = (int(2 * a0) * inv2 + int(2 * a1) * inv2 * s) % r
        lamc = (int(2 * a0) * inv2 - int(2 * a1) * inv2 * s) % r
        E, Q, j = E0, P, roots[0]
        chain = []
        for pd, x in zip(primes, e):
            if x == 0:
                continue
            for _ in range(abs(x)):
                E, maps, j = step(E, j, pd, 1 if x > 0 else -1)
                Q = pari.ellisogenyapply(maps, Q)
                chain.append((pd["ell"], j))
        jN = int(pari.lift(pari("(E)->E.j")(E)))
        assert jN == roots[0], f"loop did not close: j = {jN}"
        # isomorphism back to E0: A_N = u^4 A_0, B_N = u^6 B_0, (x,y) -> (x/u^2, y/u^3)
        A0, B0 = pari.ellinit(E0)[3], pari.ellinit(E0)[4]
        AN, BN = pari.ellinit(E)[3], pari.ellinit(E)[4]
        u2 = (BN * A0) / (AN * B0)
        u = pari.sqrt(u2)
        assert pari.lift(u) == pari.lift(u) and (u2 ** ((p - 1) // 2) == 1), "u not rational: composite not F_p-rational"
        R = [Q[0] / u2, Q[1] / (u2 * u)]
        assert int(pari.ellisoncurve(E0, R)) == 1
        match = None
        for name, l in (("lambda", lam), ("-lambda", r - lam), ("lambda_conj", lamc), ("-lambda_conj", r - lamc)):
            if pari.ellmul(E0, P, l) == R:
                match = name
                break
        results.append(dict(exponents={pd["ell"]: x for pd, x in zip(primes, e) if x}, chain=chain,
                            alpha=f"({2*a0} + {2*a1}*sqrt({DK}))/2", degree=norm, trace=tr,
                            lambda_=lam, lambda_conj=lamc, action=match,
                            cost_M=sum(pd["cost"] * abs(x) for pd, x in zip(primes, e))))
        log.append(f"loop {results[-1]['exponents']}: degree {norm}, alpha = {results[-1]['alpha']}, acts as {match}")
    # the generator is defined up to the unit -1 (the sign of the final isomorphism u),
    # so the composite acts as +lambda or -lambda; a conjugate would mean a wrong direction
    ok = all(rr["action"] in ("lambda", "-lambda") for rr in results)
    return dict(mode="selftest", DK=DK, p=p, trace=t, f_pi=fpi, r=r, h=ctx.hK,
                primes=[dict(ell=pd["ell"], frob_eig=pd["frobenius_eigenvalue_mod_l"], order=pd["order_in_ClK"]) for pd in primes],
                loops=results, all_loops_act_as_lambda=ok, log=log)


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------
def md_curve_report(out):
    cv = out["curve"]
    lines = [f"### {cv['name']} ({cv['family']}, q bits {int(cv['q']).bit_length() if isinstance(cv['q'], (int, str)) else '?'}, log2 r = {out.get('budget', {}).get('log2_r')})"]
    fr = out.get("frobenius")
    if fr:
        lines.append(f"- Delta = t^2 - 4q: {fr['delta_bits']} bits, D_K = {fr['DK']} ({fr['DK_bits']} bits), conductor f_pi = {fr['conductor_f']}"
                     + (f" (f_pi = {fr['conductor_factorization']})" if fr['conductor_factorization'] else ""))
        if fr.get("assumption"):
            lines.append(f"- assumption: {fr['assumption']}")
    if "skipped" in out:
        lines.append(f"- {out['skipped']}")
        lines.append("")
        lines.append("| prime bound | split primes | log2 h (est.) | shortest L2 (est.) | L1 (est.) | cost lower bound M | typical M | saving M | ratio (lower bound) | ratio (typical) |")
        lines.append("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
        for h in out["heuristic"]:
            lines.append(f"| {h['prime_bound']} | {h['split_primes']} | {h['log2_h_est']} | {h['shortest_L2_est']} | {h['shortest_L1_est']} | {h['cost_lower_bound_M']} | {h['cost_typical_M']} | {h['saving_M']} | {h['ratio_lower_bound']} | {h['ratio_typical']} |")
        return "\n".join(lines)
    b = out["budget"]
    lines.append(f"- GLV budget: the split saves ~{b['half_bits']} doublings = {b['doublings_saved_cost_M']} M ({b['cost_model']} step model, doubling {b['doubling_cost_M']} M)")
    for od in out["orders"]:
        if "error" in od:
            lines.append(f"- order of conductor {od['conductor']}: {od['error']}")
            continue
        cg = od["class_group"]
        lines.append(f"- order O_{od['conductor']} (D = {od['D']}): h = {cg['h']}, Cl(O_K) = {cg['cyc_K']} (bnfinit {cg['bnf_seconds']} s); "
                     f"{od['relation_search']['rank']} split primes <= bound; loops within {od['relation_search']['radius_M']:.0f} M: {od['relation_search']['loops_within_radius']} (exhaustive: {od['relation_search']['exhaustive_within_radius']}); LLL upper bound {od['relation_search']['lll_upper_bound_M']:.0f} M")
        if od.get("note"):
            lines.append(f"- {od['note']}")
        lines.append("")
        lines.append("| rank | cost M | cost / saving | loop (prime: exponent) | chain length | degree bits | alpha | Frobenius power / in Z[pi] | lambda is ±1 | GLV max |k_i| bits |")
        lines.append("|--:|--:|--:|:--|--:|--:|:--|:--|:--|--:|")
        for i, lp in enumerate(od["loops"]):
            g = lp.get("glv", {})
            if lp.get("frobenius_power") is not None:
                m_ = lp["frobenius_power"]
                zp = f"pi^{m_}" if m_ > 0 else f"pi-bar^{-m_}"
                zp += f" (acts as {lp['acts_as_integer'][:12]}{'…' if len(lp['acts_as_integer']) > 12 else ''})"
            elif lp.get("Z_pi_is_OK"):
                zp = "Z[pi] = O_K"
            else:
                zp = f"yes ({lp['acts_as_integer'][:12]}{'…' if len(lp['acts_as_integer']) > 12 else ''})" if lp.get("in_Z_pi") else "no"
            lines.append(f"| {i+1} | {lp['cost_M']:.1f} | {lp['ratio_to_saving']} | {lp['exponents']} | {lp['chain_length']} | {lp['degree_bits']} | {lp['alpha']} | {zp} | {lp.get('lambda_is_pm1')} | {g.get('max_abs_entry_bits', '')} |")
        gl = od.get("cheapest_non_frobenius_loop")
        if gl:
            lines.append("")
            lines.append(f"cheapest loop that is not a Frobenius/Verschiebung power: {gl['cost_M']:.1f} M ({gl['ratio_to_saving']}x saving), {gl['exponents']}, chain length {gl['chain_length']}"
                         + ("" if gl["within_radius"] else " — above the enumeration radius, from the LLL basis (upper bound, not necessarily the cheapest)"))
        sp = od["single_prime_loops"][:3]
        if sp:
            lines.append("")
            lines.append("single-prime loops (paper Section 3): " + "; ".join(f"l = {s['ell']}: length {s['loop_length']}, {s['cost_M']:.0f} M ({s['ratio_to_saving']}x saving)" for s in sp))
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["curve", "disc", "scan", "selftest"])
    ap.add_argument("names", nargs="*", help="curve aliases / registry slugs / built-ins (curve mode), or discriminants (disc mode)")
    ap.add_argument("--p", help="prime field (explicit curve)")
    ap.add_argument("--m", help="binary extension degree (explicit curve)")
    ap.add_argument("--trace", help="Frobenius trace (explicit curve)")
    ap.add_argument("--order", help="#E(F_q) (explicit curve; alternative to --trace)")
    ap.add_argument("--r", help="prime subgroup order (explicit curve; default: largest prime factor of #E)")
    ap.add_argument("--conductor", help="comma-separated conductors g | f_pi to analyse O_g (default: maximal order only)")
    ap.add_argument("--prime-bound", type=int, default=100, help="use split primes l <= this bound (default 100)")
    ap.add_argument("--cost-model", choices=["projective", "affine"], default="projective")
    ap.add_argument("--doubling-cost", type=float, default=8.0, help="field multiplications per doubling (paper: 8)")
    ap.add_argument("--bits", type=int, default=256, help="log2 r for disc/scan modes (saving = doubling_cost * bits/2)")
    ap.add_argument("--char2", action="store_true", help="disc/scan: binary-field step costs (Frobenius for l = 2)")
    ap.add_argument("--max-disc-bits", type=int, default=140, help="skip the class-group computation above this |D_K| (heuristic estimate instead)")
    ap.add_argument("--factor-seconds", type=int, default=60)
    ap.add_argument("--enum-cap", type=int, default=200000)
    ap.add_argument("--radius-factor", type=float, default=1.0, help="exhaustively enumerate loops costing up to this multiple of the GLV saving (default 1: every loop that beats the doublings it saves)")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--scan-from", type=int, default=3)
    ap.add_argument("--scan-to", type=int, default=2000)
    ap.add_argument("--selftest-disc", type=int, default=-71)
    ap.add_argument("--json", help="write the full result to this file")
    args = ap.parse_args()

    t0 = time.time()
    if args.mode == "curve":
        names = args.names or ["explicit"]
        results = []
        for nm in names:
            cv = resolve_curve(nm, args)
            sys.stderr.write(f"[{time.strftime('%H:%M:%S')}] {nm}: q {cv['q'].bit_length()} bits, |t| {abs(cv['trace']).bit_length()} bits\n")
            try:
                res = analyse_curve(cv, args)
            except (PariError, AssertionError) as e:
                res = dict(curve=dict(name=nm), error=f"{type(e).__name__}: {e}")
            results.append(res)
            print(md_curve_report(res) if "error" not in res or "curve" in res and "frobenius" in res else f"### {nm}\n- error: {res['error']}")
            print()
        payload = dict(mode="curve", args=vars(args), seconds=round(time.time() - t0, 1), results=results)
    elif args.mode == "disc":
        results = [analyse_disc(int(d), args) for d in args.names]
        for r in results:
            print(json.dumps(r, indent=1, default=str))
        payload = dict(mode="disc", args=vars(args), results=results)
    elif args.mode == "scan":
        payload = scan(args)
        print(f"| D_K | h | Cl | d_min | best loop cost M | cost/saving | loop | degree |")
        print("|--:|--:|:--|--:|--:|--:|:--|--:|")
        for row in payload["rows"][: args.top * 10]:
            if "error" in row:
                continue
            print(f"| {row['DK']} | {row['h']} | {row['cyc']} | {row['dmin']} | {row['best_cost_M']:.1f} | {row['ratio']} | {row['exponents']} | {row['degree']} |")
    else:
        payload = selftest(args)
        for l in payload["log"]:
            print(l)
        print("ALL LOOPS ACT AS [+-lambda] (never the conjugate):", payload["all_loops_act_as_lambda"])
        if not payload["all_loops_act_as_lambda"]:
            sys.exit(1)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(payload, fh, indent=1, default=str)
        sys.stderr.write(f"wrote {args.json}\n")


if __name__ == "__main__":
    main()
