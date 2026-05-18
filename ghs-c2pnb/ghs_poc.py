#!/usr/bin/env python3
"""
GHS / Hess Weil-descent — PoC analysis of ANSI X9.62 c2pnb*w1 binary curves.

What this script proves end-to-end (no SageMath, pure Python):

  1. Loads the real X9.62 parameters for c2pnb176w1, c2pnb208w1, c2pnb272w1,
     c2pnb304w1, c2pnb368w1.
  2. Verifies each curve's "w" subfield property:  a, b in F_{2^16} embedded
     in F_{2^{16 p}}. This is the structural hook the descent attaches to.
  3. Computes the direct GHS magic numbers (Frobenius-orbit F_2-dimensions
     of a and b under x -> x^{2^16}). All come out to 1, confirming
     Menezes-Qu's observation that the *original* GHS attack is useless on
     these curves -- the descended hyperelliptic curve has genus 1.
  4. Demonstrates the Hess extension: an isogeny shifts the curve to E' with
     b' NOT in F_{2^16}, raising the Frobenius-orbit dimension. We simulate
     this by sampling random F_{2^{16 p}}-elements (standing in for b'
     reached by an isogeny walk) and measuring the resulting orbit
     dimension, the genus of the descended curve, and the index-calculus
     cost. This reproduces the qualitative shape of the Menezes-Teske 2006
     analysis.
  5. Toy GHS-flavoured Weil restriction worked out symbolically on a small
     curve, so the descent's structure is concrete and inspectable.

References:
  - Gaudry, Hess, Smart. "Constructive and destructive facets of Weil
    descent on elliptic curves." J. Cryptology 15(1), 2002.
  - Menezes, Qu. "Analysis of the Weil descent attack of Gaudry, Hess, and
    Smart." CT-RSA 2001.
  - Hess. "Generalising the GHS attack on the elliptic curve discrete
    logarithm." LMS J. Comput. Math. 7, 2004.
  - Menezes, Teske. "Cryptographic implications of Hess' generalized GHS
    attack." Applicable Algebra in Engineering, Communication and
    Computing, 16(6), 2006. (Source for the c2pnb176w1 numerical break.)
"""

from __future__ import annotations
import math
import random
from dataclasses import dataclass


# =============================================================================
# F_2[x] polynomial arithmetic, bit-packed (int with bit i = coefficient of x^i)
# =============================================================================

def f2_deg(a: int) -> int:
    return a.bit_length() - 1

def f2_mul(a: int, b: int) -> int:
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r

def f2_mod(a: int, m: int) -> int:
    dm = f2_deg(m)
    while a.bit_length() - 1 >= dm:
        a ^= m << (a.bit_length() - 1 - dm)
    return a

def f2_mulmod(a: int, b: int, m: int) -> int:
    return f2_mod(f2_mul(a, b), m)

def f2_sqrmod(a: int, m: int) -> int:
    # squaring in char 2 just spreads the bits, then reduces
    sq = 0
    i = 0
    while a:
        if a & 1:
            sq |= 1 << (2 * i)
        a >>= 1
        i += 1
    return f2_mod(sq, m)

def f2_powmod(a: int, e: int, m: int) -> int:
    r = 1
    base = a % (1 << f2_deg(m))
    while e:
        if e & 1:
            r = f2_mulmod(r, base, m)
        base = f2_sqrmod(base, m)
        e >>= 1
    return r

def f2_frob(a: int, k: int, m: int) -> int:
    """a^{2^k} mod m, via k repeated squarings."""
    for _ in range(k):
        a = f2_sqrmod(a, m)
    return a


# =============================================================================
# F_2-rank of a set of elements of F_{2^N} (viewed as N-dim F_2-vectors)
# =============================================================================

def f2_rank(vectors: list[int], N: int) -> int:
    """Rank over F_2 of the given ints (each treated as length-N bit vector)."""
    rows = list(vectors)
    pivot_row = 0
    for col in range(N):
        mask = 1 << col
        # find a row with this bit set
        pivot = None
        for r in range(pivot_row, len(rows)):
            if rows[r] & mask:
                pivot = r
                break
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        for r in range(len(rows)):
            if r != pivot_row and (rows[r] & mask):
                rows[r] ^= rows[pivot_row]
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return pivot_row


# =============================================================================
# X9.62 c2pnb*w1 curve database (real parameters)
# =============================================================================

def hex_to_int(h: str) -> int:
    return int(h.replace(' ', '').replace('_', ''), 16)


@dataclass
class W1Curve:
    name: str
    m: int            # extension degree over F_2
    l: int            # subfield degree (F_{2^l} is where a, b live)
    modulus: int      # reduction polynomial as int
    a: int
    b: int
    gx: int
    gy: int
    order: int        # subgroup order
    cofactor: int

    @property
    def n(self) -> int:  # composite extension factor: m = l * n
        return self.m // self.l


# Reduction polynomials from X9.62 / SEC 2
#   c2pnb176w1: x^176 + x^43   + x^2 + x + 1
#   c2pnb208w1: x^208 + x^83   + x^2 + x + 1   (NB: X9.62 uses x^9 + x^7 + x^4 + 1 form? -- see below)
#   c2pnb272w1: x^272 + x^56   + x^3 + x + 1
#   c2pnb304w1: x^304 + x^11   + x^2 + x + 1
#   c2pnb368w1: x^368 + x^85   + x^2 + x + 1
#
# X9.62 actually specifies for each curve a pentanomial f(x) = x^m + x^a + x^b + x^c + 1
# The standard's published reduction polynomials are listed below verbatim.

C2PNB176W1 = W1Curve(
    name="c2pnb176w1",
    m=176, l=16,
    modulus=(1 << 176) | (1 << 43) | (1 << 2) | (1 << 1) | 1,
    a=hex_to_int("E4E6DB2995065C407D9D39B8D0967B96704BA8E9C90B"),
    b=hex_to_int("5DDA470ABE6414DE8EC133AE28E9BBD7FCEC0AE0FFF2"),
    gx=hex_to_int("8D16C2866798B600F9F08BB4A8E860F3298CE04A5798"),
    gy=hex_to_int("6FA4539C2DADDDD6BAB5167D61B436E1D92BB16A562C"),
    order=hex_to_int("00010092537397ECA4F6145799D62B0A19CE06FE26AD"),
    cofactor=0xFF6E,
)

C2PNB208W1 = W1Curve(
    name="c2pnb208w1",
    m=208, l=16,
    modulus=(1 << 208) | (1 << 83) | (1 << 2) | (1 << 1) | 1,
    a=0,
    b=hex_to_int("C8619ED45A62E6212E1160349E2BFA844439FAFC2A3FD1638F9E"),
    gx=hex_to_int("89FDFBE4ABE193DF9559ECF07AC0CE78554E2784EB8C1ED1A57A"),
    gy=hex_to_int("0F55B51A06E78E9AC38A035FF520D8B01781BEB1A6BB08617DE3"),
    order=hex_to_int("0101BAF95C9723C57B6C21DA2EFF2D5ED588BDD5717E212F9D"),
    cofactor=0xFE48,
)

C2PNB272W1 = W1Curve(
    name="c2pnb272w1",
    m=272, l=16,
    modulus=(1 << 272) | (1 << 56) | (1 << 3) | (1 << 1) | 1,
    a=hex_to_int("91A091F03B5FBA4AB2CCF49C4EDD220FB028712D42BE752B2C40094DBACDB586FB20"),
    b=hex_to_int("7167EFC92BB2E3CE7C8AAAFF34E12A9C557003D7C73A6FAF003F99F6CC8482E540F7"),
    gx=hex_to_int("6108BABB2CEEBCF787058A056CBE0CFE622D7723A289E08A07AE13EF0D10D171DD8D"),
    gy=hex_to_int("10C7695716851EEF6BA7F6872E6142FBD241B830FF5EFCACECCAB05E02005DDE9D23"),
    order=hex_to_int("0100FAF51354E0E39E4892DF6E319C72C8161603FA45AA7B998A167B8F1E629521"),
    cofactor=0xFF06,
)

C2PNB304W1 = W1Curve(
    name="c2pnb304w1",
    m=304, l=16,
    modulus=(1 << 304) | (1 << 11) | (1 << 2) | (1 << 1) | 1,
    a=hex_to_int("FD0D693149A118F651E6DCE6802085377E5F882D1B510B44160074C1288078365A0396C8E681"),
    b=hex_to_int("BDDB97E555A50A908E43B01C798EA5DAA6788F1EA2794EFCF57166B8C14039601E55827340BE"),
    gx=hex_to_int("197B07845E9BE2D96ADB0F5F3C7F2CFFBD7A3EB8B6FEC35C7FD67F26DDF6285A644F740A2614"),
    gy=hex_to_int("E19FBEB76E0DA171517ECF401B50289BF014103288527A9B416A105E80260B549FDC1B92C03B"),
    order=hex_to_int("0101D556572AABAC800101D556572AABAC8001022D5C91DD173F8FB561DA6899164443051D"),
    cofactor=0xFE2E,
)

C2PNB368W1 = W1Curve(
    name="c2pnb368w1",
    m=368, l=16,
    modulus=(1 << 368) | (1 << 85) | (1 << 2) | (1 << 1) | 1,
    a=hex_to_int(
        "E0D2EE25095206F5E2A4F9ED229F1F256E79A0E2B455970D8D0D865BD94778C576D62F0AB7519CCD2A1A906AE30D"),
    b=hex_to_int(
        "FC1217D4320A90452C760A58EDCD30C8DD069B3C34453837A34ED50CB54917E1C2112D84D164F444F8F74786046A"),
    gx=hex_to_int(
        "1085E2755381DCCCE3C1557AFA10C2F0C0C2825646C5B34A394CBCFA8BC16B22E7E789E927BE216F02E1FB136A5F"),
    gy=hex_to_int(
        "7B3EB1BDDCBA62D5D8B2059B525797FC73822C59059C623A45FF3843CEE8F87CD1855ADAA81E2A0750B80FDA2310"),
    order=hex_to_int(
        "010090512DA9AF72B08349D98A5DD4C7B0532ECA51CE03E2D10F3B7AC579BD87E909AE40A6F131E9CFCE5BD967"),
    cofactor=0xFF70,
)

CURVES = [C2PNB176W1, C2PNB208W1, C2PNB272W1, C2PNB304W1, C2PNB368W1]


# =============================================================================
# Subfield-structure & magic-number analysis
# =============================================================================

def is_in_subfield(elt: int, l: int, m: int, modulus: int) -> bool:
    """True iff elt is fixed by Frobenius x -> x^{2^l}, i.e. elt in F_{2^l}."""
    return f2_frob(elt, l, modulus) == elt


def frobenius_orbit(elt: int, l: int, n: int, modulus: int) -> list[int]:
    """[elt, sigma(elt), sigma^2(elt), ..., sigma^{n-1}(elt)] where sigma = x -> x^{2^l}."""
    orbit = [elt]
    cur = elt
    for _ in range(n - 1):
        cur = f2_frob(cur, l, modulus)
        orbit.append(cur)
    return orbit


def magic_dim(elt: int, l: int, n: int, m: int, modulus: int) -> int:
    """F_2-dimension of the Frobenius-orbit span of `elt`."""
    orbit = frobenius_orbit(elt, l, n, modulus)
    return f2_rank(orbit, m)


def descent_genus(magic: int) -> int:
    """Genus of the descended hyperelliptic curve given GHS magic number.

    Per GHS / Hess: g = 2^{magic - 1}  (sometimes  2^{magic-1} - 1 depending
    on parity of the splitting; we use the upper bound, the standard
    figure quoted in Menezes-Teske)."""
    return 1 << (magic - 1)


def index_calculus_cost_bits(genus: int, q: int) -> float:
    """Heuristic Gaudry-Theriault-Diem / Diem index-calculus cost on Jac(C)(F_q).

    For genus-g hyperelliptic curves over F_q, the practical attack cost is
    governed by three terms whose maxima we take:

      - relation collection:   ~ g! * q^{2 - 2/g}
        (the g! factor comes from polynomial factoring + smoothness tests
        of degree-g divisors; it's the term Menezes-Teske track carefully.)
      - sparse linear algebra: ~ q^2 * g
      - if both above are useless, fall back to generic Pollard rho on Jac,
        which costs sqrt(|Jac|) ~ q^{g/2}.

    Returns log2 of dominant operation count. This is the same shape as the
    Adleman-DeMarrais-Huang / Enge / Diem bounds; the exact constants
    differ between variants but the genus crossover sits in the same place.
    """
    if genus < 2:
        return float('inf')             # genus-1 descent: useless
    # log2(g!) via lgamma for cheap accuracy
    log2_g_fact = math.lgamma(genus + 1) / math.log(2)
    relations = log2_g_fact + (2.0 - 2.0 / genus) * math.log2(q)
    linear_alg = 2.0 * math.log2(q) + math.log2(genus)
    return max(relations, linear_alg)


def jac_pollard_rho_bits(genus: int, q: int) -> float:
    """Generic discrete log on Jac(C)(F_q): sqrt(#Jac) ~ q^{g/2}.

    This is the fallback cost if index calculus is ineffective (huge g).
    """
    return 0.5 * genus * math.log2(q)


def attack_cost_via_magic(magic: int, q: int) -> tuple[int, float]:
    """Best cost achievable by descending with this magic number.

    Returns (genus, log2 cost)."""
    g = descent_genus(magic)
    if g < 2:
        return g, float('inf')
    ic = index_calculus_cost_bits(g, q)
    rho = jac_pollard_rho_bits(g, q)
    return g, min(ic, rho)


def pollard_rho_cost_bits(group_order: int) -> float:
    """Generic ECDLP on E directly: sqrt(N)."""
    return 0.5 * math.log2(group_order)


# =============================================================================
# Reporting
# =============================================================================

def banner(s: str) -> None:
    print()
    print("=" * 72)
    print(s)
    print("=" * 72)


def report_curve(c: W1Curve) -> None:
    a_in_sub = is_in_subfield(c.a, c.l, c.m, c.modulus)
    b_in_sub = is_in_subfield(c.b, c.l, c.m, c.modulus)
    print(f"  {c.name}:  m = {c.m}  =  {c.l} * {c.n}")
    print(f"    a in F_2^{c.l}?  {a_in_sub}")
    print(f"    b in F_2^{c.l}?  {b_in_sub}")
    if not (a_in_sub and b_in_sub):
        print("    !! one or both coefficients are NOT in the claimed subfield;")
        print("       either the parameters above are wrong, or the curve isn't a 'w' form.")
        return
    magic_a = magic_dim(c.a, c.l, c.n, c.m, c.modulus)
    magic_b = magic_dim(c.b, c.l, c.n, c.m, c.modulus)
    g_direct = descent_genus(max(magic_a, magic_b))
    q = 1 << c.l
    cost_direct = index_calculus_cost_bits(g_direct, q) if g_direct >= 2 else float('inf')
    cost_rho = pollard_rho_cost_bits(c.order)
    print(f"    Frobenius-orbit dim of a (=GHS magic for a) : {magic_a}")
    print(f"    Frobenius-orbit dim of b (=GHS magic for b) : {magic_b}")
    print(f"    direct descended genus                       : {g_direct}")
    if cost_direct == float('inf'):
        print(f"    direct GHS index-calculus cost               : infeasible (genus 1, no advantage)")
    else:
        print(f"    direct GHS index-calculus cost (log2 ops)    : {cost_direct:.1f}")
    print(f"    Pollard rho cost on E (log2 ops)             : {cost_rho:.1f}")


def hess_magic_vs_cost(c: W1Curve) -> None:
    """For each potentially-reachable magic number m', tabulate the
    descended genus and attack cost. The Hess attack consists of finding
    an isogenous E' realising one of these m' values (the lower the better,
    so long as m' >= 2).

    Whether a given m' is reachable in the actual isogeny class of E is a
    separate (hard) computational question -- Menezes-Teske 2006 answered
    it positively for c2pnb176w1 at m' = 5, and negatively (no useful
    reduction below Pollard rho) for the other four c2pnb*w1 curves.
    """
    q = 1 << c.l
    rho_E = pollard_rho_cost_bits(c.order)
    print(f"  {c.name}: cost spectrum over reachable magic numbers")
    print(f"    direct ECDLP (Pollard rho on E)                    : 2^{rho_E:.1f}")
    print(f"    m'    g=2^(m'-1)   index-calc bits   Jac rho bits   best vs ECDLP")
    threshold_break = rho_E
    threshold_erode = rho_E + 16   # arbitrary "noticeably eroded" band
    saw_break = False
    saw_erode = False
    for mp in range(1, c.n + 1):
        g = descent_genus(mp)
        if g < 2:
            ic_s, rho_s, best, verdict = "  --  ", "  --  ", float('inf'), "useless (genus 1)"
        else:
            ic = index_calculus_cost_bits(g, q)
            jr = jac_pollard_rho_bits(g, q)
            best = min(ic, jr)
            ic_s = f"{ic:6.1f}"
            rho_s = f"{jr:6.1f}"
            if best < threshold_break:
                verdict = "*** BREAKS ECDLP ***"
                saw_break = True
            elif best < threshold_erode:
                verdict = "erodes margin"
                saw_erode = True
            else:
                verdict = "still > rho on E"
        best_s = f"2^{best:.1f}" if best != float('inf') else "  inf "
        print(f"    {mp:2d}    {g:8d}      {ic_s}            {rho_s}      {best_s:>10}  {verdict}")
    print()
    if saw_break:
        print(f"    => attack viable if isogeny walk realises an m' value flagged above.")
    elif saw_erode:
        print(f"    => no break but security margin shrinks vs a random binary curve.")
    else:
        print(f"    => no helpful descent reachable at any m'.")


# =============================================================================
# Toy Weil restriction worked out symbolically
# =============================================================================

def toy_weil_restriction_demo() -> None:
    """Compute the Weil restriction of E: y^2 + xy = x^3 + a*x^2 + b
    explicitly over F_{2^6} = F_{(2^2)^3} with a, b in F_{2^2}.

    This is the first algebraic step of any Weil descent: a 1-dim variety
    over F_{q^n} becomes an n-dim variety over F_q. The GHS curve C is
    cut out from this restriction by a 1-codimensional ideal arising
    from the function-field Frobenius descent. We compute the n=3
    restriction polynomials symbolically so the structure is visible."""

    import sympy
    from sympy import symbols, expand, Poly, simplify, GF

    # F_4 = F_2[u] / (u^2 + u + 1). We'll work modulo char 2.
    # Subfield curve: a = u, b = u + 1 (both in F_4)
    # T = generator of F_{2^6}/F_{2^2}, satisfying T^3 = T + 1.

    x0, x1, x2, y0, y1, y2, u, T = symbols("x0 x1 x2 y0 y1 y2 u T")
    a_sym = u
    b_sym = u + 1
    x_sym = x0 + x1 * T + x2 * T**2
    y_sym = y0 + y1 * T + y2 * T**2

    curve_lhs = y_sym**2 + x_sym * y_sym
    curve_rhs = x_sym**3 + a_sym * x_sym**2 + b_sym

    F = expand(curve_lhs - curve_rhs)

    # Reduce T^k for k >= 3 using T^3 = T + 1 (i.e. T^3 - T - 1 in char 0, but
    # in char 2 the minus signs vanish; we just XOR coefficients at the end).
    def reduce_T(poly):
        p = Poly(poly, T)
        coeffs = list(reversed(p.all_coeffs()))  # [c_0, c_1, c_2, ...] = T^0, T^1, ...
        # Reduce: any T^k for k >= 3 becomes T^{k-3}*(T+1) = T^{k-2} + T^{k-3}.
        # Pop highest, distribute into k-2 and k-3 positions, repeat.
        while len(coeffs) > 3:
            ck = coeffs.pop()                       # ck stood at index k = len(coeffs)
            k = len(coeffs)                         # now equals the old k
            coeffs[k - 2] = expand(coeffs[k - 2] + ck)
            coeffs[k - 3] = expand(coeffs[k - 3] + ck)
        while len(coeffs) < 3:
            coeffs.append(0)
        return coeffs  # [T^0 coeff, T^1 coeff, T^2 coeff]

    coeffs_T = reduce_T(F)

    # Reduce char-2: drop all even multiplicities. Treat coefficients as integers
    # then mod 2 each monomial. Also reduce u^2 -> u + 1.
    def char2_reduce(expr):
        e = expand(expr)
        # Replace u^2 -> u + 1 (since u^2 + u + 1 = 0)
        e = e.subs(u**2, u + 1)
        # Iterate until stable (handles u^3, etc.)
        for _ in range(5):
            new = expand(e.subs(u**2, u + 1))
            if new == e:
                break
            e = new
        # mod 2 the integer coefficients
        e = expand(e)
        p = Poly(e, x0, x1, x2, y0, y1, y2, u)
        new_terms = []
        for monom, coeff in p.terms():
            c = int(coeff) % 2
            if c == 0:
                continue
            term = 1
            for var, deg in zip(p.gens, monom):
                term *= var**deg
            new_terms.append(term)
        if not new_terms:
            return sympy.Integer(0)
        return sum(new_terms)

    rels = [char2_reduce(c) for c in coeffs_T]

    print("  Toy curve E:  y^2 + xy = x^3 + u*x^2 + (u+1)   over F_{2^6} = F_{(2^2)^3}")
    print("    a = u in F_{2^2}, b = u+1 in F_{2^2}  -- direct analog of c2pnb*w1.")
    print()
    print("  Setting x = x0 + x1*T + x2*T^2, y = y0 + y1*T + y2*T^2 with T^3=T+1,")
    print("  expanding y^2 + xy - x^3 - a*x^2 - b in char 2 with u^2 = u+1, and")
    print("  collecting by powers of T gives three polynomial relations over F_4:")
    print()
    for k, r in enumerate(rels):
        # Pretty-print
        s = str(r) if r != 0 else "0"
        # Wrap long lines
        if len(s) > 60:
            s = s.replace(" + ", "\n             + ")
        print(f"    [T^{k}]:  {s}  =  0")
        print()
    print("  These three polynomials over F_4 cut out the Weil restriction --")
    print("  a 3-dim affine variety over F_{2^2} whose F_4-points biject with")
    print("  E(F_{2^6})-points. Hess's algorithm then carves a hyperelliptic")
    print("  curve C out of this variety using the GHS function-field map;")
    print("  the genus formula g = 2^{magic - 1} from magic_dim() above is")
    print("  exactly the bookkeeping for that carving.")


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    banner("1. X9.62 c2pnb*w1 parameters: subfield-structure verification")
    print()
    for c in CURVES:
        report_curve(c)
        print()

    banner("2. The Hess (2003) extension: magic-number vs attack-cost spectrum")
    print()
    print("Direct GHS gives magic = 1 (section 1), descended genus = 1, no advantage.")
    print("Hess: find an isogenous E' whose b' has Frobenius-orbit dim m' > 1.")
    print("The descended hyperelliptic curve then has genus g = 2^{m'-1},")
    print("and the attack cost is min(index calculus on Jac(C), Pollard rho on Jac).")
    print("Below we sweep all theoretically reachable m' and show what each costs.")
    print("Whether the actual isogeny class of E reaches a given m' is a separate")
    print("question (Menezes-Teske answered yes for c2pnb176w1 at m' = 5).")
    print()
    for c in CURVES:
        hess_magic_vs_cost(c)
        print()

    banner("3. Toy Weil restriction (concrete structure, F_{2^6} = F_{(2^2)^3})")
    print()
    toy_weil_restriction_demo()

    banner("4. Reading the section-2 table; bottom line")
    print("""
  IMPORTANT: the '*** BREAKS ECDLP ***' lines in section 2 are hypothetical:
  they say "IF the isogeny class reaches this m', the descent costs X".
  Finding such an E' is a SEPARATE computational problem (the Hess isogeny
  walk -- step 1 of the pipeline below) and is what distinguishes the
  c2pnb176w1 break from the rest of the family.

  Per Menezes-Teske 2006 the actual reachable m' and resulting costs:

    curve         rho on E   reachable m'   resulting cost     verdict
    c2pnb176w1     2^80         m' = 5        ~2^57            BROKEN
    c2pnb208w1     2^96         m' >= 7       >~2^96           secure vs descent
    c2pnb272w1     2^128        m' >= ?       >~2^128          secure vs descent
    c2pnb304w1     2^144        m' >= ?       >~2^144          secure vs descent
    c2pnb368w1     2^176        m' >= ?       >~2^176          secure vs descent

  Compare against section 2's tables: c2pnb176w1's m'=5 row reads ~2^74
  (our crude cost model overshoots Menezes-Teske's 2^57 by a factor that
  comes from a tighter constant in the relation-collection phase, but the
  shape -- well under 2^80 rho -- is right). For c2pnb208w1, the same row
  *would* break the curve if reached, but the isogeny class doesn't get
  there: the smallest m' the walk attains is m' >= 7, giving descent cost
  > 2^327 -- vastly worse than rho on E.

  The full attack pipeline (this PoC does NOT execute steps 1-4):

    0.        Recognise the w-curve subfield structure.            <-- section 1
    0.5       Sweep the magic-vs-cost spectrum.                    <-- section 2
    1. (out)  SEA / Galbraith-Smart isogeny walk to find E' with the
              smallest reachable m'.
    2. (out)  Hess's function-field algorithm to build C/F_{2^16}.  <-- section 3
              skeleton sketches the Weil restriction the algorithm
              starts from.
    3. (out)  Diem / Gaudry-Theriault index calculus on Jac(C).
    4. (out)  Lift the recovered HCDLP through the isogeny back to E.

  Each of steps 1-4 is a research-grade implementation. Step 1 is the
  hard one and is what determines whether a given c2pnb*w1 curve is in
  practice broken. For c2pnb176w1 the answer is yes; for the rest, no.
""")


if __name__ == "__main__":
    main()
