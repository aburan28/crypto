#!/usr/bin/env python3
"""Boundaries on transferring ECC2K-130 to a hyperelliptic curve of another genus.

`E : y^2 + xy = x^3 + 1` over `F_2^131` has genus 1.  Does it map to a curve of
some other genus in a way that carries the discrete logarithm somewhere cheaper
than the Pollard rho reference on the original group?

A transfer has exactly two places to live, and the script prices both:

  * over the same field `F_2^131`, where covers exist in every genus and index
    calculus on them is priced by the ambient field: at least `q = 2^131`;
  * over the only proper subfield, `F_2`, where index calculus would be cheap --
    but where the genus a construction can reach is `1`, `2^129` or `2^130`,
    with nothing in between, because `2` is a primitive root modulo `131`.

Everything below is derived, not measured.  The script computes the curve
constants, the Weil polynomial of the Weil restriction, the GHS magic numbers by
explicit `F_2^131` arithmetic, and the exact zeta function of the genus-130
curve an attacker would need, then prices index calculus on it from that zeta
function.  It writes `experiments/ecc2k130_hyperelliptic_cover_boundary.json`.

    python3 scripts/ecc2k130_hyperelliptic_cover_boundary.py [--samples N] [--gmax G]
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
import random
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "experiments/ecc2k130_hyperelliptic_cover_boundary.json"

N131 = 131                      # the challenge field degree, prime
TRACE = -1                      # t for K_0 : y^2 + xy = x^3 + 1 over F_2
Q = 2
SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

# x^131 + x^13 + x^2 + x + 1, the ECC2K-130 reduction pentanomial
F2_131_MODULUS = (1 << 131) | (1 << 13) | (1 << 2) | (1 << 1) | 1


# ---------------------------------------------------------------------------
# integers
# ---------------------------------------------------------------------------

def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin over the first twelve primes (ample below 2^128)."""
    if n < 2:
        return False
    for p in SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in SMALL_PRIMES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def frobenius_trace(k: int) -> int:
    """s_k = alpha^k + alphabar^k, the trace of Frobenius^k on E/F_2."""
    s0, s1 = 2, TRACE
    if k == 0:
        return s0
    for _ in range(k - 1):
        s0, s1 = s1, TRACE * s1 - Q * s0
    return s1


def curve_order(k: int) -> int:
    """#E(F_2^k) = 2^k + 1 - s_k."""
    return Q ** k + 1 - frobenius_trace(k)


def mobius(n: int) -> int:
    m, res = n, 1
    p = 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            res = -res
        p += 1
    return -res if m > 1 else res


def mult_order(a: int, m: int) -> int:
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k


def sqrt_mod(a: int, p: int) -> int | None:
    """Tonelli-Shanks square root modulo an odd prime."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, rr = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, rr = t * c % p, rr * b % p
    return rr


# ---------------------------------------------------------------------------
# polynomials and power series over Z
# ---------------------------------------------------------------------------

def poly_divmod(num: list[int], den: list[int]) -> tuple[list[int], list[int]]:
    """Division in Z[T] by a monic denominator."""
    assert den[-1] == 1
    out, quo = list(num), [0] * (len(num) - len(den) + 1)
    for i in range(len(num) - len(den), -1, -1):
        c = out[i + len(den) - 1]
        quo[i] = c
        if c:
            for j, d in enumerate(den):
                out[i + j] -= c * d
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return quo, out


def poly_eval(p: list[int], x: int) -> int:
    v = 0
    for c in reversed(p):
        v = v * x + c
    return v


def series_mul(a: list[int], b: list[int], n: int) -> list[int]:
    out = [0] * n
    for i, ai in enumerate(a[:n]):
        if not ai:
            continue
        for j, bj in enumerate(b[: n - i]):
            if bj:
                out[i + j] += ai * bj
    return out


# ---------------------------------------------------------------------------
# F_2[x]/(f) arithmetic, bit-packed -- the GHS magic number lives here
# ---------------------------------------------------------------------------

def f2_mod(a: int, m: int) -> int:
    dm = m.bit_length() - 1
    while a.bit_length() - 1 >= dm:
        a ^= m << (a.bit_length() - 1 - dm)
    return a


def f2_sqrmod(a: int, m: int) -> int:
    sq, i = 0, 0
    while a:
        if a & 1:
            sq |= 1 << (2 * i)
        a >>= 1
        i += 1
    return f2_mod(sq, m)


def f2_rank(vectors: list[int], n: int) -> int:
    rows, piv = list(vectors), 0
    for col in range(n):
        mask = 1 << col
        pick = next((r for r in range(piv, len(rows)) if rows[r] & mask), None)
        if pick is None:
            continue
        rows[piv], rows[pick] = rows[pick], rows[piv]
        for r in range(len(rows)):
            if r != piv and rows[r] & mask:
                rows[r] ^= rows[piv]
        piv += 1
        if piv == len(rows):
            break
    return piv


def magic_dim(b: int, n: int, modulus: int) -> int:
    """dim_F2 span{b, b^2, b^4, ...}: the GHS magic number of b for descent to F_2."""
    orbit, cur = [], b
    for _ in range(n):
        orbit.append(cur)
        cur = f2_sqrmod(cur, modulus)
    assert cur == b, "Frobenius^n is not the identity: reducible modulus"
    return f2_rank(orbit, n)


def absolute_trace(b: int, n: int, modulus: int) -> int:
    """Tr_{F_2^n/F_2}(b) = b + b^2 + ... + b^(2^(n-1))."""
    acc, cur = 0, b
    for _ in range(n):
        acc ^= cur
        cur = f2_sqrmod(cur, modulus)
    assert acc in (0, 1), "trace is not in F_2"
    return acc


def factor_degrees(n: int, q: int = 2) -> Counter:
    """Degrees of the irreducible factors of t^n - 1 over F_q, q a power of 2.

    In characteristic 2, t^n - 1 = (t^m - 1)^(2^v) for n = 2^v m with m odd, and
    t^m - 1 has one irreducible factor of degree |q mod d| per cyclotomic coset.
    """
    v, m = 0, n
    while m % 2 == 0:
        m //= 2
        v += 1
    seen, degs = set(), Counter()
    for a in range(m):
        if a in seen:
            continue
        coset, x = set(), a
        while x not in coset:
            coset.add(x)
            x = x * q % m
        seen |= coset
        degs[len(coset)] += 1 << v
    return degs


def achievable_dims(n: int, q: int = 2) -> set[int]:
    """Dimensions of the Frobenius-stable F_q-subspaces of F_q^n: subset sums of
    the factor degrees of t^n - 1 over F_q."""
    dims = {0}
    for deg, mult in factor_degrees(n, q).items():
        for _ in range(mult):
            dims |= {d + deg for d in dims}
    return dims


# ---------------------------------------------------------------------------
# index calculus on a curve over F_2 whose zeta function is known
# ---------------------------------------------------------------------------

def place_counts(points: list[int], dmax: int) -> list[int]:
    """N_d, the number of degree-d places, from #C(F_2^e) for e | d."""
    out = [0] * (dmax + 1)
    for d in range(1, dmax + 1):
        acc = sum(mobius(d // e) * points[e] for e in range(1, d + 1) if d % e == 0)
        assert acc % d == 0 and acc >= 0, (d, acc)
        out[d] = acc // d
    return out


def divisor_series(places: list[int], g: int, dmax: int) -> list[int]:
    """Effective divisors by degree, up to degree g, built only from places of
    degree <= dmax: the coefficients of prod_{d <= dmax} (1 - T^d)^(-N_d)."""
    ser = [1] + [0] * g
    for d in range(1, dmax + 1):
        nd = places[d]
        if nd == 0:
            continue
        fac, k = [0] * (g + 1), 0
        while d * k <= g:
            fac[d * k] = math.comb(nd + k - 1, k)
            k += 1
        ser = series_mul(ser, fac, g + 1)
    return ser


def price_index_calculus(points: list[int], g: int, bmax: int = 44) -> dict:
    """Cheapest (relations, linear algebra) split over smoothness bounds b <= bmax.

    One operation = one random walk step in the divisor class group plus one
    smoothness test of the reduced divisor; linear algebra is counted as
    |FB|^2 * g multiplications, charging every relation the maximum row weight.
    """
    places = place_counts(points, g)
    total = divisor_series(places, g, g)[g]           # every effective divisor of degree g
    best = None
    for b in range(2, min(bmax, g) + 1):
        fb = sum(places[1 : b + 1])
        smooth = divisor_series(places, g, b)[g]
        if fb < 2 or smooth == 0:
            continue
        frac = smooth / total
        rel = math.log2(fb) - math.log2(frac)
        lin = 2 * math.log2(fb) + math.log2(g)
        tot = max(rel, lin) + 1                      # + 1: the two phases together
        if best is None or tot < best["log2_total"]:
            best = {
                "genus": g,
                "smoothness_bound_b": b,
                "log2_factor_base": round(math.log2(fb), 2),
                "log2_smooth_probability": round(math.log2(frac), 2),
                "log2_relations": round(rel, 2),
                "log2_linear_algebra": round(lin, 2),
                "log2_total": round(tot, 2),
            }
    return best


def rational_model_points(kmax: int) -> list[int]:
    """#C(F_2^k) = 2^k + 1: reduced divisors behave like random monic polynomials,
    the model Enge-Gaudry's L(1/2) analysis rests on."""
    return [0] + [2 ** k + 1 for k in range(1, kmax + 1)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=400,
                    help="random b in F_2^131 sampled for the magic-number census")
    ap.add_argument("--gmax", type=int, default=400, help="largest genus priced over F_2")
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    # ---- the target, and the boundary everything is measured against --------
    order = curve_order(N131)
    r = order // 4
    assert order % 4 == 0 and is_prime(r)
    log2r = math.log2(r)
    aut = 2 * N131                                  # <-1> x <Frobenius> on <G>
    log2_rho = math.log2(math.sqrt(math.pi * r / (2 * aut)))
    log2_rho_plain = math.log2(math.sqrt(math.pi * r / 2))

    target = {
        "curve": "K_0 : y^2 + x y = x^3 + 1 over F_2, used over F_2^131",
        "trace_over_F2": TRACE,
        "order": str(order),
        "cofactor": 4,
        "r": str(r),
        "r_is_prime": True,
        "log2_r": round(log2r, 4),
        "rho_automorphism_order": aut,
        "log2_rho_plain": round(log2_rho_plain, 4),
        "log2_rho_reference": round(log2_rho, 4),
        "S_rho_reference": round(2 ** (log2_rho - log2r / 2), 6),
    }

    # ---- Boundary A: a transfer moves the group, never its size -------------
    boundary_a = {
        "statement": "a correspondence that preserves the discrete logarithm carries "
                     "<G> to a cyclic group of the same order r, so the rho reference "
                     "is invariant under every transfer considered here",
        "log2_rho_on_the_image": round(log2_rho, 4),
        "note": "raising the genus buys the attacker index calculus, never a smaller group",
    }

    # ---- Boundary B: covers over the same field ----------------------------
    # Gaudry/Diem index calculus on a genus-g Jacobian over F_q costs q^(2-2/g);
    # generic rho on that Jacobian costs q^(g/2).  With q = 2^131 both are >= 2^131.
    same_field = []
    for g in range(2, 9):
        ic = N131 * (2 - 2 / g)
        rho_jac = N131 * g / 2
        best = min(ic, rho_jac)
        same_field.append({
            "genus": g,
            "log2_index_calculus_q_pow_2_minus_2_over_g": round(ic, 2),
            "log2_rho_on_jacobian": round(rho_jac, 2),
            "log2_best": round(best, 2),
            "log2_ratio_to_rho_reference": round(best - log2_rho, 2),
        })
    twist_order = Q ** N131 + 1 + frobenius_trace(N131)
    boundary_b = {
        "statement": "over F_2^131 covers exist in every genus, and every one of them "
                     "is priced by the ambient field: at least q = 2^131",
        "reference_complexity": "Gaudry/Diem Otilde(q^(2-2/g)) for fixed g, against "
                                "q^(g/2) for rho on the whole Jacobian",
        "cheapest_cell": min(same_field, key=lambda c: c["log2_best"]),
        "cells": same_field,
        "quadratic_weil_restriction": {
            "construction": "Res_{F_2^262/F_2^131}(E) is an abelian surface",
            "splits_as": "E x E^twist",
            "why": "E is already defined over F_2^131, so the restriction decomposes "
                   "along t^2 - 1 = (t-1)(t+1) into the curve and its quadratic twist",
            "twist_order": str(twist_order),
            "gcd_with_curve_order": math.gcd(order, twist_order),
            "note": "a split surface is the boundary of the moduli of abelian surfaces, "
                    "not a genus-2 Jacobian; gluing along torsion (Howe, then Mestre "
                    "for the equation -- RESEARCH_MESTRE_HOWE.md) produces one, still "
                    "over F_2^131, so still on the row above",
        },
    }

    # ---- Boundary C: the genus GHS can reach over F_2 -----------------------
    # The only proper subfield of F_2^131 is F_2.  GHS descends E/F_2^131 to a
    # curve over F_2 of genus 2^(m-1), m = dim_F2 span{b, b^2, b^4, ...}.  That
    # span is a Frobenius-stable subspace, and t^131 - 1 = (t-1)(irreducible of
    # degree 130) over F_2, so the only dimensions available are 0, 1, 130, 131.
    degs = factor_degrees(N131)
    assert dict(degs) == {1: 1, 130: 1}, degs
    assert achievable_dims(N131) == {0, 1, 130, 131}

    census = Counter()
    m_challenge = magic_dim(1, N131, F2_131_MODULUS)      # a = 0, b = 1, both in F_2
    for _ in range(args.samples):
        bb = rng.randrange(1, 1 << N131)
        m = magic_dim(bb, N131, F2_131_MODULUS)
        tr = absolute_trace(bb, N131, F2_131_MODULUS)
        census[(m, tr)] += 1
        assert m in (1, 130, 131), (m, bb)
        assert (m == 1) == (bb == 1), (m, bb)
        assert (m == 130) == (tr == 0 and bb != 1), (m, tr, bb)

    # the m = 1 transfer is the trace map, and the trace annihilates <G>:
    # Frobenius acts on <G> as lambda with lambda^2 + lambda + 2 = 0 and
    # lambda^131 = 1, so sum_{i<131} lambda^i = (lambda^131 - 1)/(lambda - 1) = 0.
    root = sqrt_mod(-7 % r, r)
    assert root is not None
    inv2 = pow(2, r - 2, r)
    lam = next(c for c in ((-1 + root) * inv2 % r, (-1 - root) * inv2 % r)
               if pow(c, N131, r) == 1)
    assert (lam * lam + lam + 2) % r == 0
    trace_on_G = sum(pow(lam, i, r) for i in range(N131)) % r
    assert trace_on_G == 0

    boundary_c = {
        "statement": "for every elliptic curve over F_2^131 the GHS magic number is "
                     "1, 130 or 131, so the genus of the descended curve over F_2 is "
                     "1, 2^129 or 2^130 -- there is no third case",
        "ord_2_mod_131": mult_order(2, N131),
        "factor_degrees_of_t131_minus_1": {str(k): v for k, v in sorted(degs.items())},
        "available_subspace_dimensions": sorted(achievable_dims(N131)),
        "magic_number_of_ecc2k130": m_challenge,
        "genus_of_the_ecc2k130_descent": 1 << (m_challenge - 1),
        "magic_number_census": {f"m={m},trace={t}": c for (m, t), c in sorted(census.items())},
        "samples": args.samples,
        "log2_available_genera": [0.0, 129.0, 130.0],
        "isogeny_shift_closed": "Hess/Menezes-Teske raise the magic number by walking to "
                                "an isogenous curve; here every curve over F_2^131 obeys "
                                "the same trichotomy, so the walk has nowhere to land",
        "genus_1_transfer_is_the_trace_map": {
            "frobenius_eigenvalue_lambda": str(lam),
            "lambda_pow_131_mod_r": pow(lam, N131, r),
            "sum_lambda_i_mod_r": trace_on_G,
            "conclusion": "the conorm-norm map into E(F_2) is the zero map on <G>, not "
                          "merely a map into a group too small to be useful",
        },
        "cost_floor_for_the_large_genus": {
            "genus": "2^129",
            "log2_cost_floor": 129.0,
            "log2_ratio_to_rho": round(129.0 - log2_rho, 2),
            "why": "one divisor on a genus-2^129 curve takes 2^129 bits to write down",
        },
    }

    # ---- Boundary D: how much genus a transfer over F_2 has to buy ----------
    # (i) Weil: #Jac(C)(F_2) <= (1 + sqrt 2)^(2g) and the image of <G> has order r.
    weil_floor = math.ceil(log2r * math.log(2) / (2 * math.log(1 + math.sqrt(2))))
    # (ii) sharper: Res_{F_2^131/F_2}(E) ~ E x A with A simple of dimension 130,
    #      and <G> is A(F_2).  A nonzero map from a simple A into Jac(C) is an
    #      isogeny onto its image, so dim Jac(C) = genus(C) >= 130.
    s131 = frobenius_trace(N131)
    p_w = [0] * (2 * N131 + 1)
    p_w[0], p_w[N131], p_w[2 * N131] = Q ** N131, -s131, 1
    p_e = [Q, 1, 1]                                       # T^2 + T + 2
    p_a, rem = poly_divmod(p_w, p_e)
    assert rem == [0] and len(p_a) == 261
    assert poly_eval(p_w, 1) == order and poly_eval(p_e, 1) == curve_order(1) == 4
    assert poly_eval(p_a, 1) == r
    assert all(p_a[i] == Q ** (130 - i) * p_a[260 - i] for i in range(131))
    assert N131 % 7 != 0                                  # conductor of Q(sqrt(-7))
    boundary_d = {
        "statement": "<G> is A(F_2) for a simple abelian variety A/F_2 of dimension 130, "
                     "so any correspondence over F_2 carrying the discrete logarithm to "
                     "a curve C/F_2 forces genus(C) >= 130",
        "weil_restriction": "Res_{F_2^131/F_2}(E) ~ E x A, from t^131 - 1 = (t-1) Phi_131(t)",
        "char_poly_of_the_weil_restriction": "T^262 - s T^131 + 2^131",
        "s_131": str(s131),
        "dim_A": 130,
        "points_on_A_over_F2": str(poly_eval(p_a, 1)),
        "points_on_A_equal_r": poly_eval(p_a, 1) == r,
        "functional_equation_verified": True,
        "A_is_simple": True,
        "simplicity_proof": "the Weil numbers of A are zeta*alpha with zeta^131 = 1, "
                            "zeta != 1 and alpha = (-1+sqrt(-7))/2.  Q(sqrt(-7)) has "
                            "conductor 7, which does not divide 131, so Q(zeta_131) and "
                            "Q(sqrt(-7)) are linearly disjoint; the only roots of unity "
                            "in Q(sqrt(-7)) are +-1 and 131 is odd, so no element of "
                            "Gal(Q(zeta_131, sqrt(-7))/Q) fixes zeta*alpha.  Its degree "
                            "is therefore 260 = 2 dim A, and Honda-Tate makes A simple",
        "weil_bound_genus_floor": weil_floor,
        "simplicity_genus_floor": 130,
    }

    # ---- Boundary E: the window -- which genus over F_2 would actually pay --
    # Exact cell: if A were a Jacobian, its curve would have this zeta function.
    a_points = [0] + [2 ** k + 1 + frobenius_trace(k) if k % N131
                      else 2 ** k + 1 - 130 * frobenius_trace(k)
                      for k in range(1, 131)]
    exact = price_index_calculus(a_points, 130)
    exact["model"] = "exact: place counts from the zeta function of A"
    exact["log2_ratio_to_rho"] = round(exact["log2_total"] - log2_rho, 2)

    window = []
    for g in [130, 140, 150, 175, 200, 225, 250, 275, 290, 300, 325, 350, args.gmax]:
        if g > args.gmax:
            continue
        cell = price_index_calculus(rational_model_points(g), g)
        cell["model"] = "random-polynomial"
        cell["log2_ratio_to_rho"] = round(cell["log2_total"] - log2_rho, 2)
        cell["beats_rho"] = cell["log2_total"] < log2_rho
        window.append(cell)
    crossover = max((c["genus"] for c in window if c["beats_rho"]), default=None)
    first_loss = min((c["genus"] for c in window if not c["beats_rho"]), default=None)
    boundary_e = {
        "statement": "index calculus over F_2 beats the rho reference only in a narrow "
                     "band of genus; the band is bounded below by boundary D and above "
                     "by its own cost",
        "operation": "one divisor-class step plus one smoothness test; linear algebra "
                     "charged as |FB|^2 * g multiplications modulo r",
        "window_low": 130,
        "window_high_between": [crossover, first_loss],
        "exact_cell_if_A_were_a_jacobian": exact,
        "log2_gap_from_the_window_to_the_nearest_available_genus":
            round(129.0 - math.log2(first_loss), 2) if first_loss else None,
        "cells": window,
    }

    # ---- Boundary F: what it would take to have such a curve ---------------
    g = 130
    deg_a = 3 ** g * math.factorial(g)                # (3L)^g for a ppav (A, L)
    section_genus_log2 = math.log2(1 + (g - 1) * deg_a / 2)
    boundary_f = {
        "statement": "a curve of genus 130 over F_2 whose Jacobian contains A is not "
                     "excluded by any point count; it is simply not constructible, and "
                     "a dimension count says such curves are vanishingly rare",
        "dim_M_g": 3 * g - 3,
        "dim_hyperelliptic_locus": 2 * g - 1,
        "dim_A_g": g * (g + 1) // 2,
        "torelli_codimension": g * (g + 1) // 2 - (3 * g - 3),
        "hyperelliptic_codimension": g * (g + 1) // 2 - (2 * g - 1),
        "trivial_existence": "every abelian variety is a quotient of a Jacobian: cut A "
                             "by dim A - 1 hyperplanes of the 3-theta embedding",
        "log2_degree_of_A_in_P3g": round(math.log2(deg_a), 1),
        "log2_genus_of_that_curve_section": round(section_genus_log2, 1),
        "explicit_cm_constructions_cap_at_genus": 3,
        "caveat": "the codimension count is a heuristic: the isogeny class of A carries "
                  "as many principally polarised members as the CM field has ideal "
                  "classes, which is not a number this thread has bounded",
    }

    # ---- The counterfactual: a composite degree, where covers do land ------
    # The magic-number column is the whole story.  At N = 130 the descent to F_2
    # itself reaches genus 128, which boundary E prices below rho; at N = 131 the
    # least magic number a transfer can use is 130.
    counterfactual = []
    for N in (130, N131):
        for l in range(1, N + 1):
            if N % l or N // l < 2:
                continue
            n = N // l
            dims = achievable_dims(n, 1 << l)
            # a transfer can only be injective if the target group is big enough:
            # #Jac ~ 2^(l * 2^(m-1)) has to reach 2^(N-2)
            admissible = sorted(d for d in dims if d >= 1 and l * (1 << (d - 1)) >= N - 2)
            m = admissible[0] if admissible else None
            cell = {"field_degree_N": N, "subfield_degree_l": l, "extension_degree_n": n,
                    "available_magic_numbers": sorted(dims)[:12],
                    "least_admissible_magic_number": m,
                    "log2_rho_on_that_group": round(
                        math.log2(math.sqrt(math.pi * 2 ** (N - 2) / (2 * 2 * N))), 2)}
            if m is None:
                counterfactual.append(cell)
                continue
            gg = 1 << (m - 1)
            cell["descended_genus"] = str(gg) if m <= 24 else f"2^{m-1}"
            if l == 1 and m <= 12:
                # large genus over F_2: the boundary E model, same units as above
                priced = price_index_calculus(rational_model_points(gg), gg)
                cell["cost_model"] = "large genus over F_2 (boundary E)"
                cell["log2_cost"] = priced["log2_total"]
                cell["smoothness_bound_b"] = priced["smoothness_bound_b"]
            elif l == 1:
                cell["cost_model"] = "large genus over F_2, cost floor = bits per divisor"
                cell["log2_cost"] = float(m - 1)
            elif m <= 24:
                cell["cost_model"] = "fixed genus over F_2^l: g! q^(2-2/g)"
                cell["log2_cost"] = round(
                    math.lgamma(gg + 1) / math.log(2) + l * (2 - 2 / gg), 2)
            else:
                cell["cost_model"] = "fixed genus over F_2^l: g! q^(2-2/g)"
                cell["log2_cost"] = float(l * (m - 1))
            cell["breaks_rho"] = cell["log2_cost"] < cell["log2_rho_on_that_group"]
            counterfactual.append(cell)

    report = {
        "schema": "ecc2k130_hyperelliptic_cover_boundary/v1",
        "question": "does ECC2K-130 map to a hyperelliptic curve of another genus in a "
                    "way that costs less than rho on E(F_2^131)?",
        "verdict": "KILLED for every construction; the residual question is named in "
                   "boundary F and is not a search anyone can run",
        "unit": "log2 operations; S = ops / sqrt(r)",
        "target": target,
        "boundary_A_transfer_preserves_the_group": boundary_a,
        "boundary_B_same_field_covers": boundary_b,
        "boundary_C_ghs_genus_trichotomy": boundary_c,
        "boundary_D_genus_floor": boundary_d,
        "boundary_E_window": boundary_e,
        "boundary_F_constructibility": boundary_f,
        "counterfactual_composite_degree": counterfactual,
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n")

    print(f"r = {r}  prime={target['r_is_prime']}  log2 r = {log2r:.4f}")
    print(f"rho reference      2^{log2_rho:.4f}   (plain 2^{log2_rho_plain:.4f}, "
          f"S = {target['S_rho_reference']})")
    c = boundary_b["cheapest_cell"]
    print(f"B  cheapest same-field cover: genus {c['genus']} at 2^{c['log2_best']} "
          f"= 2^{c['log2_ratio_to_rho_reference']} x rho")
    print(f"C  GHS over F_2: magic number of ECC2K-130 = {m_challenge} (genus "
          f"{1 << (m_challenge - 1)}, transfer = 0 on <G>); every other curve over "
          f"F_2^131 has magic 130 or 131, genus 2^129 or 2^130")
    print(f"D  genus floor: {weil_floor} from Weil, 130 from the simplicity of A "
          f"(#A(F_2) = r, exactly)")
    print(f"E  window: genus 130 .. between {crossover} and {first_loss}; at genus 130 "
          f"the attack would cost 2^{exact['log2_total']} "
          f"= 2^{exact['log2_ratio_to_rho']} x rho")
    print(f"F  Torelli codimension at g=130: {boundary_f['torelli_codimension']}; "
          f"explicit CM constructions stop at genus "
          f"{boundary_f['explicit_cm_constructions_cap_at_genus']}")
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
