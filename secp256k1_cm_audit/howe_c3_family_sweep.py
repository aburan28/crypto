#!/usr/bin/env python3
"""
howe_c3_family_sweep.py  --  Thread 22 (autolab 2026-07-29)

QUESTION
--------
Thread 18 (log 2026-07-21) found that 5 of the 15 pairs of j=0 sextic twists of
secp256k1 satisfy the Howe gluing conditions (H1)(H2)(H3).  Thread 2 (log
2026-07-26) then showed that the *Rosenhain* route to writing down the glued
genus-2 curve is blocked for secp256k1, because it needs F_p-rational 2-torsion
(equivalently: -b a cubic residue), which secp256k1 does not have.

That result is about a *formula*, not about *existence*.  This script asks the
existence question directly, over small proxy primes:

    For p = 1 mod 6 and E_k : y^2 = x^3 + g^k  (k=0..5, g a primitive root,
    the six sextic-twist classes), does there EXIST a genus-2 curve C/F_p with

        charpoly(Frob | Jac C) = (T^2 - t_i T + p)(T^2 - t_j T + p)

    for a Howe-qualifying pair (i,j)?

This script answers it inside the zeta_3-equivariant family

        C_{a,c} : y^2 = lc * (x^6 + a x^3 + c),   lc in {1, non-residue}

which is exactly the family of genus-2 curves carrying the order-3 automorphism
(x,y) -> (zeta_3 x, y) inherited from the zeta_3 automorphisms of the two j=0
factors.  The family is enumerated EXHAUSTIVELY: all (a,c) in F_p x F_p*,
both quadratic twists.  |family| = 2 p (p-1).

Note y^2 = (x^3 + b_i)(x^3 + b_j) is the member (a,c) = (b_i+b_j, b_i b_j), so
the "naive cover" of RESEARCH_MESTRE_HOWE.md / chlrs_igusa_formula.gp is a
single point of the search space swept here.

METHOD
------
hyperellcharpoly() in PARI is far too slow for ~10^5 curves (measured: >120 s
for p=43 alone).  We instead recover the degree-4 Weil polynomial

    P(T) = T^4 - c1 T^3 + c2 T^2 - p c1 T + p^2

from #C(F_p) and #C(F_p^2) only, vectorised over the whole family with numpy.

Write S1 = sum_{x in F_p} chi_p(f(x)),  S2 = sum_{X in F_p^2} chi_q(f(X)),
with f(x) = x^6 + a x^3 + c (monic, so 2 points at infinity in both fields):

    #C(F_p)   = p + 2 + S1        =>  c1 = p + 1 - #C(F_p)   = -1 - S1
    #C(F_p^2) = p^2 + 2 + S2      =>  c2 = (c1^2 + 1 + S2) / 2

Two tricks make S2 cheap.  Fix F_p^2 = F_p[t]/(t^2 - r), r a non-residue.
 (i) For w = A + B t, w is a square in F_p^2 iff its norm A^2 - r B^2 is a
     square in F_p.  So chi_q(A,B) = chi_p(A^2 - r B^2): no F_p^2 arithmetic.
 (ii) With U = X^3 precomputed once per prime, f(X) = U^2 + a U + c, so the
     per-curve work is two multiply-adds over the p^2 grid.

The quadratic twist lc = r is free: chi_p(r)= -1 gives S1 -> -S1, while
chi_q(r) = chi_p(r^2) = +1 gives S2 -> S2, and the infinity count flips 2 -> 0
over F_p and stays 2 over F_p^2.  Net effect (c1, c2) -> (-c1, c2), the
expected action of quadratic twisting on the Weil polynomial.

The Weil-polynomial arithmetic is checked against PARI hyperellcharpoly on a
sample of curves by --selftest.

USAGE
-----
    python3 howe_c3_family_sweep.py                # default prime list
    python3 howe_c3_family_sweep.py --selftest     # cross-check vs PARI
    python3 howe_c3_family_sweep.py -p 43,61,67
"""

import argparse
import subprocess
import sys

import numpy as np

# --------------------------------------------------------------------------
# small helpers over F_p
# --------------------------------------------------------------------------


def legendre_table(p):
    """chi[v] = Legendre symbol (v|p) for v in 0..p-1, as int8."""
    chi = np.full(p, -1, dtype=np.int8)
    chi[0] = 0
    sq = (np.arange(1, p, dtype=np.int64) ** 2) % p
    chi[sq] = 1
    return chi


def primitive_root(p):
    """Smallest primitive root mod p (p prime, small)."""
    fac = set()
    m = p - 1
    d = 2
    while d * d <= m:
        while m % d == 0:
            fac.add(d)
            m //= d
        d += 1
    if m > 1:
        fac.add(m)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g
    raise ValueError(f"no primitive root for {p}")


def non_residue(p):
    """Smallest quadratic non-residue mod p."""
    for r in range(2, p):
        if pow(r, (p - 1) // 2, p) == p - 1:
            return r
    raise ValueError(f"no non-residue for {p}")


def ec_trace(p, b):
    """Trace of Frobenius of y^2 = x^3 + b over F_p, by naive point count."""
    chi = legendre_table(p)
    x = np.arange(p, dtype=np.int64)
    npts = 1 + p + int(chi[(x * x % p * x + b) % p].sum())  # +1 for infinity
    return p + 1 - npts


def cubic_split_pattern(p, b):
    """Degree pattern of x^3 + b over F_p: (1,1,1) if it splits, (3,) if not.

    Over F_p with 3 | p-1 there is no partial splitting: x^3 + b either has
    three roots (iff -b is a cubic residue) or none.
    """
    roots = [x for x in range(p) if (x * x % p * x + b) % p == 0]
    return (1, 1, 1) if roots else (3,)


# --------------------------------------------------------------------------
# the six sextic twists and their Howe-qualifying pairs
# --------------------------------------------------------------------------


def sextic_twists(p):
    """The 6 sextic-twist classes E_k : y^2 = x^3 + g^k, with traces / orders.

    For p = 1 mod 6 the classes of j=0 curves are b * (F_p*)^6, and g^k for
    k = 0..5 hits each class exactly once (g a primitive root).
    """
    g = primitive_root(p)
    out = []
    for k in range(6):
        b = pow(g, k, p)
        t = ec_trace(p, b)
        out.append(
            {
                "k": k,
                "b": b,
                "t": t,
                "n": p + 1 - t,
                "tors2": cubic_split_pattern(p, b),
            }
        )
    return out


def howe_pairs(twists, p):
    """All 15 pairs with the (H1)(H2)(H3) verdicts of Thread 18.

    H1: #E_i != #E_j          (the two factors are not F_p-isogenous)
    H2: same 2-torsion Galois-module type (same degree pattern of x^3+b)
    H3: gcd(#E_i, #E_j) = 1   (no shared rational torsion to obstruct gluing)
    """
    from math import gcd

    rows = []
    for i in range(6):
        for j in range(i + 1, 6):
            ei, ej = twists[i], twists[j]
            h1 = ei["n"] != ej["n"]
            h2 = ei["tors2"] == ej["tors2"]
            h3 = gcd(ei["n"], ej["n"]) == 1
            rows.append(
                {
                    "pair": (i, j),
                    "H1": h1,
                    "H2": h2,
                    "H3": h3,
                    "ok": h1 and h2 and h3,
                    # target Weil polynomial coefficients
                    "c1": ei["t"] + ej["t"],
                    "c2": 2 * p + ei["t"] * ej["t"],
                }
            )
    return rows


# --------------------------------------------------------------------------
# exhaustive sweep of the zeta_3-equivariant family
# --------------------------------------------------------------------------


def c3_family_weil(p, chunk=256):
    """Weil coefficients (c1, c2) of y^2 = x^6 + a x^3 + c for all (a,c).

    Returns (avals, cvals, c1, c2, ok) as flat arrays over the grid
    a in 0..p-1, c in 1..p-1, with `ok` False where the curve is singular
    (i.e. x^6 + a x^3 + c not squarefree).

    The quadratic twist is not listed separately: it is (-c1, c2).
    """
    chi = legendre_table(p).astype(np.int64)
    r = non_residue(p)

    # --- grid of curves -----------------------------------------------
    A_, C_ = np.meshgrid(
        np.arange(p, dtype=np.int64), np.arange(1, p, dtype=np.int64), indexing="ij"
    )
    avals = A_.ravel()
    cvals = C_.ravel()
    ncur = avals.size

    # --- singularity test ---------------------------------------------
    # x^6 + a x^3 + c = (x^3 - r1)(x^3 - r2) with r1 r2 = c, r1 + r2 = -a.
    # c != 0 already, so 0 is not a root.  The curve is singular iff r1 == r2,
    # i.e. iff the discriminant a^2 - 4c vanishes.  (r1 != r2 with both
    # non-zero always gives six distinct roots in F_p-bar, since x^3 - r1 and
    # x^3 - r2 are separable and coprime.)
    ok = ((avals * avals - 4 * cvals) % p) != 0

    # --- S1: sum over F_p ----------------------------------------------
    # u = x^3 ranges over F_p with multiplicity; just evaluate at every x.
    x = np.arange(p, dtype=np.int64)
    U1 = (x * x % p) * x % p  # x^3
    U1sq = U1 * U1 % p
    # f = U1^2 + a U1 + c, shape (p, ncur) built in chunks
    S1 = np.zeros(ncur, dtype=np.int64)
    for s in range(0, ncur, chunk):
        e = min(s + chunk, ncur)
        fv = (U1sq[:, None] + U1[:, None] * avals[None, s:e] + cvals[None, s:e]) % p
        S1[s:e] = chi[fv].sum(axis=0)

    # --- S2: sum over F_p^2 = F_p[t]/(t^2 - r) --------------------------
    # X = x1 + x2 t;  X^3 = (x1^3 + 3 r x1 x2^2) + (3 x1^2 x2 + r x2^3) t
    x1 = np.repeat(np.arange(p, dtype=np.int64), p)
    x2 = np.tile(np.arange(p, dtype=np.int64), p)
    Ua = (x1 * x1 % p * x1 + 3 * r % p * x1 % p * (x2 * x2 % p)) % p
    Ub = (3 * (x1 * x1 % p) % p * x2 + r * (x2 * x2 % p) % p * x2) % p
    # U^2 = (Ua^2 + r Ub^2) + (2 Ua Ub) t
    Usq_a = (Ua * Ua + r * (Ub * Ub % p)) % p
    Usq_b = (2 * Ua % p) * Ub % p

    S2 = np.zeros(ncur, dtype=np.int64)
    for s in range(0, ncur, chunk):
        e = min(s + chunk, ncur)
        Ra = (Usq_a[:, None] + Ua[:, None] * avals[None, s:e] + cvals[None, s:e]) % p
        Rb = (Usq_b[:, None] + Ub[:, None] * avals[None, s:e]) % p
        # chi_{F_p^2}(Ra + Rb t) = chi_{F_p}(norm) = chi_p(Ra^2 - r Rb^2)
        nrm = (Ra * Ra - r * (Rb * Rb % p)) % p
        S2[s:e] = chi[nrm].sum(axis=0)

    c1 = -1 - S1
    c2 = (c1 * c1 + 1 + S2) // 2
    return avals, cvals, c1, c2, ok


# --------------------------------------------------------------------------
# PARI cross-check
# --------------------------------------------------------------------------


def pari_charpoly(p, a, c, lc=1):
    """(c1, c2) from PARI hyperellcharpoly, for cross-validation."""
    script = (
        f"cp = hyperellcharpoly(Mod(1,{p})*({lc}*(x^6+{a}*x^3+{c})));"
        f'print(polcoef(cp,3)," ",polcoef(cp,2));'
    )
    out = subprocess.run(
        ["gp", "-q"], input=script, capture_output=True, text=True, timeout=300
    )
    if out.returncode != 0:
        raise RuntimeError(out.stderr)
    # gp prints stack-resize warnings on stdout; the answer is the last line.
    lines = [ln for ln in out.stdout.strip().splitlines() if ln.strip()]
    a3, a2 = lines[-1].split()
    return -int(a3), int(a2)  # P = T^4 - c1 T^3 + c2 T^2 - ...


def selftest(p=31, nsample=12):
    print(f"--- selftest: numpy point-count vs PARI hyperellcharpoly, p={p} ---")
    avals, cvals, c1, c2, ok = c3_family_weil(p)
    idx = np.nonzero(ok)[0]
    rng = np.random.default_rng(20260729)
    pick = rng.choice(idx, size=min(nsample, idx.size), replace=False)
    r = non_residue(p)
    bad = 0
    for i in pick:
        a, c = int(avals[i]), int(cvals[i])
        got = (int(c1[i]), int(c2[i]))
        want = pari_charpoly(p, a, c, lc=1)
        twist_got = (-int(c1[i]), int(c2[i]))
        twist_want = pari_charpoly(p, a, c, lc=r)
        s1 = "OK " if got == want else "BAD"
        s2 = "OK " if twist_got == twist_want else "BAD"
        bad += (got != want) + (twist_got != twist_want)
        print(
            f"  a={a:3d} c={c:3d}  lc=1  ours={got} pari={want} {s1}"
            f"   |  lc={r}  ours={twist_got} pari={twist_want} {s2}"
        )
    print(f"--- selftest {'PASSED' if bad == 0 else f'FAILED ({bad} mismatches)'} ---\n")
    return bad == 0


# --------------------------------------------------------------------------
# per-prime driver
# --------------------------------------------------------------------------


def run_prime(p, verbose=True):
    tw = sextic_twists(p)
    pairs = howe_pairs(tw, p)
    avals, cvals, c1, c2, ok = c3_family_weil(p)

    # both quadratic twists of every family member
    fam_c1 = np.concatenate([c1[ok], -c1[ok]])
    fam_c2 = np.concatenate([c2[ok], c2[ok]])
    fam_a = np.concatenate([avals[ok], avals[ok]])
    fam_c = np.concatenate([cvals[ok], cvals[ok]])
    fam_lc = np.concatenate(
        [np.ones(int(ok.sum()), dtype=np.int64), np.full(int(ok.sum()), non_residue(p))]
    )

    if verbose:
        print(f"===== p = {p}   (g = {primitive_root(p)}, "
              f"non-residue = {non_residue(p)}) =====")
        print("  sextic twists:")
        for e in tw:
            typ = "II(split)" if e["tors2"] == (1, 1, 1) else "I (irred)"
            print(f"    k={e['k']}  b={e['b']:4d}  t={e['t']:5d}  "
                  f"#E={e['n']:5d}  2-tors {typ}")
        print(f"  family size (nonsingular, both twists) = {fam_c1.size}")

    results = []
    for row in pairs:
        hit = np.nonzero((fam_c1 == row["c1"]) & (fam_c2 == row["c2"]))[0]
        rec = dict(row)
        rec["nhit"] = int(hit.size)
        rec["example"] = (
            None
            if hit.size == 0
            else (int(fam_lc[hit[0]]), int(fam_a[hit[0]]), int(fam_c[hit[0]]))
        )
        results.append(rec)

    if verbose:
        print("  pair | H1 H2 H3 | Howe | target (c1,c2)      | C3-hits | example "
              "lc,a,c")
        for r_ in results:
            i, j = r_["pair"]
            mark = "YES " if r_["ok"] else "no  "
            ex = "" if r_["example"] is None else str(r_["example"])
            print(
                f"  ({i},{j}) |  {int(r_['H1'])}  {int(r_['H2'])}  {int(r_['H3'])} "
                f"| {mark} | ({r_['c1']:6d},{r_['c2']:8d}) | {r_['nhit']:7d} | {ex}"
            )
        print()
    return tw, results


DEFAULT_PRIMES = [7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103, 109, 127]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--primes", default=None,
                    help="comma-separated primes = 1 mod 6")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        if not selftest():
            sys.exit(1)

    primes = (
        DEFAULT_PRIMES
        if args.primes is None
        else [int(s) for s in args.primes.split(",")]
    )

    summary = []
    for p in primes:
        if p % 6 != 1:
            print(f"skip p={p}: not 1 mod 6 (no 6 distinct sextic twists)")
            continue
        tw, res = run_prime(p)
        q = [r for r in res if r["ok"]]
        summary.append(
            {
                "p": p,
                "n_qualifying": len(q),
                "n_realised": sum(1 for r in q if r["nhit"] > 0),
                "detail": [(r["pair"], r["nhit"]) for r in q],
            }
        )

    print("================ SUMMARY ================")
    print(" p    | Howe-qualifying pairs | realised in C3 family | detail")
    for s in summary:
        det = " ".join(f"{pr}:{n}" for pr, n in s["detail"])
        print(f" {s['p']:4d} | {s['n_qualifying']:21d} | {s['n_realised']:21d} | {det}")


if __name__ == "__main__":
    main()
