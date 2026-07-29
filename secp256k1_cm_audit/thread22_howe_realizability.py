#!/usr/bin/env python3
"""
Thread 22: Is the Howe condition (H1)+(H2)+(H3) actually PREDICTIVE of
(2,2)-glueability over the base field F_p?

Background
----------
Thread 18 (2026-07-21, commit a287abc) applied the codebase's Howe test to the
6 sextic twists of secp256k1 and found 5/15 pairs satisfy (H1)+(H2)+(H3), as
implemented in secp256k1_cm_audit/sextic_twist_howe_check.gp:

    H1 : n_i != n_j                        (distinct group orders)
    H2 : x^3 + b*g^i and x^3 + b*g^j have the same factorisation type over F_p
    H3 : gcd(n_i, n_j) = 1

Thread 2 (2026-07-26, commit c591765) showed the *explicit* Rosenhain
construction is blocked for secp256k1 by a cubic-residue obstruction, and
proposed Thread 22: search directly for a genus-2 curve C/F_p with

    charpoly(Frob | Jac C) = P_{E_i}(T) * P_{E_j}(T).

This script runs that search EXHAUSTIVELY over small proxy primes, where a
complete enumeration of genus-2 curves is feasible, and asks:

    Does (H1)+(H2)+(H3) predict realisability?

Falsifier (stated in advance):
  * A pair PASSES H1+H2+H3 but NO genus-2 curve over F_p realises P_i*P_j
    -> the conditions are NOT sufficient (false positive).
  * A pair FAILS H1+H2+H3 but a genus-2 curve DOES realise P_i*P_j
    -> the conditions are NOT necessary (false negative), i.e. Thread 18's
    "5/15" would be an undercount of the true cover-attack surface.

Method
------
1.  Proxy prime p = 1 (mod 3), fix b, build the six sextic twists
    E_k : y^2 = x^3 + b*g^k, g a primitive root mod p.  (g^k rather than the
    order-6 element h^k used in sextic_twist_howe_check.gp: g^k is
    unconditionally a transversal of F_p^* / (F_p^*)^6.)
    Traces t_k by direct point count; n_k = p + 1 - t_k.
2.  Evaluate H1/H2/H3 for all 15 pairs.
3.  EXHAUSTIVELY enumerate genus-2 curves y^2 = f(x) over F_p:
      deg 5 : f = x^5 + a3 x^3 + a2 x^2 + a1 x + a0            (p^4 models)
      deg 6 : f = c6 x^6 + c4 x^4 + c3 x^3 + c2 x^2 + c1 x + c0,
              c6 in {1, nu}                                     (2 p^5 models)
    Together these cover every F_p-isomorphism class of genus-2 curve: a
    rational Weierstrass point can be moved to infinity (deg-5 case), the
    leading coefficient normalised (x -> a x with y -> c y realises every
    scaling in degree 5 since <5th powers, squares> = F_p^*, and every square
    class in degree 6), and the subleading coefficient killed by a translation
    (needs p > 3).
4.  For each model compute
        s1 = p + 1 - #C(F_p)          via a sum of Legendre symbols
        s2  from #C(F_p^2)            via a sum of quadratic characters of F_p^2
    so that charpoly = T^4 - s1 T^3 + s2 T^2 - p s1 T + p^2.
    Stage 1 (s1) runs on the full grid; stage 2 (s2) only on s1-survivors.
5.  Pair (i,j) is REALISED iff some SQUAREFREE f attains
        s1 = t_i + t_j ,  s2 = t_i t_j + 2p .
    Witnesses are squarefree-checked via deg gcd(f, f') and printed so that
    thread22_pari_crosscheck.gp can re-derive their charpoly independently.

Output: secp256k1_cm_audit/thread22_howe_realizability_output.txt
"""

import sys
import numpy as np
from itertools import combinations

MAX_WITNESS_SCAN = 400          # models inspected per class looking for squarefree

# ---------------------------------------------------------------- F_p basics


def prime_factors(n):
    fs, d = set(), 2
    while d * d <= n:
        while n % d == 0:
            fs.add(d)
            n //= d
        d += 1
    if n > 1:
        fs.add(n)
    return fs


def legendre_table(p):
    """chi[a] = Legendre symbol (a|p), chi[0] = 0."""
    chi = np.zeros(p, dtype=np.int8)
    a = np.arange(1, p, dtype=np.int64)
    chi[(a * a) % p] = 1
    chi[chi == 0] = -1
    chi[0] = 0
    return chi


def primitive_root(p):
    fac = prime_factors(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g
    raise RuntimeError("no primitive root")


def non_residue(p, chi):
    for a in range(2, p):
        if chi[a] == -1:
            return a
    raise RuntimeError("no non-residue")


# ------------------------------------------- F_p^2 = F_p[u]/(u^2 - nu)


def chi2_table(p, nu):
    """chi2[w0, w1] = quadratic character of w0 + w1*u in F_{p^2}."""
    t = np.zeros((p, p), dtype=np.int8)
    z0 = np.arange(p, dtype=np.int64).reshape(p, 1)
    z1 = np.arange(p, dtype=np.int64).reshape(1, p)
    s0 = (z0 * z0 + nu * z1 * z1) % p          # (z0 + z1 u)^2
    s1 = (2 * z0 * z1) % p
    t[s0.ravel(), s1.ravel()] = 1
    t[t == 0] = -1
    t[0, 0] = 0
    return t


def fp2_powers(p, nu, z0, z1, kmax):
    """[(z^k)_0, (z^k)_1] for k = 0..kmax where z = z0 + z1*u."""
    out = [(1, 0)]
    c0, c1 = 1, 0
    for _ in range(kmax):
        c0, c1 = (c0 * z0 + nu * c1 * z1) % p, (c0 * z1 + c1 * z0) % p
        out.append((c0, c1))
    return out


# --------------------------------------------------------- elliptic-curve side


PATTERN = {0: "[3] irred", 1: "[1,2]", 3: "[1,1,1] split"}


def ec_trace(p, b, chi):
    x = np.arange(p, dtype=np.int64)
    return int(-chi[(x * x % p * x + b) % p].sum())


def cubic_root_count(p, c):
    x = np.arange(p, dtype=np.int64)
    return int((((x * x % p) * x + c) % p == 0).sum())


def sextic_twist_data(p, b, chi):
    g = primitive_root(p)
    ts, ns, pats, bs = [], [], [], []
    gk = 1
    for _ in range(6):
        bk = (b * gk) % p
        bs.append(bk)
        ts.append(ec_trace(p, bk, chi))
        ns.append(p + 1 - ts[-1])
        pats.append(cubic_root_count(p, bk))
        gk = (gk * g) % p
    return g, bs, ts, ns, pats


def howe_table(p, ts, ns, pats):
    rows = []
    for i, j in combinations(range(6), 2):
        h1 = ns[i] != ns[j]
        h2 = pats[i] == pats[j]
        gc = int(np.gcd(ns[i], ns[j]))
        rows.append(dict(pair=(i, j), H1=bool(h1), H2=bool(h2), H3=(gc == 1),
                         ALL=bool(h1 and h2 and gc == 1), gcd=gc,
                         s1=ts[i] + ts[j], s2=ts[i] * ts[j] + 2 * p))
    return rows


# ------------------------------------------------- genus-2 exhaustive search


def _grid(p, ndim):
    out = []
    for d in range(ndim):
        shape = [1] * ndim
        shape[d] = p
        out.append(np.arange(p, dtype=np.int64).reshape(shape))
    return out


def _stage1_deg5(p, chi):
    A3, A2, A1, A0 = _grid(p, 4)
    S = np.zeros((p, p, p, p), dtype=np.int16)
    for x in range(p):
        x2 = x * x % p
        x3 = x2 * x % p
        x5 = x3 * x2 % p
        S += chi[(x5 + A3 * x3 + A2 * x2 + A1 * x + A0) % p]
    return (-S.astype(np.int64))                      # s1


def _stage1_deg6(p, chi, c6):
    C4, C3, C2, C1, C0 = _grid(p, 5)
    S = np.zeros((p,) * 5, dtype=np.int16)
    for x in range(p):
        x2 = x * x % p
        x3 = x2 * x % p
        x4 = x2 * x2 % p
        x6 = x3 * x3 % p
        S += chi[(c6 * x6 + C4 * x4 + C3 * x3 + C2 * x2 + C1 * x + C0) % p]
    return -(S.astype(np.int64) + int(chi[c6]))       # s1


def _stage2(p, chi2, nu, coeffs, s1v, deg, c6):
    """s2 for the surviving models, from #C(F_{p^2}).

    Conjugate field elements contribute equally (f has F_p coefficients and
    chi2(w^p) = chi2(w)), so only z1 <= (p-1)/2 is swept, doubling z1 > 0.
    """
    n = s1v.shape[0]
    S2 = np.zeros(n, dtype=np.int32)
    coeffs = [c.astype(np.int32) for c in coeffs]
    half = (p - 1) // 2
    for z1 in range(0, half + 1):
        mult = 1 if z1 == 0 else 2
        for z0 in range(p):
            pw = fp2_powers(p, nu, z0, z1, 6)
            if deg == 5:
                a3, a2, a1, a0 = coeffs
                f0 = (pw[5][0] + a3 * pw[3][0] + a2 * pw[2][0]
                      + a1 * pw[1][0] + a0) % p
                f1 = (pw[5][1] + a3 * pw[3][1] + a2 * pw[2][1]
                      + a1 * pw[1][1]) % p
            else:
                c4, c3, c2, c1, c0 = coeffs
                f0 = (c6 * pw[6][0] + c4 * pw[4][0] + c3 * pw[3][0]
                      + c2 * pw[2][0] + c1 * pw[1][0] + c0) % p
                f1 = (c6 * pw[6][1] + c4 * pw[4][1] + c3 * pw[3][1]
                      + c2 * pw[2][1] + c1 * pw[1][1]) % p
            if mult == 1:
                S2 += chi2[f0, f1]
            else:
                S2 += 2 * chi2[f0, f1].astype(np.int32)
    extra = 0 if deg == 5 else 1          # chi2(c6) = 1 for c6 in F_p^*
    T2 = -(S2.astype(np.int64) + extra)   # = s1^2 - 2 s2
    num = s1v * s1v - T2
    assert not (num & 1).any(), "s1^2 - T2 must be even"
    return num // 2


def sweep(p, chi, chi2, nu, target_s1, deg, c6=None):
    """Return (class -> witness-index-array, coefficient arrays) for models
    whose s1 lies in target_s1."""
    s1 = _stage1_deg5(p, chi) if deg == 5 else _stage1_deg6(p, chi, c6)
    mask = np.isin(s1, target_s1)
    if not mask.any():
        return {}, None, None
    idx = np.nonzero(mask)
    coeffs = [a.astype(np.int64) for a in idx]   # index == coefficient value
    s1v = s1[mask].astype(np.int64)
    del s1, mask
    s2v = _stage2(p, chi2, nu, coeffs, s1v, deg, c6)
    return coeffs, s1v, s2v


def witness_poly(deg, c6, cs):
    """coefficient list, lowest degree first."""
    if deg == 5:
        a3, a2, a1, a0 = cs
        return [a0, a1, a2, a3, 0, 1]
    c4, c3, c2, c1, c0 = cs
    return [c0, c1, c2, c3, c4, 0, c6]


# ---------------------------------------------------- squarefree witness check


def poly_gcd_deg(p, f):
    """deg gcd(f, f') over F_p; 0 means squarefree, -1 means gcd is 0."""
    def trim(a):
        while a and a[-1] % p == 0:
            a.pop()
        return a

    def rem(a, b):
        a = a[:]
        inv = pow(b[-1], p - 2, p)
        while len(a) >= len(b) and a:
            c = a[-1] * inv % p
            sh = len(a) - len(b)
            for i, bc in enumerate(b):
                a[i + sh] = (a[i + sh] - c * bc) % p
            trim(a)
        return a

    f = trim([c % p for c in f])
    d = trim([(i * c) % p for i, c in enumerate(f)][1:])
    a, b = f, d
    while b:
        a, b = b, rem(a, b)
    return len(a) - 1 if a else -1


def find_squarefree(p, deg, c6, coeffs, sel):
    """First squarefree model among the selected indices (capped scan)."""
    scanned = 0
    for k in sel[:MAX_WITNESS_SCAN]:
        cs = tuple(int(a[k]) for a in coeffs)
        poly = witness_poly(deg, c6, cs)
        scanned += 1
        if poly_gcd_deg(p, poly) == 0:
            return cs, scanned
    return None, scanned


# ------------------------------------------------------------------ driver


def run_prime(p, b, out, full=True):
    chi = legendre_table(p)
    nu = non_residue(p, chi)
    chi2 = chi2_table(p, nu)
    g, bs, ts, ns, pats = sextic_twist_data(p, b, chi)
    rows = howe_table(p, ts, ns, pats)
    target_s1 = sorted({r["s1"] for r in rows})

    families = [(5, None)] + ([(6, 1), (6, nu)] if full else [])
    results = []
    for deg, c6 in families:
        print(f"    sweep deg={deg} c6={c6}", file=sys.stderr)
        results.append((deg, c6) + sweep(p, chi, chi2, nu, target_s1, deg, c6))

    # realisability per pair
    realised = {}
    for r in rows:
        key = (r["s1"], r["s2"])
        if key in realised:
            continue
        hit = None
        nmatch = 0
        for deg, c6, coeffs, s1v, s2v in results:
            if coeffs is None:
                continue
            sel = np.nonzero((s1v == r["s1"]) & (s2v == r["s2"]))[0]
            nmatch += int(sel.shape[0])
            if hit is None and sel.shape[0]:
                cs, _ = find_squarefree(p, deg, c6, coeffs, sel)
                if cs is not None:
                    hit = (deg, c6, cs)
        realised[key] = (hit, nmatch)

    w = out.write
    w("=" * 78 + "\n")
    w(f"p = {p}   b = {b}   primitive root g = {g}   non-residue nu = {nu}\n")
    w(f"p mod 3 = {p % 3},  p mod 4 = {p % 4},  p mod 12 = {p % 12}\n")
    w("=" * 78 + "\n")
    w("sextic twists  E_k : y^2 = x^3 + b*g^k\n")
    w(f"{'k':>2} {'b*g^k':>7} {'t_k':>6} {'n_k':>6}  2-torsion pattern\n")
    for k in range(6):
        w(f"{k:>2} {bs[k]:>7} {ts[k]:>6} {ns[k]:>6}  {PATTERN[pats[k]]}\n")
    w(f"trace sums: t0+t2+t4 = {ts[0]+ts[2]+ts[4]}, "
      f"t1+t3+t5 = {ts[1]+ts[3]+ts[5]}\n")
    cover = ("COMPLETE (deg 5 + deg 6, all F_p-iso classes)" if full
             else "PARTIAL (deg 5 only; classes realised only by curves with "
                  "no rational Weierstrass point may be missed)")
    w(f"enumeration coverage: {cover}\n\n")

    w("Howe conditions vs. exhaustive genus-2 realisability:\n")
    w(f"{'pair':>5} {'H1':>2} {'H2':>2} {'H3':>2} {'ALL':>4} {'gcd':>4} "
      f"{'s1':>4} {'s2':>5} {'#models':>8} {'real':>5}  verdict\n")
    fp = fn = 0
    for r in rows:
        hit, nmatch = realised[(r["s1"], r["s2"])]
        real = hit is not None
        verdict = ""
        if r["ALL"] and not real:
            verdict = "FALSE POSITIVE (H passes, not realised)"
            fp += 1
        elif (not r["ALL"]) and real:
            verdict = "FALSE NEGATIVE (H fails, realised anyway)"
            fn += 1
        i, j = r["pair"]
        w(f"({i},{j}) {'Y' if r['H1'] else 'N':>2} {'Y' if r['H2'] else 'N':>2} "
          f"{'Y' if r['H3'] else 'N':>2} {'YES' if r['ALL'] else 'no':>4} "
          f"{r['gcd']:>4} {r['s1']:>4} {r['s2']:>5} {nmatch:>8} "
          f"{'YES' if real else 'no':>5}  {verdict}\n")

    nq = sum(1 for r in rows if r["ALL"])
    nr = sum(1 for r in rows if realised[(r["s1"], r["s2"])][0] is not None)
    w(f"\nHowe-qualifying : {nq}/15    realised : {nr}/15    "
      f"false-pos : {fp}    false-neg : {fn}\n")

    w("\nsquarefree witnesses (y^2 = f(x)); feed to hyperellcharpoly:\n")
    seen = set()
    gp_rows = []
    for r in rows:
        key = (r["s1"], r["s2"])
        hit, _ = realised[key]
        if hit is None or key in seen:
            continue
        seen.add(key)
        deg, c6, cs = hit
        poly = witness_poly(deg, c6, cs)
        terms = " + ".join(f"{c}*x^{k}" for k, c in enumerate(poly) if c)
        w(f"  pair {r['pair']} s1={r['s1']:>4} s2={r['s2']:>5}  "
          f"p={p}  f = {terms}\n")
        gp_rows.append((p, r["pair"][0], r["pair"][1], r["s1"], r["s2"], poly))
    w("\n")
    out.flush()
    return dict(p=p, b=b, qualifying=nq, realised=nr, fp=fp, fn=fn,
                full=full, rows=rows, realmap=realised, gp_rows=gp_rows)


def main():
    # p = 1 mod 3 so that six distinct sextic twists exist.
    # full=True -> complete enumeration (deg 5 AND deg 6).
    cases = [(7, 1, True), (13, 7, True), (19, 7, True), (31, 7, False),
             (37, 7, False), (43, 7, False)]
    path = "secp256k1_cm_audit/thread22_howe_realizability_output.txt"
    summ = []
    with open(path, "w") as out:
        out.write("Thread 22: Howe (H1)+(H2)+(H3) vs. exhaustive genus-2 "
                  "realisability over small proxy primes\n\n")
        for p, b, full in cases:
            print(f"[thread22] p={p} b={b} full={full}", file=sys.stderr)
            summ.append(run_prime(p, b, out, full=full))
        out.write("=" * 78 + "\nSUMMARY\n" + "=" * 78 + "\n")
        out.write(f"{'p':>5} {'b':>3} {'coverage':>9} {'H-qual':>7} "
                  f"{'realised':>9} {'falsePos':>9} {'falseNeg':>9}\n")
        for s in summ:
            out.write(f"{s['p']:>5} {s['b']:>3} "
                      f"{'full' if s['full'] else 'deg5':>9} "
                      f"{s['qualifying']:>7} {s['realised']:>9} "
                      f"{s['fp']:>9} {s['fn']:>9}\n")
        out.write(f"\ntotal false positives (H passes, no genus-2 curve): "
                  f"{sum(s['fp'] for s in summ)}\n")
        out.write(f"total false negatives (H fails, curve exists)     : "
                  f"{sum(s['fn'] for s in summ)}\n")
        out.write("\nNOTE: false negatives on deg5-only rows are sound "
                  "(a witness was found);\n"
                  "      false positives on deg5-only rows are NOT "
                  "conclusive (coverage is partial).\n")

    # witness table for the independent PARI cross-check
    with open("secp256k1_cm_audit/thread22_witnesses.gp", "w") as gp:
        gp.write("\\\\ auto-generated by thread22_howe_realizability.py\n")
        gp.write("\\\\ [p, i, j, s1, s2, f coefficients low->high]\n")
        gp.write("{W = [\n")
        items = [row for s in summ for row in s["gp_rows"]]
        for k, (p, i, j, s1, s2, poly) in enumerate(items):
            sep = "," if k + 1 < len(items) else ""
            gp.write(f"  [{p}, {i}, {j}, {s1}, {s2}, "
                     f"[{', '.join(str(c) for c in poly)}]]{sep}\n")
        gp.write("];}\n")
    print(open(path).read())


if __name__ == "__main__":
    main()
