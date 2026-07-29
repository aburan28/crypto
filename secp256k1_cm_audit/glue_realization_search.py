"""
Thread 22: Exhaustive glue-realization search over small proxy primes.

QUESTION
--------
Thread 18 (2026-07-21) found that 5 of the 15 pairs of j=0 sextic twists of
secp256k1 satisfy Howe's conditions (H1) non-isogenous, (H2) matching 2-torsion
Galois structure, (H3) gcd of orders = 1.  Howe (1996) then guarantees that
A = (E_i x E_j)/Gamma_alpha is the Jacobian of a smooth genus-2 curve C/F_p.

Thread 2 (2026-07-26) showed the Rosenhain/CHLRS route to writing C down
explicitly is BLOCKED over F_p for secp256k1 (cubic-residue obstruction:
x^3 + 7 has no root in F_p, so E[2] is not F_p-rational).

This script answers the *realization* question empirically at proxy scale:

    For a small prime p = 1 mod 6, does there EXIST a genus-2 curve C/F_p with
        charpoly_Frob(Jac C) = P_i(T) * P_j(T),
        P_k(T) = T^2 - t_k T + p,
    for each of the 15 pairs (i,j) of j=0 sextic twists?

We settle it by EXHAUSTIVE enumeration of all genus-2 curves over F_p (up to
F_p-isomorphism), so a negative answer is a proof of non-realization at that p,
not a search failure.

NOTE ON THREAD 7.  Thread 7 (2026-07-09) closed the depth-2 Richelot question
via Honda-Tate: P_i*P_j is reducible, so every abelian surface in that isogeny
class is non-simple.  That does NOT answer the present question: a Jacobian of
a *smooth* genus-2 curve can perfectly well be isogenous (not isomorphic) to a
product of elliptic curves -- that is exactly what "gluing" produces.  So the
realization question is independent of Thread 7 and is settled here.

COMPLETENESS OF THE ENUMERATION
-------------------------------
Every genus-2 curve over F_p (p > 5) has a model y^2 = f(x) with deg f in {5,6}.

  Case A (deg 5).  y^2 = c*g(x), g monic quintic.  The substitution
  x = c*X, y = c^3*Y turns it into Y^2 = X^5 + (b4/c)X^4 + (b3/c^2)X^3 + ...,
  i.e. every odd-degree model is isomorphic to a MONIC one (in particular the
  quadratic twist is already inside the monic family -- no separate twist
  needed).  Then x -> x - b4/5 kills the x^4 term (p != 5).
      => family A = { x^5 + a3 x^3 + a2 x^2 + a1 x + a0 },  p^4 models.

  Case B (deg 6).  y^2 = c*g(x), g monic sextic.  x -> u*x, y -> u^3*y leaves
  the leading coefficient c fixed, so c matters exactly modulo squares:
  c in {1, nu} with nu a fixed non-residue.  Then x -> x - b5/6 kills the x^5
  term (p != 2,3).
      => family B = { c*(x^6 + b4 x^4 + b3 x^3 + b2 x^2 + b1 x + b0) },
         2*p^5 models.

Family A u family B hits every F_p-isomorphism class of genus-2 curve at least
once (with redundancy, which is harmless).

Family B is the essential one here: for the secp256k1-like twists x^3+b is
irreducible, so the glued curve's Weierstrass points fall in Galois orbits of
size 3 or 6 and C has NO rational Weierstrass point -- family A would miss it.

ALGORITHM
---------
1.  Enumerate families A and B, computing N1 = #C(F_p) with a Legendre table
    (numpy-vectorised over the two lowest coefficients).
2.  Keep the curves whose s1 = p+1-N1 equals t_i+t_j for some pair (i,j).
3.  For those survivors compute N2 = #C(F_{p^2}) and hence the full Frobenius
    charpoly  T^4 - e1 T^3 + e2 T^2 - p e1 T + p^2  with
        e1 = s1,  e2 = (s1^2 - s2)/2,  s2 = p^2 + 1 - N2.
4.  Discard singular models (gcd(f, f') != 1).
5.  Report, per pair, whether P_i*P_j is realized and by which curve.

Run:  python3 glue_realization_search.py [p ...]      (default: 7 13 19)
"""

import sys
import numpy as np

# ---------------------------------------------------------------------------
# F_p helpers
# ---------------------------------------------------------------------------


def legendre_table(p):
    """chi[a] = 0, 1, -1 for a = 0, QR, non-QR."""
    chi = -np.ones(p, dtype=np.int64)
    chi[0] = 0
    for a in range(1, p):
        chi[(a * a) % p] = 1
    return chi


def non_residue(p, chi):
    for a in range(2, p):
        if chi[a] == -1:
            return a
    raise RuntimeError("no non-residue")


def poly_gcd_deg(f, p):
    """Degree of gcd(f, f') over F_p; f is a coefficient list, low -> high."""
    fp = [(i * c) % p for i, c in enumerate(f)][1:]

    def norm(a):
        while a and a[-1] % p == 0:
            a.pop()
        return a

    a, b = norm(list(f)), norm(list(fp))
    while b:
        inv = pow(b[-1], -1, p)
        while len(a) >= len(b):
            fac = (a[-1] * inv) % p
            shift = len(a) - len(b)
            for i, c in enumerate(b):
                a[i + shift] = (a[i + shift] - fac * c) % p
            a = norm(a)
            if not a:
                break
        a, b = b, a
    return len(a) - 1 if a else -1


# ---------------------------------------------------------------------------
# j = 0 sextic twists
# ---------------------------------------------------------------------------


def primitive_root(p):
    fac = []
    m, d = p - 1, 2
    while d * d <= m:
        if m % d == 0:
            fac.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        fac.append(m)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g
    raise RuntimeError("no primitive root")


def sextic_twists(p, chi):
    """Return (bs, traces, root_patterns) for the 6 sextic twists y^2 = x^3 + b."""
    g = primitive_root(p)
    bs = [pow(g, k, p) for k in range(6)]
    traces, pats = [], []
    xs = np.arange(p, dtype=np.int64)
    cubes = (xs * xs % p) * xs % p
    for b in bs:
        vals = (cubes + b) % p
        # affine points + 1 point at infinity
        n_aff = int(p + chi[vals].sum())
        traces.append(p + 1 - (n_aff + 1))
        pats.append(int(np.count_nonzero((cubes + b) % p == 0)))
    return bs, traces, pats


# ---------------------------------------------------------------------------
# F_{p^2} arithmetic, index form  e = a + b*p  with  u^2 = nu
# ---------------------------------------------------------------------------


class Fp2:
    def __init__(self, p, nu):
        self.p, self.nu = p, nu
        q = p * p
        idx = np.arange(q, dtype=np.int64)
        self.A = idx % p
        self.B = idx // p
        # square-class table over F_{p^2}
        sa = (self.A * self.A + nu * self.B * self.B) % p
        sb = (2 * self.A * self.B) % p
        sq = np.zeros(q, dtype=np.int64)
        sq[sa + p * sb] = 1
        sq[0] = 0
        self.chi2 = np.where(sq == 1, 1, -1)
        self.chi2[0] = 0

    def mul(self, a1, b1, a2, b2):
        p, nu = self.p, self.nu
        return ((a1 * a2 + nu * b1 * b2) % p, (a1 * b2 + a2 * b1) % p)

    def horner(self, coeffs_hi_to_lo):
        """Evaluate a batch of polynomials at every point of F_{p^2}.

        coeffs_hi_to_lo: (nbatch, ncoef) int array of F_p coefficients,
        highest degree first.  Returns (ra, rb) of shape (nbatch, p^2).
        """
        p = self.p
        n = coeffs_hi_to_lo.shape[0]
        xa = self.A[None, :]
        xb = self.B[None, :]
        ra = np.zeros((n, p * p), dtype=np.int64)
        rb = np.zeros((n, p * p), dtype=np.int64)
        for k in range(coeffs_hi_to_lo.shape[1]):
            ra, rb = self.mul(ra, rb, xa, xb)
            ra = (ra + coeffs_hi_to_lo[:, k, None]) % p
        return ra, rb


# ---------------------------------------------------------------------------
# charpoly of a batch of survivors
# ---------------------------------------------------------------------------


def charpolys(p, f2, coeffs, deg, lead, s1):
    """coeffs: (n, deg+1) hi->lo, monic part; lead: (n,) leading coeff in F_p.

    Returns array of (e1, e2) for each curve.
    """
    ra, rb = f2.horner(coeffs)
    la = lead[:, None] % p
    ra, rb = f2.mul(ra, rb, la, np.zeros_like(la))
    chi = f2.chi2[ra + p * rb]
    n_aff = p * p + chi.sum(axis=1)
    # points at infinity over F_{p^2}: every element of F_p* is a square in
    # F_{p^2}, so a sextic model always has 2, a quintic model always has 1.
    n_inf = 2 if deg == 6 else 1
    n2 = n_aff + n_inf
    s2 = p * p + 1 - n2
    e1 = s1
    e2 = (s1 * s1 - s2) // 2
    return e1, e2


# ---------------------------------------------------------------------------
# main sweep
# ---------------------------------------------------------------------------


def sweep(p, verbose=True, early_exit=False, only_diff1=False):
    assert p % 6 == 1, "need p = 1 mod 6 for six distinct sextic twists"
    chi = legendre_table(p)
    nu = non_residue(p, chi)
    bs, traces, pats = sextic_twists(p, chi)
    orders = [p + 1 - t for t in traces]

    print("=" * 72)
    print(f"p = {p}   (p mod 6 = {p % 6}),  non-residue nu = {nu}")
    print("=" * 72)
    print("j=0 sextic twists  y^2 = x^3 + b:")
    for k in range(6):
        print(
            f"  k={k}: b={bs[k]:3d}  t={traces[k]:4d}  n={orders[k]:4d}  "
            f"#roots(x^3+b)={pats[k]}"
        )
    # the 6 traces are +-(2a-b), +-(a-2b), +-(a+b) with the three squares
    # summing to 6p, so the sum over all six twists is 12p.
    ss = sum(t * t for t in traces)
    print(f"  sum t_k = {sum(traces)}   sum t_k^2 = {ss}  (= 12p? {ss == 12 * p})")

    # ---- Howe conditions per pair -------------------------------------
    pairs = []
    for i in range(6):
        for j in range(i + 1, 6):
            h1 = traces[i] != traces[j]
            h2 = pats[i] == pats[j]
            from math import gcd

            h3 = gcd(orders[i], orders[j]) == 1
            pairs.append(
                {
                    "i": i,
                    "j": j,
                    "H1": h1,
                    "H2": h2,
                    "H3": h3,
                    "howe": h1 and h2 and h3,
                    "e1": traces[i] + traces[j],
                    "e2": traces[i] * traces[j] + 2 * p,
                }
            )

    # --only-diff1 restricts the (expensive) search to the pairs the
    # |t_i - t_j| = 1 criterion predicts are NOT realizable.  Everything else
    # is reported as "not searched" rather than as a negative.
    for pr in pairs:
        pr["searched"] = (not only_diff1) or abs(
            traces[pr["i"]] - traces[pr["j"]]
        ) == 1
    searched = [pr for pr in pairs if pr["searched"]]

    targets = {}
    for pr in searched:
        targets.setdefault(pr["e1"], set()).add(pr["e2"])
    target_s1 = np.array(sorted(targets), dtype=np.int64)
    if verbose:
        print(f"\ntarget s1 values (= t_i+t_j): {sorted(targets)}")

    # encode the 15 target charpolys (e1, e2) as single integers so the
    # expensive per-curve squarefree test only runs on actual hits.
    SHIFT = 8 * p + 1

    def code(e1, e2):
        return (e1 + 4 * p) * SHIFT + (e2 + 4 * p)

    target_codes = np.array(
        sorted({code(pr["e1"], pr["e2"]) for pr in searched}), dtype=np.int64
    )

    xs = np.arange(p, dtype=np.int64)
    x2 = xs * xs % p
    x3 = x2 * xs % p
    x4 = x3 * xs % p
    x5 = x4 * xs % p
    x6 = x5 * xs % p
    s1_ok = np.zeros(4 * p + 1, dtype=bool)  # index s1 + 2p
    for v in target_s1:
        s1_ok[v + 2 * p] = True

    found = {}  # (e1, e2) -> model description
    found_coefs = {}  # (e1, e2) -> (lead, [coeffs hi->lo]) for PARI re-verification

    def record(e1, e2, desc, lead=None, cf=None):
        key = (int(e1), int(e2))
        if key not in found:
            found[key] = desc
            found_coefs[key] = (int(lead), [int(v) for v in cf])

    f2 = Fp2(p, nu)
    B1, B0 = np.meshgrid(xs, xs, indexing="ij")  # (b1, b0)
    B1f, B0f = B1.ravel(), B0.ravel()

    # S[b1, b0] = sum_x chi(base[x] + b1*x + b0) is the inner-loop hot spot.
    # Writing H[b1, v] = #{x : base[x] + b1*x = v} turns it into the matrix
    # product  S = H @ CHIMAT  with CHIMAT[v, b0] = chi(v + b0), i.e. a BLAS
    # p x p x p matmul instead of a p^3 int64 gather.  Counts and chi values
    # are tiny, so float64 is exact here.
    CHIMAT = chi[(xs[:, None] + xs[None, :]) % p].astype(np.float64)
    ROWOFF = (xs * p)[:, None]  # for the bincount trick
    XB1 = np.outer(xs, xs) % p  # XB1[b1, x] = b1 * x

    def inner_S(base):
        G = (base[None, :] + XB1) % p  # G[b1, x]
        H = np.bincount((G + ROWOFF).ravel(), minlength=p * p).reshape(p, p)
        return np.rint(H.astype(np.float64) @ CHIMAT).astype(np.int64)

    # ---------------- family B: y^2 = c*(x^6 + b4x^4 + b3x^3 + b2x^2 + b1x + b0)
    # Early exit once all 15 target charpolys are hit: exhaustiveness is only
    # needed to *prove* a non-realization, so it costs nothing when all 15 are
    # realized.  `complete` records whether the sweep actually ran to the end.
    ntarget = len({(pr["e1"], pr["e2"]) for pr in searched})
    complete = True
    nsurv = 0
    for b4 in range(p):
        if early_exit and len(found) == ntarget:
            complete = False
            break
        for b3 in range(p):
            if early_exit and len(found) == ntarget:
                complete = False
                break
            for b2 in range(p):
                base = (x6 + b4 * x4 + b3 * x3 + b2 * x2) % p
                S = inner_S(base)  # S[b1, b0]
                for c, sgn, ninf in ((1, 1, 2), (nu, -1, 0)):
                    n1 = p + sgn * S + ninf
                    s1 = p + 1 - n1
                    m = s1_ok[s1 + 2 * p]
                    if not m.any():
                        continue
                    idx = np.flatnonzero(m.ravel())
                    b1v, b0v = B1f[idx], B0f[idx]
                    coeffs = np.stack(
                        [
                            np.ones(len(idx), dtype=np.int64),
                            np.zeros(len(idx), dtype=np.int64),
                            np.full(len(idx), b4, dtype=np.int64),
                            np.full(len(idx), b3, dtype=np.int64),
                            np.full(len(idx), b2, dtype=np.int64),
                            b1v,
                            b0v,
                        ],
                        axis=1,
                    )
                    s1v = (p + 1 - (p + sgn * S + ninf)).ravel()[idx]
                    lead = np.full(len(idx), c, dtype=np.int64)
                    e1v, e2v = charpolys(p, f2, coeffs, 6, lead, s1v)
                    nsurv += len(idx)
                    hits = np.flatnonzero(np.isin(code(e1v, e2v), target_codes))
                    for t in hits:
                        # only now pay for the squarefree (non-singular) test
                        if poly_gcd_deg(
                            [int(b0v[t]), int(b1v[t]), b2, b3, b4, 0, 1], p
                        ):
                            continue
                        record(
                            e1v[t],
                            e2v[t],
                            f"y^2 = {c}*(x^6 + {b4}x^4 + {b3}x^3 + {b2}x^2 "
                            f"+ {b1v[t]}x + {b0v[t]})",
                            lead=c,
                            cf=coeffs[t],
                        )

    # ---------------- family A: y^2 = x^5 + a3x^3 + a2x^2 + a1x + a0
    for a3 in range(p):
        if early_exit and len(found) == ntarget:
            complete = False
            break
        for a2 in range(p):
            base = (x5 + a3 * x3 + a2 * x2) % p
            S = inner_S(base)
            n1 = p + S + 1
            s1 = p + 1 - n1
            m = s1_ok[s1 + 2 * p]
            if not m.any():
                continue
            idx = np.flatnonzero(m.ravel())
            a1v, a0v = B1f[idx], B0f[idx]
            coeffs = np.stack(
                [
                    np.ones(len(idx), dtype=np.int64),
                    np.zeros(len(idx), dtype=np.int64),
                    np.full(len(idx), a3, dtype=np.int64),
                    np.full(len(idx), a2, dtype=np.int64),
                    a1v,
                    a0v,
                ],
                axis=1,
            )
            s1v = s1.ravel()[idx]
            lead = np.ones(len(idx), dtype=np.int64)
            e1v, e2v = charpolys(p, f2, coeffs, 5, lead, s1v)
            nsurv += len(idx)
            hits = np.flatnonzero(np.isin(code(e1v, e2v), target_codes))
            for t in hits:
                if poly_gcd_deg([int(a0v[t]), int(a1v[t]), a2, a3, 0, 1], p):
                    continue
                record(
                    e1v[t],
                    e2v[t],
                    f"y^2 = x^5 + {a3}x^3 + {a2}x^2 + {a1v[t]}x + {a0v[t]}",
                    lead=1,
                    cf=coeffs[t],
                )

    print(f"\nmodels whose s1 = t_i+t_j for some pair (charpoly computed): {nsurv}")
    print(f"distinct target charpolys hit by a smooth model: {len(found)}/{ntarget}")
    print(
        "enumeration: EXHAUSTIVE (negatives are proofs)"
        if complete
        else "enumeration: early-exited after all targets hit (no negatives to prove)"
    )

    # ---- report --------------------------------------------------------
    print("\npair  H1 H2 H3  Howe   t_i   t_j  |ti-tj|  (e1,e2)      realized?")
    nreal_howe = nreal_non = 0
    rows = []
    for pr in pairs:
        key = (pr["e1"], pr["e2"])
        hit = found.get(key) if pr["searched"] else None
        pr["realized"] = hit is not None
        pr["model"] = hit
        if not pr["searched"]:
            pass
        elif pr["howe"]:
            nreal_howe += bool(hit)
        else:
            nreal_non += bool(hit)
        pr["p"] = p
        pr["t_i"] = traces[pr["i"]]
        pr["t_j"] = traces[pr["j"]]
        pr["coefs"] = found_coefs.get(key)
        rows.append(pr)
        dt = abs(traces[pr["i"]] - traces[pr["j"]])
        print(
            f"({pr['i']},{pr['j']})   {'Y' if pr['H1'] else 'N'}  "
            f"{'Y' if pr['H2'] else 'N'}  {'Y' if pr['H3'] else 'N'}   "
            f"{'YES' if pr['howe'] else ' no'}  "
            f"{traces[pr['i']]:4d} {traces[pr['j']]:4d}   "
            f"{dt:4d}   ({pr['e1']:4d},{pr['e2']:5d})   "
            + (
                ("REALIZED" if hit else "NOT REALIZED")
                if pr["searched"]
                else "(not searched)"
            )
        )

    nhowe = sum(1 for pr in searched if pr["howe"])
    nnon = len(searched) - nhowe
    print(f"\n(searched pairs: {len(searched)}/15)")
    print(f"Howe-qualifying pairs: {nhowe},  realized: {nreal_howe}/{nhowe}")
    print(f"non-qualifying pairs : {nnon},  realized: {nreal_non}/{nnon}")

    # the |t_i - t_j| = 1 criterion (see the module docstring / Thread 22 note)
    d1 = [pr for pr in searched if abs(traces[pr["i"]] - traces[pr["j"]]) == 1]
    dg = [pr for pr in searched if abs(traces[pr["i"]] - traces[pr["j"]]) != 1]
    print(
        f"|t_i-t_j| = 1 : {len(d1):2d} pairs, realized "
        f"{sum(1 for pr in d1 if pr['realized'])}"
    )
    print(
        f"|t_i-t_j| > 1 : {len(dg):2d} pairs, realized "
        f"{sum(1 for pr in dg if pr['realized'])}"
    )
    print("\nexample models for realized pairs:")
    for pr in rows:
        if pr["realized"]:
            tag = "HOWE" if pr["howe"] else "non-Howe"
            print(f"  ({pr['i']},{pr['j']}) [{tag}]: {pr['model']}")
    print()
    return rows


if __name__ == "__main__":
    flags = {"--early-exit", "--only-diff1"}
    argv = [a for a in sys.argv[1:] if a not in flags]
    early = "--early-exit" in sys.argv
    d1only = "--only-diff1" in sys.argv
    ps = [int(a) for a in argv] or [7, 13, 19]
    out = []
    for p in ps:
        for pr in sweep(p, early_exit=early, only_diff1=d1only):
            if pr["realized"]:
                lead, cf = pr["coefs"]
                out.append(
                    "%d %d %d %d %d %d %d %d %d %s"
                    % (
                        pr["p"],
                        pr["i"],
                        pr["j"],
                        pr["t_i"],
                        pr["t_j"],
                        pr["e1"],
                        pr["e2"],
                        int(pr["howe"]),
                        lead,
                        ",".join(str(v) for v in cf),
                    )
                )
    if d1only:
        print("(--only-diff1: no models file written)")
        sys.exit(0)
    path = "secp256k1_cm_audit/glue_realization_models.txt"
    with open(path, "w") as fh:
        fh.write("# p i j t_i t_j e1 e2 howe lead coeffs(hi->lo)\n")
        fh.write("\n".join(out) + "\n")
    print(f"wrote {len(out)} models to {path}")
