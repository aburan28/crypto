"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector is
reachable.  Follow-up to 2026-07-29 (Thread 20), which established:

  * the planted vector is NEVER lambda_1 (the trivial vector n*S_D*e_m is
    always shorter), so recovery is a BDD/coset condition, not SVP;
  * the viability wall is a K1 wall (eff = K1*K2/n), not a lambda wall;
  * no curve-level invariant separates success from failure.

This run tests two reformulations against a quantitative predictor.

PREDICTOR (Gaussian heuristic).  With column scalings S_K1 = n/K1,
S_K2 = n/K2, S_D = 1, S_KANNAN = n, and eff = K1*K2/n:

    det(L) = (n*S_K1)^m * S_D * S_K2^m * S_KANNAN = n^(2m+1) * eff^(-m)
    E||v_planted||^2 = m*(K1^2/a)*S_K1^2 + m*(K2^2/a)*S_K2^2
                       + (n^2/3)*S_D^2 + S_KANNAN^2
    gap = GH(det, D) / ||v_planted||          (recovery plausible iff gap > 1)

where a = 3 for uncentered nonce limbs (k in [0,K)) and a = 12 for centered
limbs (u = k - K/2 in [-K/2, K/2)).  Asymptotically in m,

    gap -> sqrt(a) / (sqrt(2*pi*e) * sqrt(eff))

so the ASYMPTOTIC wall sits at eff = a / (2*pi*e):  0.176 uncentered,
0.703 centered.  The uncentered value reproduces the 2026-07-29 T3 sweep
(eff 0.05 -> 19/20, 0.15 -> 3/20, 0.25 -> 0/20) with no fitted parameter.

gap(m) approaches that asymptote from BELOW and converges slowly (the
det^(1/D) factor carries n^(-1/D)), so at the m = 12 used throughout Phase 2
the realised wall is well short of it: gap = 1 lands near eff = 0.07
uncentered and eff = 0.21 centered.  Extra signatures therefore do help, but
only logarithmically -- which is why the T4b sweep (m = 8..32) read as a flat
wall.  See check C of glv_hnp_phase2_centered_addendum.py.

REFORMULATIONS TESTED
  cent  — centre the nonce limbs (absorb K1/2 + lam*K2/2 into the A_i
          constants).  Predicted: ~2.7x more headroom in eff at m = 12,
          rising to 4x as m -> infinity.
  proj  — delete the d column (project along e_m), killing the trivial
          vector n*S_D*e_m and the n^2/3 that d contributes to ||v||^2.
          d is then recovered algebraically from (k1_0, k2_0).  This is the
          literal Thread-23 proposal.  MEASURED: a bit-exact no-op (E2b),
          because LLL already isolates the trivial vector and reduces inside
          its orthogonal complement.  Its GH gap is correspondingly
          over-optimistic -- always score this lattice with the UNPROJECTED
          gap (E2).

Run: python3 glv_hnp_phase2_centered.py
Depends: fpylll, sympy.
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) — verbatim from
# glv_hnp_phase2_lambda_threshold.py so results stay comparable.
# ---------------------------------------------------------------------------

def modinv(a, m):
    return pow(a, -1, m)

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        s = 3 * x1 * x1 * modinv(2 * y1, p) % p
    else:
        s = (y2 - y1) * modinv(x2 - x1, p) % p
    x3 = (s * s - x1 - x2) % p
    y3 = (s * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(P, k, p):
    if k == 0: return None
    R, Q = None, P
    while k > 0:
        if k & 1: R = ec_add(R, Q, p)
        Q = ec_add(Q, Q, p)
        k >>= 1
    return R

def tonelli_shanks(n, p):
    n %= p
    if n == 0: return 0
    if pow(n, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    mm, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0 and ec_mul((x, y), n, p) is None:
            return (x, y)
    return None

# ---------------------------------------------------------------------------
# j=0 CM structure
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in (a + s, a - s):
            if num % 2 == 0:
                b = num // 2
                if b >= 0 and a * a - a * b + b * b == p:
                    return (a, b)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    sq = tonelli_shanks((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None
    return (min(r1, r2), max(r1, r2))

def lam_star(lam, n):
    return min(lam, n - lam) / n

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        if ec_mul((x, y), n, p) is None:
            G = find_generator(p, b, n, seed=seed)
            if G is not None:
                return (p, b, n, None, G)
    return None

def collect_curves(lo, hi, want):
    """Collect j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3."""
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    n_cand = p + 1 - t
                    if n_cand < 2 or n_cand % 3 != 1 or not sympy.isprime(n_cand):
                        continue
                    roots = glv_roots(n_cand)
                    if roots is None:
                        continue
                    cur = build_curve(p, n_cand)
                    if cur is None:
                        continue
                    out.append((p, cur[1], n_cand, roots[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out

# ---------------------------------------------------------------------------
# Lattice construction — baseline (2026-06-15 scaling) plus the two variants
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return n // k1_bound, 1, max(1, n // k2_bound), n

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42):
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    while len(sigs) < m and attempts < 200000:
        attempts += 1
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        k_full = (k1 + lam * k2) % n
        if k_full == 0: continue
        R = ec_mul(G, k_full, p)
        if R is None: continue
        r = R[0] % n
        if r == 0: continue
        h = rng.randint(0, n - 1)
        s = modinv(k_full, n) * (h + d_secret * r) % n
        if s == 0: continue
        s_inv = modinv(s, n)
        A = h * s_inv % n
        B = r * s_inv % n
        assert (A + B * d_secret) % n == k_full
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2})
    return sigs

def build_lattice(sigs, n, lam, k1_bound, k2_bound, centered, projected):
    """
    Columns (uncentred layout, dim 2m+2):
        0..m-1     k1 limbs   (scale S_K1)
        m          d          (scale S_D)      <- deleted when projected
        m+1..2m    k2 limbs   (scale S_K2)
        2m+1       Kannan     (scale S_KANNAN)
    When `centered`, the constants A_i are shifted by -(K1/2 + lam*K2/2) so
    the planted limbs are u_i = k1_i - K1//2 and v_i = k2_i - K2//2.
    """
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    shift = (k1_bound // 2 + lam * (k2_bound // 2)) % n if centered else 0

    dim = 2 * m + 2
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1                       # modular reduction rows
    for i in range(m):                           # d row
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):                           # k2 rows
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):                           # Kannan row
        M[2 * m + 1][i] = (sigs[i]['A'] - shift) % n * S_K1
    M[2 * m + 1][dim - 1] = S_KAN

    if projected:
        # Delete column m (project along e_m).  The d row survives as a
        # generator; the lattice drops to rank 2m+1 in dim 2m+1, and the
        # trivial vector n*S_D*e_m is gone.
        M = [row[:m] + row[m + 1:] for row in M]
    return M, (S_K1, S_D, S_K2, S_KAN)

# ---------------------------------------------------------------------------
# Predictor
# ---------------------------------------------------------------------------

def gh_unit(D):
    """Gaussian-heuristic constant: lambda_1 ~ gh_unit(D) * det^(1/D)."""
    return math.exp((math.lgamma(D / 2.0 + 1) / D)) / math.sqrt(math.pi)

def predict(m, n, k1_bound, k2_bound, centered, projected):
    """Return (gap, log2 det, E||v||, D). Recovery plausible iff gap > 1."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    a = 12.0 if centered else 3.0
    log_det = (m * math.log(n * S_K1) + m * math.log(S_K2) + math.log(S_KAN))
    vsq = (m * (k1_bound ** 2 / a) * S_K1 ** 2
           + m * (k2_bound ** 2 / a) * S_K2 ** 2
           + S_KAN ** 2)
    if projected:
        D = 2 * m + 1
    else:
        D = 2 * m + 2
        log_det += math.log(S_D)
        vsq += (n ** 2 / 3.0) * S_D ** 2
    lam1 = gh_unit(D) * math.exp(log_det / D)
    return lam1 / math.sqrt(vsq), log_det / math.log(2), math.sqrt(vsq), D

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def recover(reduced, sigs, n, lam, k1_bound, k2_bound, centered, projected,
            d_secret):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in reduced:
        if abs(row[-1]) != S_KAN:
            continue
        sign = 1 if row[-1] > 0 else -1
        if not projected:
            d_cand = (sign * row[m]) % n
        else:
            # d is not in the lattice any more: read limb 0 back out and
            # solve A_0 + B_0*d = k1_0 + lam*k2_0 (mod n) for d.
            if sign * row[0] % S_K1 or sign * row[m] % S_K2:
                continue
            u0 = sign * row[0] // S_K1
            v0 = sign * row[m] // S_K2
            k1_0 = u0 + (k1_bound // 2 if centered else 0)
            k2_0 = v0 + (k2_bound // 2 if centered else 0)
            k_full = (k1_0 + lam * k2_0) % n
            if sigs[0]['B'] % n == 0:
                continue
            d_cand = (k_full - sigs[0]['A']) * modinv(sigs[0]['B'], n) % n
        if d_cand and d_cand == d_secret:
            return True
    return False

def run_once(curve, m, d_secret, k1_bound, centered, projected, seed=42,
             use_bkz=False, bkz_beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M, _ = build_lattice(sigs, n, lam, k1_bound, k2_bound, centered, projected)
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    reduced = [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]
    return recover(reduced, sigs, n, lam, k1_bound, k2_bound, centered,
                   projected, d_secret)

def rate(curve, m, k1_bound, seeds, centered, projected, **kw):
    n = curve[2]
    w = t = 0
    for s in seeds:
        d = random.Random(s ^ 0xABCD).randint(1, n - 1)
        r = run_once(curve, m, d, k1_bound, centered, projected, seed=s, **kw)
        if r is None:
            continue
        t += 1
        w += bool(r)
    return w, t

VARIANTS = [('base', False, False), ('proj', False, True),
            ('cent', True, False), ('cent+proj', True, True)]

SEEDS = [42, 1234, 9999, 555, 31337]
M_SIGS = 12

def main():
    print("=" * 78)
    print("Thread 23 — centering + projection of the GLV-HNP Phase-2 lattice")
    print("=" * 78)
    print(f"m = {M_SIGS} signatures, {len(SEEDS)} seeds/point, LLL unless stated")
    print(f"predicted eff wall:  uncentred {3/(2*math.pi*math.e):.4f}   "
          f"centred {12/(2*math.pi*math.e):.4f}")

    print("\nCollecting 17-bit j=0 GLV curves ...")
    CURVES = collect_curves(2 ** 16, 2 ** 17, 20)
    print(f"  got {len(CURVES)} curves, "
          f"n in [{min(c[2] for c in CURVES)}, {max(c[2] for c in CURVES)}]")

    # ===========================================================================
    # E0 — sanity: predicted vs measured planted-vector norm
    # ===========================================================================
    print("\n" + "=" * 78)
    print("E0 — planted-norm model check (measured vs predicted, ratio should be ~1)")
    print("=" * 78)
    print(f"{'n':>7} {'K1':>4} {'variant':>10} {'pred':>12} {'meas':>12} {'ratio':>6}")
    for (p, b, n, lam, G) in CURVES[:2]:
        k2 = math.isqrt(n) + 1
        k1 = max(2, int(0.10 * n / k2))
        d = random.Random(7).randint(1, n - 1)
        sigs = gen_signatures(G, d, M_SIGS, n, lam, p, k1, k2, 42)
        S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2)
        for name, cen, proj in VARIANTS:
            vsq = 0
            for s in sigs:
                u = s['k1'] - (k1 // 2 if cen else 0)
                v = s['k2'] - (k2 // 2 if cen else 0)
                vsq += (u * S_K1) ** 2 + (v * S_K2) ** 2
            if not proj:
                vsq += (d * S_D) ** 2
            vsq += S_KAN ** 2
            pred = predict(M_SIGS, n, k1, k2, cen, proj)[2]
            meas = math.sqrt(vsq)
            print(f"{n:>7} {k1:>4} {name:>10} {pred:>12.3e} {meas:>12.3e} "
                  f"{meas/pred:>6.3f}")

    # ===========================================================================
    # E1 — eff sweep, all four variants
    # ===========================================================================
    print("\n" + "=" * 78)
    print("E1 — success rate vs eff = K1*K2/n, all four variants")
    print("=" * 78)
    EFFS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.55, 0.70, 0.85)
    pts = []          # (variant, projected, gap, wins, trials) per curve x eff
    print(f"{'eff':>6} " + " ".join(f"{v[0]:>12}" for v in VARIANTS))
    print(f"{'':>6} " + " ".join(f"{'wins    gap':>12}" for _ in VARIANTS))
    for eff_t in EFFS:
        cells = []
        for name, cen, proj in VARIANTS:
            w = t = 0
            gaps = []
            for curve in CURVES:
                n = curve[2]
                k2 = math.isqrt(n) + 1
                k1 = max(2, int(eff_t * n / k2))
                g = predict(M_SIGS, n, k1, k2, cen, proj)[0]
                ww, tt = rate(curve, M_SIGS, k1, SEEDS, cen, proj)
                pts.append((name, proj, g, ww, tt))
                gaps.append(g); w += ww; t += tt
            cells.append(f"{w:>3}/{t:<3} {sum(gaps)/len(gaps):>5.2f}")
        print(f"{eff_t:>6.2f} " + " ".join(f"{c:>12}" for c in cells))

    # ===========================================================================
    # E2 — is the gap>1 rule predictive?  (pooled over E1)
    # ===========================================================================
    print("\n" + "=" * 78)
    print("E2 — calibration of the gap rule (per curve x eff, pooled over seeds)")
    print("=" * 78)
    buckets = [(0, .5), (.5, .7), (.7, .85), (.85, 1.0), (1.0, 1.2), (1.2, 1.5),
               (1.5, 2.0), (2.0, 1e9)]
    for label, keep in (("unprojected variants (base, cent)", False),
                        ("projected variants (proj, cent+proj)", True)):
        sel_all = [(g, w, t) for (nm, pr, g, w, t) in pts if pr == keep]
        print(f"\n  {label}")
        print(f"  {'gap range':>14} {'pts':>5} {'seed-level success':>19} "
              f"{'full 5/5 rate':>14}")
        for lo, hi in buckets:
            sel = [(w, t) for g, w, t in sel_all if lo <= g < hi]
            if sel:
                sw = sum(w for w, _ in sel); st = sum(t for _, t in sel)
                full = sum(1 for w, t in sel if w == t and t > 0) / len(sel)
                print(f"  {f'[{lo:.2f},{hi:.2f})':>14} {len(sel):>5} "
                      f"{sw/st:>19.2f} {full:>14.2f}")

    # ===========================================================================
    # E2b — does the projection change ANYTHING?  paired comparison
    # ===========================================================================
    print("\n" + "=" * 78)
    print("E2b — projection effect, paired per (curve, eff, centering)")
    print("=" * 78)
    by_key = {}
    for (nm, pr, g, w, t) in pts:
        by_key.setdefault(nm, []).append((w, t))
    for a, b in (('base', 'proj'), ('cent', 'cent+proj')):
        xs, ys = by_key[a], by_key[b]
        diff = [x[0] - y[0] for x, y in zip(xs, ys)]
        print(f"  {a:>10} vs {b:<10} "
              f"total {sum(x[0] for x in xs):>4} vs {sum(y[0] for y in ys):>4} ; "
              f"paired points differing: {sum(1 for d in diff if d):>3}/{len(diff)}"
              f" ; max |delta| = {max(abs(d) for d in diff)}")

    # ===========================================================================
    # E3 — direct replication of the 2026-07-29 T4 K1 sweep
    # ===========================================================================
    print("\n" + "=" * 78)
    print("E3 — T4 replication: K1 sweep on the two 12-bit curves from 2026-07-29")
    print("=" * 78)
    T4 = []
    for n_want in (2557, 2677):
        p = int(sympy.nextprime(n_want))
        found = None
        pp = int(sympy.nextprime(2000))
        while pp < 4000 and found is None:
            if pp % 3 == 1:
                eis = eisenstein_decompose(pp)
                if eis:
                    for t in j0_traces(*eis):
                        if pp + 1 - t == n_want and sympy.isprime(n_want):
                            roots = glv_roots(n_want)
                            if roots:
                                cur = build_curve(pp, n_want)
                                if cur:
                                    found = (pp, cur[1], n_want, roots[0], cur[4])
                            break
            pp = int(sympy.nextprime(pp))
        if found:
            T4.append(found)

    K1S = (2, 3, 4, 6, 8, 12, 16, 24, 32, 48)
    for curve in T4:
        p, b, n, lam, G = curve
        k2 = math.isqrt(n) + 1
        print(f"\ncurve p={p} n={n} lam={lam} lam*={lam_star(lam,n):.3f} K2={k2}")
        print(f"{'K1':>4} {'eff':>6} " +
              " ".join(f"{v[0]:>11}" for v in VARIANTS))
        for k1 in K1S:
            eff = k1 * k2 / n
            cells = []
            for name, cen, proj in VARIANTS:
                w, t = rate(curve, M_SIGS, k1, SEEDS, cen, proj)
                g = predict(M_SIGS, n, k1, k2, cen, proj)[0]
                cells.append(f"{w}/{t} {g:>5.2f}")
            print(f"{k1:>4} {eff:>6.3f} " + " ".join(f"{c:>11}" for c in cells))

    # ===========================================================================
    # E4 — m-independence of the wall (T4b replication, centred variant)
    # ===========================================================================
    print("\n" + "=" * 78)
    print("E4 — does more data move the wall?  (fixed eff, sweep m)")
    print("=" * 78)
    curve = CURVES[0]
    n = curve[2]
    k2 = math.isqrt(n) + 1
    print(f"curve n={n}, K2={k2}")
    print(f"{'eff':>6} {'K1':>4} {'m':>4} " + " ".join(f"{v[0]:>11}" for v in VARIANTS))
    for eff_t in (0.25, 0.55):
        k1 = max(2, int(eff_t * n / k2))
        for m in (8, 12, 20, 32):
            cells = []
            for name, cen, proj in VARIANTS:
                w, t = rate(curve, m, k1, SEEDS, cen, proj)
                g = predict(m, n, k1, k2, cen, proj)[0]
                cells.append(f"{w}/{t} {g:>5.2f}")
            print(f"{k1*k2/n:>6.3f} {k1:>4} {m:>4} " +
                  " ".join(f"{c:>11}" for c in cells))

    # ===========================================================================
    # E5 — BKZ on the centred lattice past its wall
    # ===========================================================================
    print("\n" + "=" * 78)
    print("E5 — BKZ(20) vs LLL on the centred lattice at/above its wall")
    print("=" * 78)
    print(f"{'eff':>6} {'cent LLL':>10} {'cent BKZ20':>12} {'gap':>6}")
    for eff_t in (0.70, 0.85):
        wl = tl = wb = tb = 0
        gs = []
        for curve in CURVES[:8]:
            n = curve[2]
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(eff_t * n / k2))
            gs.append(predict(M_SIGS, n, k1, k2, True, False)[0])
            a, c = rate(curve, M_SIGS, k1, SEEDS, True, False)
            wl += a; tl += c
            a, c = rate(curve, M_SIGS, k1, SEEDS, True, False,
                        use_bkz=True, bkz_beta=20)
            wb += a; tb += c
        print(f"{eff_t:>6.2f} {f'{wl}/{tl}':>10} {f'{wb}/{tb}':>12} "
              f"{sum(gs)/len(gs):>6.2f}")

    print("\ndone.")



if __name__ == "__main__":
    main()
