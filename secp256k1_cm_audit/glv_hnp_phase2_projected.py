"""
Thread 23 — GLV-HNP Phase 2: remove the trivial vector from the lattice.

Background (RESEARCH_AUTOLAB_LOG.md 2026-07-29, exp T5)
-------------------------------------------------------
The Phase-2 lattice A (glv_hnp_phase2_20bit.py:263, build_glv_lattice) is
(2m+2)-dimensional with scalings

    S_K1 = n/K1,  S_D = 1,  S_K2 = n/K2,  S_KANNAN = n

and planted vector

    v = (k1_i*S_K1 | d*S_D | k2_i*S_K2 | S_KANNAN),   ||v|| ~ n*sqrt(2m/3 + 4/3).

T5 showed that on every curve tested LLL's shortest vector is *not* v but the
trivial vector

    t0 = n*S_D*e_m = (0,...,0 | n | 0,...,0 | 0),     ||t0|| = n*S_D = n,

which is in A because n*(d-row) reduces to 0 in the first m coordinates.  t0
carries no information: it is exactly the ambiguity d <-> d+n, and d is only
defined mod n.  So v is never lambda_1 and recovery is a coset/BDD condition,
not an SVP condition.

This script implements the fix proposed as Thread 23: quotient the trivial
direction out.  Since t0 spans A cap <e_m>, the quotient A/<t0> is the
orthogonal projection of A along e_m, i.e. simply *delete the d-column*:

    Lattice B (dim 2m+1), generators (2m+2 of them, rank 2m+1):
      n*S_K1*e_i                      i < m
      b_d   = (B_i*S_K1)_i | 0_m | 0
      b_k2i = -lam*S_K1*e_i + S_K2*e_{m+i}
      b_kan = (A_i*S_K1)_i | 0_m | S_KANNAN

    planted   v'' = (k1_i*S_K1 | k2_i*S_K2 | S_KANNAN)

d is no longer a lattice coordinate; it is recovered algebraically from the
recovered (k1_0, k2_0):

    d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  (mod n).

Lattice C interpolates: lattice A with S_K1, S_K2, S_KANNAN multiplied by T
(equivalently S_D = 1/T).  T = 1 is A; T -> infinity approaches B.

Falsifier stated in the 2026-07-29 log entry:
  * if sv/pv rises above 1 (planted becomes lambda_1) AND the K1 wall of exp T4
    moves outward on the lam*=0.07 curve (currently K1 ~ 4-6), the
    reformulation is a real improvement;
  * if the wall stays at K1 ~ 4-6, the wall is information-theoretic and
    Phase 2 is at its ceiling.

Run: python3 glv_hnp_phase2_projected.py
Deps: fpylll, sympy  (cysignals is pulled in by fpylll at import time)
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a = 0) -- verbatim from
# glv_hnp_phase2_20bit.py so the comparison to prior runs is exact.
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
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while True:
        if t == 0: return 0
        if t == 1: return r
        i, tmp = 0, t
        while tmp != 1: tmp = tmp * tmp % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(10000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def eisenstein_decompose(p):
    for a in range(1, 2 * math.isqrt(p // 3) + 3):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            break
        s = math.isqrt(disc)
        if s * s != disc:
            continue
        for num in [a + s, a - s]:
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
    return min(lam % n, n - (lam % n)) / n

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

# ---------------------------------------------------------------------------
# Signature generation -- verbatim from glv_hnp_phase2_20bit.py:236
# ---------------------------------------------------------------------------

def gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed=42):
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
        sigs.append({'A': A, 'B': B, 'k1': k1, 'k2': k2, 'k_full': k_full})
    return sigs

# ---------------------------------------------------------------------------
# Lattice A -- the 2026-06-15 construction, with an added scale knob T.
#   T = 1 reproduces build_glv_lattice() in glv_hnp_phase2_20bit.py:263 exactly.
# ---------------------------------------------------------------------------

def build_lattice_A(sigs, n, lam, k1_bound, k2_bound, T=1):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1 = (n // k1_bound) * T
    S_D = 1
    S_K2 = max(1, n // k2_bound) * T
    S_KANNAN = n * T

    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    for i in range(m):
        M[2 * m + 1][i] = sigs[i]['A'] * S_K1
    M[2 * m + 1][dim - 1] = S_KANNAN
    return M, dict(S_K1=S_K1, S_D=S_D, S_K2=S_K2, S_KAN=S_KANNAN, m=m, dim=dim)

def planted_A(sigs, d, n, sc):
    m, dim = sc['m'], sc['dim']
    v = [0] * dim
    for i in range(m):
        v[i] = sigs[i]['k1'] * sc['S_K1']
        v[m + 1 + i] = sigs[i]['k2'] * sc['S_K2']
    v[m] = d * sc['S_D']
    v[dim - 1] = sc['S_KAN']
    return v

def recover_A(rows, sigs, n, sc, d_secret):
    m, dim = sc['m'], sc['dim']
    for row in rows:
        last = row[dim - 1]
        if abs(last) != sc['S_KAN']:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0:
            continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Lattice B -- lattice A projected along e_m (the d-column deleted).
# Rank 2m+1, given by 2m+2 generators (b_d is rationally in the span of the
# n-rows, which is exactly why the trivial vector disappears).
# ---------------------------------------------------------------------------

def build_lattice_B(sigs, n, lam, k1_bound, k2_bound, T=1):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1 = (n // k1_bound) * T
    S_K2 = max(1, n // k2_bound) * T
    S_KANNAN = n * T

    rows = []
    for i in range(m):                                    # n-rows
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim                                         # b_d (no d column)
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                                    # k2-rows
        r = [0] * dim; r[i] = -lam * S_K1; r[m + i] = S_K2; rows.append(r)
    r = [0] * dim                                         # Kannan row
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[dim - 1] = S_KANNAN
    rows.append(r)
    return rows, dict(S_K1=S_K1, S_K2=S_K2, S_KAN=S_KANNAN, m=m, dim=dim)

def planted_B(sigs, sc):
    m, dim = sc['m'], sc['dim']
    v = [0] * dim
    for i in range(m):
        v[i] = sigs[i]['k1'] * sc['S_K1']
        v[m + i] = sigs[i]['k2'] * sc['S_K2']
    v[dim - 1] = sc['S_KAN']
    return v

def recover_B(rows, sigs, n, lam, sc, d_secret):
    """d is not a coordinate: read (k1_0, k2_0) off the vector and solve
       k1_0 + lam*k2_0 = A_0 + B_0*d (mod n) for d."""
    m, dim = sc['m'], sc['dim']
    A0, B0 = sigs[0]['A'], sigs[0]['B']
    B0inv = modinv(B0, n) if math.gcd(B0, n) == 1 else None
    if B0inv is None:
        return None
    for row in rows:
        last = row[dim - 1]
        if abs(last) != sc['S_KAN']:
            continue
        sign = 1 if last > 0 else -1
        x0 = sign * row[0]
        y0 = sign * row[m]
        if x0 % sc['S_K1'] or y0 % sc['S_K2']:
            continue
        k1_0 = x0 // sc['S_K1']
        k2_0 = y0 // sc['S_K2']
        d_cand = (k1_0 + lam * k2_0 - A0) * B0inv % n
        if d_cand == 0:
            continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Reduction driver
# ---------------------------------------------------------------------------

def nrm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def reduce_rows(rows, use_bkz=False, beta=20):
    nr, nc = len(rows), len(rows[0])
    Aint = IntegerMatrix(nr, nc)
    for i in range(nr):
        for j in range(nc):
            Aint[i, j] = int(rows[i][j])
    if use_bkz:
        BKZ.reduction(Aint, BKZ.Param(beta))
    else:
        LLL.reduction(Aint)
    return [[Aint[i][j] for j in range(nc)] for i in range(nr)]

def shortest_nonzero(rows):
    best, bn = None, float('inf')
    for r in rows:
        if any(r):
            v = nrm(r)
            if v < bn:
                bn, best = v, r
    return best, bn

def trial(curve, m, d_secret, k1_bound, seed, lattice='A', T=1,
          use_bkz=False, beta=20, want_geo=False):
    """Returns dict(ok, pv, sv, triv) or None if signatures could not be made."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, b, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    if lattice == 'A':
        rows, sc = build_lattice_A(sigs, n, lam, k1_bound, k2_bound, T)
        pv = planted_A(sigs, d_secret, n, sc)
    else:
        rows, sc = build_lattice_B(sigs, n, lam, k1_bound, k2_bound, T)
        pv = planted_B(sigs, sc)
    red = reduce_rows(rows, use_bkz, beta)
    if lattice == 'A':
        ok = recover_A(red, sigs, n, sc, d_secret) is not None
    else:
        ok = recover_B(red, sigs, n, lam, sc, d_secret) is not None
    out = {'ok': ok}
    if want_geo:
        sv, svn = shortest_nonzero(red)
        pvn = nrm(pv)
        triv = False
        if lattice == 'A' and sv is not None:
            mm, dim = sc['m'], sc['dim']
            triv = (abs(sv[mm]) == n * sc['S_D'] and
                    all(sv[j] == 0 for j in range(dim) if j != mm))
        out.update(pv=pvn, sv=svn, ratio=(svn / pvn if pvn else float('nan')),
                   triv=triv, sc=sc)
    return out

def grid(curve, m, k1_bound, seeds, lattice='A', T=1, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    wins, tot = 0, 0
    for s in seeds:
        d_trial = random.Random(s + 7777).randint(1, n - 1)
        r = trial(curve, m, d_trial, k1_bound, s, lattice, T, use_bkz, beta)
        if r is None:
            continue
        tot += 1
        wins += bool(r['ok'])
    return wins, tot

# ===========================================================================

print("=" * 78)
print("Thread 23 — projecting out the trivial vector in the GLV-HNP Phase-2 lattice")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,          p,    b, n,    lam,  K1, m
    ("8-bit/199",     211,  2, 199,  106,  2,  6),
    ("12-bit/2557",   2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677",   2677, 2, 2647, 185,  8,  10),
]
hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None and (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G), k1, m))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P0: is the planted vector lambda_1 after the projection?")
print("-" * 78)
print("sv/pv = ||shortest vector after LLL|| / ||planted vector||.")
print("T5 (2026-07-29) measured sv/pv in [0.34, 0.61] in lattice A, with the")
print("shortest vector always the trivial n*S_D*e_m.  A value >= 1 means the")
print("planted vector is (at least) as short as anything LLL found.\n")
print(f"{'curve':<14} {'K1':>3} {'m':>3} | {'A: sv/pv':>9} {'triv':>5} {'ok':>3} "
      f"| {'B: sv/pv':>9} {'ok':>3}")
print("-" * 78)
p0 = []
for label, curve, k1, m in hist:
    n = curve[2]
    d_trial = random.Random(42 + 7777).randint(1, n - 1)
    ra = trial(curve, m, d_trial, k1, 42, 'A', want_geo=True)
    rb = trial(curve, m, d_trial, k1, 42, 'B', want_geo=True)
    p0.append((label, ra, rb))
    print(f"{label:<14} {k1:>3} {m:>3} | {ra['ratio']:>9.3f} "
          f"{str(ra['triv']):>5} {str(ra['ok']):>3} "
          f"| {rb['ratio']:>9.3f} {str(rb['ok']):>3}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P1: T4 grid replication (2026-07-29) — does the K1 wall move?")
print("-" * 78)
print("Same setup as T4: K2 = isqrt(n)+1 = 52, 5 seeds, m as in HIST.")
print("T4 baseline (lattice A):")
print("  12-bit/2557 (lam*=0.340): 5 5 5 5 5 4 1 0  for K1 = 2 3 4 6 8 12 16 24")
print("  12-bit/2677 (lam*=0.070): 5 5 5 2 0 0 0 0\n")
K1S = [2, 3, 4, 6, 8, 12, 16, 24]
p1 = {}
for label, curve, _k1, m in hist[1:]:
    n, lam = curve[2], curve[3]
    print(f"{label}  (lam*={lam_star(lam, n):.3f}, m={m})")
    for lat in ('A', 'B'):
        row = []
        for k1 in K1S:
            w, t = grid(curve, m, k1, SEEDS, lattice=lat)
            row.append(w)
        p1[(label, lat)] = row
        print(f"   lattice {lat}:  K1= " +
              "  ".join(f"{k:>2}:{w}/5" for k, w in zip(K1S, row)))
    print()

# ---------------------------------------------------------------------------
print("-" * 78)
print("EXP P2: T-interpolation on the 12-bit/2677 wall (lattice C = A scaled)")
print("-" * 78)
print("Lattice A with S_K1,S_K2,S_KANNAN multiplied by T (i.e. S_D = 1/T).")
print("T -> infinity is lattice B.  K1 = 6 and 8 straddle the T4 wall.\n")
label, curve, _k1, m = hist[2]
print(f"{'K1':>3} | " + " ".join(f"{'T='+str(t):>7}" for t in [1, 2, 4, 16, 256, 4096])
      + f" |{'B':>5}")
print("-" * 78)
p2 = {}
for k1 in (4, 6, 8, 12):
    cells = []
    for T in [1, 2, 4, 16, 256, 4096]:
        w, t = grid(curve, m, k1, SEEDS, lattice='A', T=T)
        cells.append(w)
    wB, _ = grid(curve, m, k1, SEEDS, lattice='B')
    p2[k1] = (cells, wB)
    print(f"{k1:>3} | " + " ".join(f"{str(c)+'/5':>7}" for c in cells) +
          f" |{str(wB)+'/5':>5}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P3: does more data rescue lattice B at the wall?  (T4b analogue)")
print("-" * 78)
print("2026-07-29 T4b: lattice A at K1=8 on 12-bit/2677 gave 0,0,1,0,1 of 5")
print("for m = 8,12,16,24,32 — more signatures did not help.\n")
label, curve, _k1, _m = hist[2]
print(f"{'m':>3} | {'A':>5} {'B':>5}")
print("-" * 78)
p3 = {}
for mm in (8, 12, 16, 24, 32):
    wa, _ = grid(curve, mm, 8, SEEDS, lattice='A')
    wb, _ = grid(curve, mm, 8, SEEDS, lattice='B')
    p3[mm] = (wa, wb)
    print(f"{mm:>3} | {str(wa)+'/5':>5} {str(wb)+'/5':>5}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P4: fresh 17-bit curve sweep, A vs B at the hard efficiencies")
print("-" * 78)
print("T3 (2026-07-29) baseline on 20 fresh 17-bit curves, m=12, 5 seeds:")
print("  eff=0.05 -> 19/20 curves recover 5/5;  eff=0.15 -> 3/20;  eff=0.25 -> 0/20")
print("If B is a real improvement the eff=0.15 and 0.25 columns must rise.\n")

def search_curves_flat(lo, hi, want):
    """j=0 GLV curves with p in [lo,hi), n prime, n = 1 mod 3, lam* spread."""
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    nc = p + 1 - t
                    if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    roots = glv_roots(nc)
                    if roots is None:
                        continue
                    cur = build_curve(p, nc)
                    if cur is None:
                        continue
                    out.append((p, cur[1], nc, roots[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out

curves17 = search_curves_flat(2**16, 2**17, 20)
print(f"collected {len(curves17)} fresh 17-bit curves "
      f"(lam* in [{min(lam_star(c[3], c[2]) for c in curves17):.4f}, "
      f"{max(lam_star(c[3], c[2]) for c in curves17):.4f}])\n")
print(f"{'eff':>6} | {'A: curves 5/5':>14} {'A: trials':>11} "
      f"| {'B: curves 5/5':>14} {'B: trials':>11}")
print("-" * 78)
p4 = {}
for eff in (0.05, 0.15, 0.25):
    stats = {}
    for lat in ('A', 'B'):
        full, tw, tt = 0, 0, 0
        for c in curves17:
            n = c[2]
            k2 = math.isqrt(n) + 1
            k1 = max(2, int(round(eff * n / k2)))
            w, t = grid(c, 12, k1, SEEDS, lattice=lat)
            tw += w; tt += t
            if t and w == t:
                full += 1
        stats[lat] = (full, tw, tt)
    p4[eff] = stats
    fa, wa, ta = stats['A']; fb, wb, tb = stats['B']
    print(f"{eff:>6.2f} | {str(fa)+'/'+str(len(curves17)):>14} {str(wa)+'/'+str(ta):>11} "
          f"| {str(fb)+'/'+str(len(curves17)):>14} {str(wb)+'/'+str(tb):>11}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P5: WHY is B identical to A?  The trivial vector is inert.")
print("-" * 78)
print("Claim: t0 = n*S_D*e_m is axis-aligned, so LLL splits A as {t0} + LLL(B).")
print("(a) per-trial agreement over the whole P1 grid;")
print("(b) structural: strip t0 from LLL(A), delete column m, compare the")
print("    lattice generated with LLL(B) via row-style Hermite normal form.\n")

def hnf(rows):
    from sympy import Matrix
    from sympy.matrices.normalforms import hermite_normal_form
    nz = [r for r in rows if any(r)]
    return hermite_normal_form(Matrix(nz).T)

agree = dis = 0
for label, curve, _k1, m in hist[1:]:
    n = curve[2]
    for k1 in K1S:
        for s in SEEDS:
            d_trial = random.Random(s + 7777).randint(1, n - 1)
            ra = trial(curve, m, d_trial, k1, s, 'A')
            rb = trial(curve, m, d_trial, k1, s, 'B')
            if ra is None or rb is None:
                continue
            if ra['ok'] == rb['ok']:
                agree += 1
            else:
                dis += 1
                print(f"    DISAGREE {label} K1={k1} seed={s}: A={ra['ok']} B={rb['ok']}")
print(f"(a) per-trial: {agree} agree, {dis} disagree "
      f"({100.0*agree/(agree+dis):.1f}% identical)")

print(f"\n(b) {'curve':<14} {'K1':>3} | {'t0 in LLL(A)':>12} {'HNF equal':>10} "
      f"{'norm profile equal':>19}")
p5b = []
for label, curve, k1, m in hist[1:]:
    p, b, n, lam, G = curve
    d_trial = random.Random(42 + 7777).randint(1, n - 1)
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_trial, m, n, lam, p, b, k1, k2b, 42)
    rA, scA = build_lattice_A(sigs, n, lam, k1, k2b)
    rB, scB = build_lattice_B(sigs, n, lam, k1, k2b)
    redA = reduce_rows(rA); redB = reduce_rows(rB)
    dimA = scA['dim']
    t0_rows = [r for r in redA
               if abs(r[m]) == n * scA['S_D']
               and all(r[j] == 0 for j in range(dimA) if j != m)]
    rest = [r for r in redA if r not in t0_rows and any(r)]
    proj = [[r[j] for j in range(dimA) if j != m] for r in rest]
    same_hnf = (hnf(proj) == hnf(redB))
    na = sorted(round(nrm(r), 6) for r in proj)
    nb = sorted(round(nrm(r), 6) for r in redB if any(r))
    same_nrm = (na == nb)
    p5b.append((label, len(t0_rows) == 1, same_hnf, same_nrm))
    print(f"    {label:<14} {k1:>3} | {str(len(t0_rows) == 1):>12} "
          f"{str(same_hnf):>10} {str(same_nrm):>19}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("EXP P6: what DOES set the wall?  Gaussian-heuristic ratio of lattice B.")
print("-" * 78)
print("det(B) = (n*S_K1)^m * S_K2^m * S_KAN / n = n^(2m) / eff^m,  N = 2m+1")
print("GH(B)  = sqrt(N/(2*pi*e)) * det^(1/N)")
print("||pv|| = n * sqrt(2m/3 + 1)   (k1,k2 uniform; d is gone in B)")
print("R = ||pv|| / GH(B).  R < 1 => planted vector is below the heuristic")
print("lambda_1 and unique-SVP/BDD is plausible; R >> 1 => it is not.\n")

def gh_ratio(n, m, eff):
    N = 2 * m + 1
    log_det = 2 * m * math.log(n) - m * math.log(eff)
    gh = math.sqrt(N / (2 * math.pi * math.e)) * math.exp(log_det / N)
    pv = n * math.sqrt(2 * m / 3.0 + 1)
    return pv / gh

print(f"{'curve':<14} {'m':>3} {'K1':>3} {'eff':>7} {'R':>7} {'wins':>6}")
print("-" * 78)
p6 = []
for label, curve, _k1, m in hist[1:]:
    n = curve[2]; k2b = math.isqrt(n) + 1
    for k1, w in zip(K1S, p1[(label, 'B')]):
        eff = k1 * k2b / n
        R = gh_ratio(n, m, eff)
        p6.append((label, m, k1, eff, R, w))
        print(f"{label:<14} {m:>3} {k1:>3} {eff:>7.4f} {R:>7.3f} {str(w)+'/5':>6}")

succ = [R for *_x, R, w in [(a, b, c, d, e, f) for a, b, c, d, e, f in p6] if w == 5]
fail = [R for *_x, R, w in [(a, b, c, d, e, f) for a, b, c, d, e, f in p6] if w == 0]
print(f"\nR over 5/5 cells: [{min(succ):.3f}, {max(succ):.3f}]  (n={len(succ)})")
print(f"R over 0/5 cells: [{min(fail):.3f}, {max(fail):.3f}]  (n={len(fail)})")

print(f"\n17-bit sweep, same statistic (m=12, {len(curves17)} curves):")
print(f"{'eff':>6} {'R (mean)':>9} {'curves 5/5':>11}")
for eff in (0.05, 0.15, 0.25):
    Rs = [gh_ratio(c[2], 12, eff) for c in curves17]
    print(f"{eff:>6.2f} {sum(Rs)/len(Rs):>9.3f} "
          f"{str(p4[eff]['B'][0])+'/'+str(len(curves17)):>11}")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("SUMMARY")
print("=" * 78)
for label, ra, rb in p0:
    print(f"P0 {label:<14} sv/pv  A={ra['ratio']:.3f} (trivial={ra['triv']})  "
          f"B={rb['ratio']:.3f}")
for (label, lat), row in sorted(p1.items()):
    print(f"P1 {label:<14} {lat}: " + " ".join(f"{k}:{w}" for k, w in zip(K1S, row)))
for k1, (cells, wB) in sorted(p2.items()):
    print(f"P2 K1={k1:<3} T=1,2,4,16,256,4096 -> {cells}   B -> {wB}")
for mm, (wa, wb) in sorted(p3.items()):
    print(f"P3 m={mm:<3} A={wa}/5  B={wb}/5")
for eff, st in sorted(p4.items()):
    print(f"P4 eff={eff:.2f} A={st['A'][0]}/{len(curves17)} curves, "
          f"{st['A'][1]}/{st['A'][2]} trials   "
          f"B={st['B'][0]}/{len(curves17)} curves, {st['B'][1]}/{st['B'][2]} trials")
print(f"P5 per-trial agreement A vs B: {agree}/{agree+dis}")
for label, one_t0, same_hnf, same_nrm in p5b:
    print(f"P5 {label:<14} unique t0={one_t0}  HNF(proj LLL(A))==HNF(LLL(B)): "
          f"{same_hnf}  norms equal: {same_nrm}")
print(f"P6 R(5/5 cells) in [{min(succ):.3f}, {max(succ):.3f}]; "
      f"R(0/5 cells) in [{min(fail):.3f}, {max(fail):.3f}]")
print("=" * 78)
