"""
GLV-HNP Phase 2 / Thread 23 — reformulate the lattice so the planted vector
is the shortest vector.

Motivation (RESEARCH_AUTOLAB_LOG.md 2026-07-29, exp T5): in the Phase-2 lattice
of `glv_hnp_phase2_20bit.py:262` the LLL shortest vector is ALWAYS the trivial
vector n*S_D*e_m (100% of its energy in the d-column, |sv[m]|/n = 1.0000 on every
curve tested, sv/pv in [0.337, 0.368]).  Recovery is therefore not an SVP
condition but a BDD/coset condition, which retro-explains why six consecutive
curve-level invariants failed to separate the success/failure classes.

Thread 23 as proposed: quotient out the trivial n*e_m direction and check
whether the planted vector becomes lambda_1.

Construction.  Only row m of the (2m+2)-dim basis has a nonzero d-column entry
(S_D = 1), so the coordinate projection pi that deletes column m has

    ker(pi|_L) = L cap span(e_m) = n*Z*e_m      (n | B_i fails generically)

hence pi(L) has rank 2m+1 and every fibre of pi|_L is a coset of n*Z*e_m: the
d-coordinate of a preimage is still uniquely determined mod n.  No information
is lost.  d is read back algebraically from (k1_i, k2_i):

    k1_i = A_i + d*B_i - lam*k2_i  (mod n)   =>   d = (k1_i - A_i + lam*k2_i)/B_i

Norms: ||pi(v_planted)||^2 = ||v_planted||^2 - d^2, and the trivial vector maps
to 0, so it is gone from the projected lattice entirely.

Falsifier (stated 2026-07-29): the reformulation is a real improvement iff
sv/pv rises to ~1 AND the K1 wall on the lam*=0.07 curve (currently K1 ~ 4-6)
moves outward.  If the wall stays put, the wall is information-theoretic and
Phase 2 is at its ceiling.

Run: python3 glv_hnp_phase2_projected.py
Requires: fpylll, sympy   (pip install fpylll cysignals sympy)
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a=0) -- verbatim from
# glv_hnp_phase2_lambda_threshold.py so numbers are comparable across sessions.
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
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Shared scaling + signatures
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    """(S_K1, S_D, S_K2, S_KANNAN) -- identical to the 2026-06-15 construction."""
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

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

def norm(v):
    return math.sqrt(sum(x * x for x in v))

# ---------------------------------------------------------------------------
# BASELINE lattice (verbatim, dim 2m+2, columns: k1[0..m-1] | d | k2 | kannan)
# ---------------------------------------------------------------------------

def build_lattice_base(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
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
    M[2 * m + 1][dim - 1] = S_KAN
    return M

def planted_base(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_base(rows, m, n, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand and d_cand == d_secret:
            return True
    return False

# ---------------------------------------------------------------------------
# PROJECTED lattice (Thread 23): delete the d-column.
# 2m+2 generators, 2m+1 columns: k1[0..m-1] | k2[m..2m-1] | kannan[2m]
# ---------------------------------------------------------------------------

def build_lattice_proj(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):                                   # n*S_K1*e_i
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim                                        # d-carrier (no d-column)
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):                                   # k2 rows
        r = [0] * dim
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    r = [0] * dim                                        # Kannan row
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[dim - 1] = S_KAN
    rows.append(r)
    return rows

def planted_proj(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_proj(rows, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read d back algebraically from (k1_i, k2_i); the d-column is gone."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        w = [sign * x for x in row]
        # every generator is a multiple of S_K1 / S_K2 in its block
        if any(w[i] % S_K1 for i in range(m)): continue
        if any(w[m + i] % S_K2 for i in range(m)): continue
        for i in range(m):
            B = sigs[i]['B'] % n
            if math.gcd(B, n) != 1: continue
            k1i = w[i] // S_K1
            k2i = w[m + i] // S_K2
            d_cand = (k1i - sigs[i]['A'] + lam * k2i) * modinv(B, n) % n
            if d_cand and d_cand == d_secret:
                return True
    return False

# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------

def lll_rows(M, beta=0):
    A = IntegerMatrix.from_matrix(M)
    if beta:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]

def shortest_nonzero(rows):
    best, bv = None, None
    for r in rows:
        nr = norm(r)
        if nr > 0 and (best is None or nr < best):
            best, bv = nr, r
    return best, bv

def run_once(curve, m, d_secret, k1_bound, seed, mode, beta=0):
    """mode in {'base','proj'}.  Returns (ok, ||planted||, ||sv||, sv_is_planted)."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if mode == 'base':
        M = build_lattice_base(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_base(sigs, d_secret, n, k1_bound, k2_bound)
        red = lll_rows(M, beta)
        ok = recover_base(red, m, n, S_KAN, d_secret)
    else:
        M = build_lattice_proj(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_proj(sigs, d_secret, n, k1_bound, k2_bound)
        red = lll_rows(M, beta)
        ok = recover_proj(red, sigs, n, lam, k1_bound, k2_bound, d_secret)
    sn, sv = shortest_nonzero(red)
    pn = norm(pv)
    is_pl = sv is not None and (sv == pv or [-x for x in sv] == pv)
    return ok, pn, sn, is_pl

def rate(curve, m, k1_bound, seeds, mode, beta=0, d_secret=None):
    p, b, n, lam, G = curve
    wins, ratios, plcnt, tot = 0, [], 0, 0
    for s in seeds:
        d = d_secret if d_secret is not None else (random.Random(s ^ 0xD).randrange(1, n))
        res = run_once(curve, m, d, k1_bound, s, mode, beta)
        if res is None:
            continue
        ok, pn, sn, is_pl = res
        tot += 1
        wins += bool(ok)
        plcnt += bool(is_pl)
        ratios.append(sn / pn if pn else float('nan'))
    mr = sum(ratios) / len(ratios) if ratios else float('nan')
    return wins, tot, mr, plcnt

# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]

# Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

print("=" * 78)
print("Thread 23 - projected Phase-2 lattice (d-column quotiented out)")
print("=" * 78)

hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, label
    assert (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G), k1, m))

# ===========================================================================
# P0 - correctness of the projected construction
# ===========================================================================
print("\n" + "-" * 78)
print("P0: the projected lattice contains pi(v_planted), and d is recoverable")
print("-" * 78)
print(f"{'curve':<18} {'m':>3} {'dim_b':>6} {'dim_p':>6} {'pi(v) in L?':>12} "
      f"{'||v_b||':>12} {'||v_p||':>12} {'drop':>7}")
allok = True
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    d = random.Random(7).randrange(1, n)
    sigs = gen_signatures(G, d, m, n, lam, p, k1, k2b, 42)
    Mp = build_lattice_proj(sigs, n, lam, k1, k2b)
    vp = planted_proj(sigs, d, n, k1, k2b)
    vb = planted_base(sigs, d, n, k1, k2b)
    # membership: solve for the integer combination directly
    S_K1, S_D, S_K2, S_KAN = scales(n, k1, k2b)
    combo = [0] * len(Mp)
    combo[m] = d                       # d-carrier
    for i in range(m):
        combo[m + 1 + i] = sigs[i]['k2']
    combo[2 * m + 1] = 1               # Kannan row
    acc = [0] * (2 * m + 1)
    for c, row in zip(combo, Mp):
        if c:
            for j, x in enumerate(row):
                acc[j] += c * x
    for i in range(m):                 # reduce mod n*S_K1 using rows 0..m-1
        acc[i] %= (n * S_K1)
    member = (acc == vp)
    allok &= member
    print(f"{label:<18} {m:>3} {2*m+2:>6} {2*m+1:>6} {str(member):>12} "
          f"{norm(vb):>12.4e} {norm(vp):>12.4e} {1-norm(vp)/norm(vb):>6.1%}")
print(f"\nP0 verdict: {'PASS' if allok else 'FAIL'}")

# ===========================================================================
# P1 - anchors: baseline vs projected, sv/pv and is-sv-the-planted-vector
# ===========================================================================
print("\n" + "-" * 78)
print("P1: anchors, LLL.  sv/pv = ||LLL b_1|| / ||v_planted||;  pl = #seeds")
print("    where b_1 IS the planted vector (up to sign).  5 seeds.")
print("-" * 78)
print(f"{'curve':<18} {'K1':>3} {'m':>3} | {'base ok':>8} {'sv/pv':>7} {'pl':>3} "
      f"| {'proj ok':>8} {'sv/pv':>7} {'pl':>3}")
for label, curve, k1, m in hist:
    wb, tb, rb, pb = rate(curve, m, k1, SEEDS, 'base')
    wp, tp, rp, pp = rate(curve, m, k1, SEEDS, 'proj')
    print(f"{label:<18} {k1:>3} {m:>3} | {wb:>4}/{tb:<3} {rb:>7.3f} {pb:>3} "
          f"| {wp:>4}/{tp:<3} {rp:>7.3f} {pp:>3}")

# ===========================================================================
# P2 - the K1 wall (T4 replication of 2026-07-29)
# ===========================================================================
print("\n" + "-" * 78)
print("P2: K1 wall.  T4 grid from 2026-07-29 re-run in both lattices.")
print("    baseline expectation: 2557 walls at K1~12-16, 2677 at K1~4-6.")
print("-" * 78)
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
for label, curve, _k1, m in hist[1:]:
    p, b, n, lam, G = curve
    print(f"\n{label}   lam*={lam_star(lam,n):.3f}  n={n}  m={m}")
    hdr = "  ".join(f"{k:>5}" for k in K1_GRID)
    print(f"  {'K1':<6}{hdr}")
    for mode in ('base', 'proj'):
        cells = []
        for k1 in K1_GRID:
            w, t, r, _ = rate(curve, m, k1, SEEDS, mode)
            cells.append(f"{w}/{t}")
        print(f"  {mode:<6}" + "  ".join(f"{c:>5}" for c in cells))
    cells = []
    for k1 in K1_GRID:
        w, t, r, _ = rate(curve, m, k1, SEEDS, 'proj', beta=20)
        cells.append(f"{w}/{t}")
    print(f"  {'pBKZ20':<6}" + "  ".join(f"{c:>5}" for c in cells))

# ===========================================================================
# P3 - does more data rescue the wall?  (T4b replication, K1 = 8)
# ===========================================================================
print("\n" + "-" * 78)
print("P3: T4b replication -- 12-bit/2677 at K1=8, m sweep.  Baseline was")
print("    0,0,1,0,1 of 5 for m = 8,12,16,24,32 (a K1 wall, not an m wall).")
print("-" * 78)
_, curve2677, _, _ = hist[2]
M_GRID = [8, 12, 16, 24, 32]
print(f"  {'m':<6}" + "  ".join(f"{x:>5}" for x in M_GRID))
for mode in ('base', 'proj'):
    cells = []
    for mm in M_GRID:
        w, t, r, _ = rate(curve2677, mm, 8, SEEDS, mode)
        cells.append(f"{w}/{t}")
    print(f"  {mode:<6}" + "  ".join(f"{c:>5}" for c in cells))

# ===========================================================================
# P4 - fresh 17-bit curves at the eff=0.15 regime (baseline scored 3/20)
# ===========================================================================
print("\n" + "-" * 78)
print("P4: fresh 17-bit j=0 GLV curves, eff = K1*K2/n = 0.15, m=12, 5 seeds.")
print("    2026-07-29 T3 baseline at this eff: 3/20 curves recovered 5/5.")
print("-" * 78)

def search_curves(lo, hi, want=12):
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
                    lam = roots[0]
                    rng = random.Random(nc)
                    bb = None
                    for _ in range(400):
                        cand = rng.randint(1, p - 1)
                        x = rng.randint(0, p - 1)
                        rhs = (pow(x, 3, p) + cand) % p
                        y = tonelli_shanks(rhs, p)
                        if y is None or y == 0:
                            continue
                        if ec_mul((x, y), nc, p) is None:
                            bb = cand
                            break
                    if bb is None:
                        continue
                    G = find_generator(p, bb, nc)
                    if G is None:
                        continue
                    out.append((p, bb, nc, lam, G))
                    break
        p = int(sympy.nextprime(p))
    return out

curves17 = search_curves(2 ** 16, 2 ** 17, want=12)
print(f"  found {len(curves17)} curves\n")
print(f"  {'p':>7} {'n':>7} {'lam*':>6} {'K1':>4} {'eff':>6} | "
      f"{'base':>6} {'sv/pv':>7} | {'proj':>6} {'sv/pv':>7} {'pl':>3}")
agg = {'base': 0, 'proj': 0, 'nb': 0}
for cur in curves17:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    k1 = max(2, round(0.15 * n / k2b))
    eff = k1 * k2b / n
    wb, tb, rb, _ = rate(cur, 12, k1, SEEDS, 'base')
    wp, tp, rp, pp = rate(cur, 12, k1, SEEDS, 'proj')
    agg['base'] += (wb == tb and tb > 0)
    agg['proj'] += (wp == tp and tp > 0)
    agg['nb'] += 1
    print(f"  {p:>7} {n:>7} {lam_star(lam,n):>6.3f} {k1:>4} {eff:>6.3f} | "
          f"{wb:>2}/{tb:<3} {rb:>7.3f} | {wp:>2}/{tp:<3} {rp:>7.3f} {pp:>3}")
print(f"\n  curves recovering 5/5:  base {agg['base']}/{agg['nb']}   "
      f"proj {agg['proj']}/{agg['nb']}")

# ===========================================================================
# P5 - what IS the shortest vector now?  (projected analogue of 2026-07-29 T5)
# ===========================================================================
print("\n" + "-" * 78)
print("P5: block-energy decomposition of LLL b_1 in the projected lattice.")
print("    Baseline T5 verdict: b_1 = n*S_D*e_m, d-block energy 1.0000 always.")
print("    'dcar' = b_1 lies in the d-carrier sublattice (k2 block and Kannan")
print("    coord both zero), i.e. b_1 is (c*B_i mod n)*S_K1 for some c.")
print("-" * 78)
print(f"{'curve':<18} {'K1':>3} {'seed':>6} | {'sv/pv':>7} {'k1blk':>6} {'k2blk':>6} "
      f"{'kan':>6} {'dcar':>5} {'ok':>3}")
for label, curve, k1_def, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_def, k2b)
    for seed in SEEDS[:3]:
        d = random.Random(seed ^ 0xD).randrange(1, n)
        sigs = gen_signatures(G, d, m, n, lam, p, k1_def, k2b, seed)
        if len(sigs) < m:
            continue
        M = build_lattice_proj(sigs, n, lam, k1_def, k2b)
        pv = planted_proj(sigs, d, n, k1_def, k2b)
        red = lll_rows(M)
        sn, sv = shortest_nonzero(red)
        e1 = sum(sv[i] ** 2 for i in range(m))
        e2 = sum(sv[m + i] ** 2 for i in range(m))
        e3 = sv[2 * m] ** 2
        tot = e1 + e2 + e3
        dcar = (e2 == 0 and e3 == 0)
        ok = recover_proj(red, sigs, n, lam, k1_def, k2b, d)
        print(f"{label:<18} {k1_def:>3} {seed:>6} | {sn/norm(pv):>7.3f} "
              f"{e1/tot:>6.3f} {e2/tot:>6.3f} {e3/tot:>6.3f} {str(dcar):>5} "
              f"{str(bool(ok)):>5}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
