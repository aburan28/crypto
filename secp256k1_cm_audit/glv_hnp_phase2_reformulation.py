"""
GLV-HNP Phase 2, Thread 23: reformulate the lattice so the planted vector is
the target of the reduction, not a mid-profile vector.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29 run, exp T5):
  The Phase-2 lattice always contains the trivial vector  n*S_D*e_m  (norm
  n*S_D), which is 2-3x SHORTER than the planted vector on every curve tested.
  That run concluded:
    (i)  "no choice of S_D removes it -- both vectors scale linearly in S_D";
    (ii) recovery is therefore a BDD/coset condition, not an SVP condition.
  and proposed Thread 23: project out the e_m direction, or replace the Kannan
  embedding by an explicit CVP call, then re-measure the K1 wall.

  Falsifier stated on 2026-07-29: if the K1 wall on the lam*=0.07 curve
  (currently K1 ~ 4-6) moves outward after the reformulation, the
  reformulation is a real improvement; if it stays, the wall is
  information-theoretic and Phase 2 is at its ceiling.

This script runs that falsifier, plus the algebra that motivates it.

  U0  ALGEBRA.  Claim (i) above is wrong.  Only ONE coordinate of the planted
      vector scales with S_D:
          ||n*S_D*e_m||^2       = (n*S_D)^2
          E||v_planted||^2      = m*(K1*S_K1)^2/3 + (n*S_D)^2/3
                                  + m*(K2*S_K2)^2/3 + S_KANNAN^2
                                ~ n^2 * (2m/3 + S_D^2/3 + 1)
      so the trivial vector is the LONGER of the two as soon as
          (2/3)*S_D^2 > 2m/3 + 1   <=>   S_D^2 > m + 3/2.
      U0 verifies this numerically and reports S_D* = ceil(sqrt(m+1.5)).

  U1  S_D SWEEP.  Does making the planted vector lambda_1 (via S_D > S_D*)
      actually change recovery?  Measures, per S_D: sv/pv, the fraction of the
      shortest vector's energy in the d column, and the LLL success rate.

  U2  d-COLUMN-FREE FORMULATION ("L2m+1").  Drop column m entirely.  d is not
      small, so it never belonged in a BDD instance; it is still uniquely
      recoverable from the k1/k2 columns via
          d = (k1_0 + lam*k2_0 - A_0) * B_0^{-1}  (mod n).
      Lattice (dim 2m+1, k1 cols | k2 cols | Kannan col), generators
          n*S_K1*e_i ;  (B_i*S_K1, 0, 0) ;  (-lam*S_K1*e_i + S_K2*e_{m+i}) ;
          (A_i*S_K1, 0, S_KANNAN).
      This removes the trivial vector n*S_D*e_m from the lattice outright.
      U2 re-runs the 2026-07-29 T4 K1 sweep on both 12-bit curves in both
      formulations -- this is the Thread 23 falsifier.

  U3  CLOSED-FORM PREDICTOR.  For the d-column-free lattice,
          det = (n*S_K1)^m * S_K2^m / n = n^(2m-1) / eff^m,   eff = K1*K2/n
          GH  = sqrt(2m/(2*pi*e)) * det^(1/(2m))
          E||e|| = n*sqrt(2m/3)                (Kannan coord excluded)
      =>  gamma := E||e|| / GH = sqrt(2*pi*e/3) * n^(1/(2m)) * sqrt(eff)
                              ~ 2.4849 * n^(1/(2m)) * sqrt(eff).
      Success predicted iff gamma < 1, i.e.
          eff < 3 / (2*pi*e * n^(1/m)) .
      Note gamma does not involve lam at all -- consistent with the 2026-07-29
      falsification of the lam* threshold.  U3 scores gamma against every
      empirical cell measured in U1/U2 plus a fresh 17-bit grid.

Run: python3 glv_hnp_phase2_reformulation.py
"""

import math
import random
import sympy
from fpylll import IntegerMatrix, LLL, BKZ

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (same as glv_hnp_phase2_lambda_threshold.py)
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
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c = i, b * b % p
        t, r = t * c % p, r * b % p

def find_generator(p, b, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(4000):
        x = rng.randrange(p)
        rhs = (x * x % p * x + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        P = (x, y)
        if ec_mul(P, n, p) is None and ec_mul(P, 1, p) is not None:
            return P
    return None

def eisenstein_decompose(p):
    """p = a^2 - a*b + b^2 for p = 1 mod 3 (norm form of Z[omega])."""
    if p % 3 != 1:
        return None
    lim = int(2 * math.isqrt(p) // 1) + 2
    for a in range(-lim, lim + 1):
        disc = 4 * p - 3 * a * a
        if disc < 0:
            continue
        r = math.isqrt(disc)
        if r * r != disc:
            continue
        if (a + r) % 2 == 0:
            return (a, (a + r) // 2)
    return None

def j0_traces(a, b):
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    neg3 = (n - 3) % n
    sq = tonelli_shanks(neg3, n)
    if sq is None:
        return None, None
    inv2 = modinv(2, n)
    r1 = (n - 1 + sq) * inv2 % n
    r2 = (n - 1 + (n - sq)) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1, r2 = r2, r1
    if (r1 * r1 + r1 + 1) % n != 0:
        return None, None
    return min(r1, r2), max(r1, r2)

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Scales, signatures
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound, S_D=1):
    """(S_K1, S_D, S_K2, S_KANNAN).  S_D is now a free parameter."""
    return (n // k1_bound, S_D, max(1, n // k2_bound), n)

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
    return math.sqrt(sum(float(x) * float(x) for x in v))

# ---------------------------------------------------------------------------
# Formulation A: the historical Phase-2 lattice, S_D generalised
# ---------------------------------------------------------------------------

def build_lattice_A(sigs, n, lam, k1_bound, k2_bound, S_D=1):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
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

def planted_A(sigs, d_secret, n, k1_bound, k2_bound, S_D=1):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_A(reduced, m, n, S_D, S_KAN, d_secret):
    dim = 2 * m + 2
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        val = sign * row[m]
        if val % S_D != 0: continue
        d_cand = (val // S_D) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Formulation B: d column dropped  (Thread 23 reformulation)
# ---------------------------------------------------------------------------

def build_lattice_B(sigs, n, lam, k1_bound, k2_bound):
    """dim = 2m+1: [k1_0..k1_{m-1} | k2_0..k2_{m-1} | Kannan].
    2m+2 generators (rank 2m+1); LLL tolerates the dependency."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    rows = []
    for i in range(m):
        r = [0] * dim; r[i] = n * S_K1; rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['B'] * S_K1
    rows.append(r)
    for i in range(m):
        r = [0] * dim
        r[i] = -lam * S_K1
        r[m + i] = S_K2
        rows.append(r)
    r = [0] * dim
    for i in range(m):
        r[i] = sigs[i]['A'] * S_K1
    r[dim - 1] = S_KAN
    rows.append(r)
    return rows

def planted_B(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_B(reduced, sigs, n, lam, k1_bound, k2_bound, d_secret):
    """Read (k1_i, k2_i) off a Kannan row, then solve for d."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'], n)
    for row in reduced:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        if (sign * row[0]) % S_K1 != 0 or (sign * row[m]) % S_K2 != 0:
            continue
        k1_0 = (sign * row[0]) // S_K1
        k2_0 = (sign * row[m]) // S_K2
        d_cand = (k1_0 + lam * k2_0 - sigs[0]['A']) * B0inv % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Reduction drivers
# ---------------------------------------------------------------------------

def reduce_rows(rows, use_bkz=False, beta=20):
    A = IntegerMatrix.from_matrix(rows)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(A.ncols)] for i in range(A.nrows)]

def run_A(curve, m, d_secret, k1_bound, seed=42, S_D=1, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m: return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2b, S_D)
    M = build_lattice_A(sigs, n, lam, k1_bound, k2b, S_D)
    red = reduce_rows(M, use_bkz, beta)
    ok = recover_A(red, m, n, S_D, S_KAN, d_secret) is not None
    pv = planted_A(sigs, d_secret, n, k1_bound, k2b, S_D)
    nz = [r for r in red if any(r)]
    sv = min(nz, key=norm)
    # fraction of shortest-vector energy sitting in the d column
    e_d = (float(sv[m]) ** 2) / max(1e-30, sum(float(x) ** 2 for x in sv))
    return dict(ok=ok, pn=norm(pv), sn=norm(sv), e_d=e_d,
                triv=abs(sv[m]) == n * S_D and all(x == 0 for j, x in enumerate(sv) if j != m))

def run_B(curve, m, d_secret, k1_bound, seed=42, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2b, seed)
    if len(sigs) < m: return None
    rows = build_lattice_B(sigs, n, lam, k1_bound, k2b)
    red = reduce_rows(rows, use_bkz, beta)
    ok = recover_B(red, sigs, n, lam, k1_bound, k2b, d_secret) is not None
    pv = planted_B(sigs, n, k1_bound, k2b)
    nz = [r for r in red if any(r)]
    sv = min(nz, key=norm)
    return dict(ok=ok, pn=norm(pv), sn=norm(sv))

def rate(fn, curve, m, k1_bound, seeds, **kw):
    p, b, n, lam, G = curve
    wins, ratios = 0, []
    for s in seeds:
        d = random.Random(s + 7777).randint(1, n - 1)
        r = fn(curve, m, d, k1_bound, s, **kw)
        if r is None: continue
        wins += bool(r['ok'])
        ratios.append(r['sn'] / r['pn'] if r['pn'] else float('nan'))
    return wins, len(seeds), (sum(ratios) / len(ratios) if ratios else float('nan'))

# ---------------------------------------------------------------------------
# Predictor
# ---------------------------------------------------------------------------

def gamma_predictor(n, m, k1_bound, k2_bound):
    """gamma = E||e|| / GH for the d-column-free lattice (formulation B)."""
    S_K1 = n // k1_bound
    S_K2 = max(1, n // k2_bound)
    dim = 2 * m
    log_det = m * math.log(n * S_K1) + m * math.log(S_K2) - math.log(n)
    log_gh = 0.5 * math.log(dim / (2 * math.pi * math.e)) + log_det / dim
    log_e = math.log(n) + 0.5 * math.log(2 * m / 3.0)
    return math.exp(log_e - log_gh)

def gamma_closed_form(n, m, k1_bound, k2_bound):
    eff = k1_bound * k2_bound / n
    return math.sqrt(2 * math.pi * math.e / 3.0) * n ** (1.0 / (2 * m)) * math.sqrt(eff)

# ---------------------------------------------------------------------------
# Curve construction
# ---------------------------------------------------------------------------

def build_curve(p, n, seed=12345):
    rng = random.Random(seed)
    for _ in range(400):
        b = rng.randint(1, p - 1)
        rhs_ok = False
        for _ in range(30):
            x = rng.randrange(p)
            rhs = (x * x % p * x + b) % p
            y = tonelli_shanks(rhs, p)
            if y is not None and y != 0:
                rhs_ok = True
                if ec_mul((x, y), n, p) is None:
                    G = find_generator(p, b, n)
                    if G is not None:
                        lam, _ = glv_roots(n)
                        if lam is None: return None
                        return (p, b, n, lam, G)
                break
        if not rhs_ok:
            continue
    return None

def search_curves(lo, hi, want=8, seed=999):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                for t in j0_traces(*eis):
                    nc = p + 1 - t
                    if nc < 4 or not sympy.isprime(nc) or nc % 3 != 1:
                        continue
                    lam, _ = glv_roots(nc)
                    if lam is None:
                        continue
                    cur = build_curve(p, nc, seed)
                    if cur is not None:
                        out.append(cur)
                        break
        p = int(sympy.nextprime(p))
    return out


# ===========================================================================
print("=" * 78)
print("Thread 23 — reformulating the Phase-2 GLV-HNP lattice")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]
hist = []
for label, p, b, n, lam, k1, m in HIST:
    G = find_generator(p, b, n)
    assert G is not None, label
    assert (lam * lam + lam + 1) % n == 0, label
    hist.append((label, (p, b, n, lam, G), k1, m))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("U0: does S_D remove the trivial vector n*S_D*e_m?  (corrects 2026-07-29 T5)")
print("-" * 78)
print("||n*S_D*e_m|| = n*S_D                       (linear in S_D)")
print("E||v_planted||^2 ~ n^2*(2m/3 + S_D^2/3 + 1)  (only ONE of 2m+2 coords scales)")
print("=> trivial is longer iff S_D^2 > m + 3/2\n")
print(f"{'curve':<18} {'m':>3} {'S_D*':>5} {'triv(S_D=1)':>12} {'pv(S_D=1)':>11} "
      f"{'triv(S_D*)':>11} {'pv(S_D*)':>10}")
u0_rows = []
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    sdstar = math.ceil(math.sqrt(m + 1.5))
    k2b = math.isqrt(n) + 1
    d0 = random.Random(42 + 7777).randint(1, n - 1)
    sg = gen_signatures(G, d0, m, n, lam, p, k1, k2b, 42)
    pv1 = norm(planted_A(sg, d0, n, k1, k2b, 1))
    pvs = norm(planted_A(sg, d0, n, k1, k2b, sdstar))
    print(f"{label:<18} {m:>3} {sdstar:>5} {float(n*1):>12.4g} {pv1:>11.4g} "
          f"{float(n*sdstar):>11.4g} {pvs:>10.4g}")
    u0_rows.append((label, m, sdstar, n * 1 < pv1, n * sdstar > pvs))
print()
for label, m, sdstar, a, bq in u0_rows:
    print(f"  {label:<18} S_D=1: trivial shorter than planted = {a}   "
          f"S_D={sdstar}: trivial longer than planted = {bq}")

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("U1: S_D sweep in formulation A — does making the planted vector lambda_1 help?")
print("-" * 78)
print("sv/pv = ||LLL shortest|| / ||planted||;  E_d = fraction of sv energy in col m")
u1_rows = []
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    sdstar = math.ceil(math.sqrt(m + 1.5))
    grid = sorted({1, 2, sdstar, 2 * sdstar, 4 * sdstar, 16 * sdstar})
    print(f"\n  {label}  (m={m}, K1={k1}, lam*={lam_star(lam,n):.4f}, S_D*={sdstar})")
    print(f"    {'S_D':>5} {'LLL':>6} {'sv/pv':>8} {'E_d':>7} {'sv==trivial':>12}")
    for S_D in grid:
        wins = 0; rr = []; ed = []; tv = 0; tot = 0
        for s in SEEDS:
            d = random.Random(s + 7777).randint(1, n - 1)
            r = run_A(curve, m, d, k1, s, S_D=S_D)
            if r is None: continue
            tot += 1
            wins += bool(r['ok']); rr.append(r['sn'] / r['pn'])
            ed.append(r['e_d']); tv += bool(r['triv'])
        print(f"    {S_D:>5} {str(wins)+'/'+str(tot):>6} "
              f"{sum(rr)/len(rr):>8.3f} {sum(ed)/len(ed):>7.3f} {str(tv)+'/'+str(tot):>12}")
        u1_rows.append((label, m, k1, S_D, wins, tot))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("U2: THREAD 23 FALSIFIER — K1 wall, formulation A (S_D=1) vs B (d col dropped)")
print("-" * 78)
K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
u2_rows = []
for label, curve, k1_hist, m_hist in hist[1:]:
    p, b, n, lam, G = curve
    m = 12
    sdstar = math.ceil(math.sqrt(m + 1.5))
    print(f"\n  {label}  (m={m}, lam*={lam_star(lam,n):.4f}, n={n})")
    hdr = "    " + f"{'formulation':<22}" + "".join(f"{('K1='+str(k)):>8}" for k in K1_GRID)
    print(hdr)
    for tag, fn, kw in [("A  S_D=1 (baseline)", run_A, dict(S_D=1)),
                        (f"A  S_D={sdstar} (U1 fix)", run_A, dict(S_D=sdstar)),
                        ("B  d-col dropped", run_B, dict())]:
        cells = []
        for K1 in K1_GRID:
            w, t, _ = rate(fn, curve, m, K1, SEEDS, **kw)
            cells.append(f"{w}/{t}")
            u2_rows.append((label, tag, K1, m, w, t,
                            gamma_closed_form(n, m, K1, math.isqrt(n) + 1)))
        print("    " + f"{tag:<22}" + "".join(f"{c:>8}" for c in cells))
    print("    " + f"{'gamma (predictor)':<22}"
          + "".join(f"{gamma_closed_form(n, m, K, math.isqrt(n)+1):>8.2f}" for K in K1_GRID))
    # is the planted vector lambda_1 in formulation B?
    cells = []
    for K1 in K1_GRID:
        _, _, mr = rate(run_B, curve, m, K1, SEEDS)
        cells.append(f"{mr:.3f}")
    print("    " + f"{'B sv/pv (1.0=planted)':<22}" + "".join(f"{c:>8}" for c in cells))

# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("U3: closed-form predictor gamma = sqrt(2*pi*e/3) * n^(1/2m) * sqrt(eff)")
print("-" * 78)
print("Sanity: numeric GH-based gamma vs closed form (should agree to ~1e-3)")
for label, curve, k1, m in hist:
    p, b, n, lam, G = curve
    k2b = math.isqrt(n) + 1
    g1 = gamma_predictor(n, m, k1, k2b)
    g2 = gamma_closed_form(n, m, k1, k2b)
    print(f"  {label:<18} m={m:>3} K1={k1:<3} numeric={g1:.4f}  closed={g2:.4f}  "
          f"rel.diff={abs(g1-g2)/g2:.2e}")

print("\nFresh 17-bit validation grid (curves independent of U1/U2):")
fresh = search_curves(2 ** 16, 2 ** 17, want=12)
print(f"  found {len(fresh)} 17-bit j=0 GLV curves")
grid_rows = []
M_SIGS = 12
for cur in fresh:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    for eff_t in [0.03, 0.05, 0.08, 0.12, 0.16, 0.20, 0.25, 0.30]:
        K1 = max(2, int(eff_t * n / k2b))
        w, t, _ = rate(run_B, cur, M_SIGS, K1, SEEDS)
        wA, tA, _ = rate(run_A, cur, M_SIGS, K1, SEEDS, S_D=1)
        g = gamma_closed_form(n, M_SIGS, K1, k2b)
        grid_rows.append((n, lam_star(lam, n), K1, K1 * k2b / n, g, w, t, wA, tA))

print(f"\n  {'n':>7} {'lam*':>7} {'K1':>5} {'eff':>7} {'gamma':>7} "
      f"{'B wins':>7} {'A wins':>7}")
for n_, ls, K1, eff, g, w, t, wA, tA in grid_rows:
    print(f"  {n_:>7} {ls:>7.4f} {K1:>5} {eff:>7.4f} {g:>7.3f} "
          f"{str(w)+'/'+str(t):>7} {str(wA)+'/'+str(tA):>7}")

# scoring: does gamma < 1 predict success?
def score(rows, wi, ti):
    tp = fp = tn = fn_ = 0
    for r in rows:
        g = r[4]; succ = r[wi] > r[ti] / 2.0
        pred = g < 1.0
        if pred and succ: tp += 1
        elif pred and not succ: fp += 1
        elif not pred and not succ: tn += 1
        else: fn_ += 1
    acc = (tp + tn) / max(1, len(rows))
    return tp, fp, tn, fn_, acc

for name, wi, ti in [("B (d-col dropped)", 5, 6), ("A (baseline)", 7, 8)]:
    tp, fp, tn, fn_, acc = score(grid_rows, wi, ti)
    print(f"\n  gamma<1 as a majority-success predictor, {name}: "
          f"TP={tp} FP={fp} TN={tn} FN={fn_}  accuracy={acc:.3f}")
    base = max(sum(1 for r in grid_rows if r[wi] > r[ti] / 2.0),
               sum(1 for r in grid_rows if r[wi] <= r[ti] / 2.0)) / max(1, len(grid_rows))
    print(f"  majority baseline = {base:.3f}")

# empirical gamma boundary
succ_g = sorted(r[4] for r in grid_rows if r[5] > r[6] / 2.0)
fail_g = sorted(r[4] for r in grid_rows if r[5] <= r[6] / 2.0)
if succ_g and fail_g:
    print(f"\n  formulation B: max gamma among successes = {max(succ_g):.3f}")
    print(f"                 min gamma among failures  = {min(fail_g):.3f}")
    print(f"                 separated = {max(succ_g) < min(fail_g)}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("U4: the m-crossing test — gamma's ONLY non-trivial prediction")
print("-" * 78)
print("gamma = 2.4849 * n^(1/2m) * sqrt(eff) decreases with m at FIXED eff, so at")
print("an eff inside the mixed band the attack must go from failure to success as")
print("m grows, crossing gamma=1 at")
print("      m* = ln(n) / (2 * ln(1 / (2.4849*sqrt(eff)))) .")
print("This is a parameter-free prediction: eff sets m*, nothing is fitted.")
print("It also decides the 2026-07-29 T4b observation (m=8..32 at K1=8 did not")
print("rescue the lam*=0.07 curve) — gamma says T4b never crossed 1 at all.\n")

def m_star(n, eff):
    c = 2.4849 * math.sqrt(eff)
    if c >= 1.0:
        return float('inf')
    return math.log(n) / (2.0 * math.log(1.0 / c))

M_GRID = [6, 8, 10, 12, 16, 20, 24, 28, 32, 40]
u4 = []
for cur in fresh[:3]:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    for eff_t in [0.05, 0.10, 0.16]:
        K1 = max(2, int(round(eff_t * n / k2b)))
        eff = K1 * k2b / n
        ms = m_star(n, eff)
        print(f"  n={n}  lam*={lam_star(lam,n):.4f}  K1={K1}  eff={eff:.4f}  "
              f"predicted m* = {ms:.1f}")
        print(f"    {'m':>5}" + "".join(f"{v:>8}" for v in M_GRID))
        cells, gams = [], []
        for m in M_GRID:
            w, t, _ = rate(run_B, cur, m, K1, SEEDS)
            cells.append(f"{w}/{t}")
            gams.append(gamma_closed_form(n, m, K1, k2b))
            u4.append((n, eff, m, w, t, gams[-1], ms))
        print(f"    {'B':>5}" + "".join(f"{c:>8}" for c in cells))
        print(f"    {'gamma':>5}" + "".join(f"{g:>8.2f}" for g in gams))
        # observed crossing: smallest m with majority success
        obs = None
        for (nn, ee, m, w, t, g, _s) in u4[-len(M_GRID):]:
            if w > t / 2.0:
                obs = m
                break
        print(f"    observed crossing m_obs = {obs}   (predicted m* = {ms:.1f})\n")

print("  summary: predicted vs observed crossing")
print(f"  {'n':>7} {'eff':>7} {'m*':>8} {'m_obs':>7} {'verdict':>10}")
seen = set()
for (n_, eff, m, w, t, g, ms) in u4:
    key = (n_, round(eff, 4))
    if key in seen: continue
    seen.add(key)
    rows = [r for r in u4 if (r[0], round(r[1], 4)) == key]
    obs = next((r[2] for r in sorted(rows, key=lambda r: r[2]) if r[3] > r[4] / 2.0), None)
    if obs is None:
        verdict = "no cross" if ms > M_GRID[-1] else "MISS"
    elif ms == float('inf'):
        verdict = "MISS"
    else:
        verdict = "ok" if abs(obs - ms) <= max(4.0, 0.5 * ms) else "MISS"
    print(f"  {n_:>7} {eff:>7.4f} {ms:>8.1f} {str(obs):>7} {verdict:>10}")


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("U5: is the residual wall reduction-quality or information-theoretic?")
print("-" * 78)
print("U4 found curves at the SAME (n, eff) with wildly different crossings")
print("(m_obs = 8 / 16 / none at eff~0.10).  Take the hardest one and push both")
print("m and the reduction strength.  If BKZ(20) at large m still fails while a")
print("sibling curve succeeds at m=8, the residue is curve-structural, not")
print("reduction-quality.\n")

deep = fresh[:3]
M_DEEP = [40, 48, 64, 80]
for cur in deep:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    K1 = max(2, int(round(0.10 * n / k2b)))
    eff = K1 * k2b / n
    print(f"  n={n}  lam*={lam_star(lam,n):.4f}  K1={K1}  eff={eff:.4f}")
    print(f"    {'m':>7}" + "".join(f"{v:>9}" for v in M_DEEP)
          + "   (dim = 2m+1)")
    for tag, kw in [("LLL", dict()), ("BKZ(20)", dict(use_bkz=True, beta=20))]:
        cells = []
        for m in M_DEEP:
            if tag != "LLL" and m > 48:
                cells.append("-")
                continue
            w, t, _ = rate(run_B, cur, m, K1, SEEDS, **kw)
            cells.append(f"{w}/{t}")
        print(f"    {tag:>7}" + "".join(f"{c:>9}" for c in cells))
    print(f"    {'gamma':>7}"
          + "".join(f"{gamma_closed_form(n, m, K1, k2b):>9.2f}" for m in M_DEEP))
    print()


# ---------------------------------------------------------------------------
print("\n" + "-" * 78)
print("U6: is the residue a CURVE invariant or a SIGNATURE-SAMPLE effect?")
print("-" * 78)
print("Six curve-level invariants failed between 2026-06-21 and 2026-06-29, plus")
print("lam* (2026-07-29 T3) and rho (T2).  Before proposing a seventh, test whether")
print("the quantity being predicted is a function of the curve AT ALL.")
print("n is prime, so every non-identity point generates; changing G re-randomises")
print("every (A_i, B_i) while leaving (p, n, lam) — every curve invariant — fixed.")
print("If the outcome moves with G, no function of the curve can predict it.\n")

WIDE = [42, 1234, 9999, 555, 31337, 7, 101, 2718, 8191, 65537,
        13, 271828, 4444, 90210, 31415, 1618, 2024, 777, 5150, 99991]
for cur in fresh[:3]:
    p, b, n, lam, G = cur
    k2b = math.isqrt(n) + 1
    K1 = max(2, int(round(0.10 * n / k2b)))
    m = 16
    w, t, _ = rate(run_B, cur, m, K1, WIDE)
    print(f"  n={n} lam*={lam_star(lam,n):.4f} K1={K1} m={m}: "
          f"{w}/{t} over 20 signature seeds, generator fixed")

print("\n  same n, same lam, DIFFERENT generator G (5 seeds each):")
print(f"  {'n':>7} {'G-seed':>8} {'G.x':>8} " + "".join(f"{'m='+str(v):>7}" for v in [12, 16, 24]))
for cur in fresh[:3]:
    p, b, n, lam, G0 = cur
    k2b = math.isqrt(n) + 1
    K1 = max(2, int(round(0.10 * n / k2b)))
    for gs in [12345, 2, 3, 4, 5]:
        G = find_generator(p, b, n, seed=gs)
        if G is None:
            continue
        cur2 = (p, b, n, lam, G)
        cells = []
        for m in [12, 16, 24]:
            w, t, _ = rate(run_B, cur2, m, K1, SEEDS)
            cells.append(f"{w}/{t}")
        print(f"  {n:>7} {gs:>8} {G[0]:>8} " + "".join(f"{c:>7}" for c in cells))

print("\n" + "=" * 78)
print("done")
print("=" * 78)
