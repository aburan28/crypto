"""
GLV-HNP Phase 2, Thread 23: does removing the trivial vector help?

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, exp T5):
  The Phase-2 Kannan lattice L (dim 2m+2) always contains the vector
      v_triv = n * S_D * e_m          (= n*row_m - sum_i B_i*row_i)
  which is 2-3x SHORTER than the planted vector on every curve tested.
  So the planted vector is never lambda_1; recovery is a BDD/coset condition
  ("shortest vector among those with last coordinate +-S_KANNAN"), not SVP.

  That run proposed Thread 23: quotient out the trivial direction and re-test.
  Falsifier as stated on 2026-07-29:
     "if sv/pv rises above 1 after the reformulation AND the K1 wall in T4
      moves outward on the lam*=0.07 curve (currently K1 ~ 4-6), the
      reformulation is a real improvement; if the wall stays at K1 ~ 4-6,
      then the wall is information-theoretic and Phase 2 is at its ceiling."

Three variants are compared, all on identical signature sets:

  A "kannan"  baseline, dim 2m+2, columns (k1_i | d | k2_i | Kannan).
              Construction verbatim from glv_hnp_phase2_20bit.py:263.

  B "proj"    dim 2m+1.  Drop the d-column.  Legitimate because
              L cap span(e_m) = Z*n*S_D*e_m, so deleting coordinate m is
              exactly the projection L -> L/(Z*v_triv); det(L') = det(L)/n.
              d carries no information anyway (it is only defined mod n and
              is recoverable from any single signature once (k1_i,k2_i) are
              known, via d = B_i^{-1}(k1_i + lam*k2_i - A_i) mod n).

  C "cvp"     dim 2m.  Drop the Kannan column too and solve the honest CVP
              with an explicit closest-vector call (fplll exact CVP, Babai
              nearest-plane fallback) against target (-A_i*S_K1 | 0).

Experiments:
  E1  geometry: ||pv||, ||sv||, sv/pv, det, Gaussian heuristic, per variant,
      on the three historical curves.  Directly tests "sv/pv rises above 1".
  E2  the T4 K1-wall grid (K1 in 2..24, 5 seeds) for all three variants on
      both historical 12-bit curves.  Directly tests "the wall moves".
  E3  m-sweep at the wall on the lam*=0.07 curve, all three variants.
  E4  Babai-only vs exact CVP vs BKZ-20/40 Kannan — formulation or enumeration?
  E5  the (K1, m) recovery frontier vs the information-theoretic m_thresh.
  E6  distinct-nonce control (draw (k1,k2) without replacement).

RESULTS (2026-08-04 run; full output in glv_hnp_phase2_thread23_output.txt):

  E1  The projection does what it was designed to do geometrically: sv/pv rises
      from 0.586 -> 0.843 (8-bit), 0.517 -> 0.532 (12-bit/2557) and
      0.407 -> 0.813 (12-bit/2677).  But it never crosses 1, so the planted
      vector is still not lambda_1.  Deleting the d-column removes v_triv and
      costs exactly log2(n) of determinant, and the two effects very nearly
      cancel: pv/gh is 1.237 -> 1.206, 1.284 -> 1.290, 1.349 -> 1.301.
      That near-cancellation is the whole story of Thread 23 — see below.

  E2  The K1 wall does NOT move.  All three variants wall at K1 <= 8 on
      12-bit/2557 and K1 <= 4 on 12-bit/2677, reproducing T4 exactly.
      By the 2026-07-29 falsifier, the projection per se is a dead end.

  E3/E4  BUT the exact-CVP variant breaks the T4b ceiling.  At K1=8 on the
      lam*=0.070 curve, where 2026-07-29 concluded "more data does not rescue
      it", cvp reaches 5/5 at m >= 24 while kannan/proj never exceed 2/5 and
      BKZ-40 does not close the gap.  E4 attributes this to the ENUMERATION,
      not the formulation: Babai nearest-plane on the same dim-2m lattice
      scores 0/5 everywhere.  Cost is negligible at this scale (0.06 s at
      m=32, dim 64) because the instance is genuine BDD, so fplll's
      enumeration radius from the Babai residual is tiny.

  E5  The gain is narrow.  min-m for >=3/5 at K1=8 is 20 for cvp vs >40 for
      kannan, but at K1 >= 12 both are >40 while m_thresh is only 6.  So the
      practical requirement diverges from the information-theoretic bound as
      soon as eff > ~0.2; the ceiling is real, exact CVP just moves it by one
      grid step in K1.

  E6  NOT a nonce-collision artifact.  Only K1*K2 = 416 nonces exist, so with
      replacement m=24 repeats one with probability 0.48 and a repeated ECDSA
      nonce leaks d with no lattice at all.  Re-run without replacement, cvp
      still scores 5/5 at m = 24/28/32 (3/5 at m=20).  The gain is real.

Run: python3 glv_hnp_phase2_thread23_projection.py
"""

import math
import random

from fpylll import IntegerMatrix, LLL, BKZ

try:
    from fpylll import CVP as _FPCVP
except ImportError:  # pragma: no cover
    _FPCVP = None

# ---------------------------------------------------------------------------
# Arithmetic helpers (verbatim from glv_hnp_phase2_lambda_threshold.py:87)
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
    for _ in range(20000):
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is not None and y != 0:
            P = (x, y)
            if ec_mul(P, n, p) is None:
                return P
    return None

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Signatures + scales (verbatim from glv_hnp_phase2_lambda_threshold.py:227)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
    return (n // k1_bound, 1, max(1, n // k2_bound), n)

def gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed=42,
                   distinct=False):
    """distinct=True draws (k1,k2) WITHOUT replacement.  Needed as a control:
    only k1_bound*k2_bound nonces exist, so for m >~ sqrt(K1*K2) the default
    sampling repeats a nonce, and a repeated ECDSA nonce leaks d algebraically
    (d = (A_j-A_i)/(B_i-B_j)) independently of any lattice."""
    rng = random.Random(seed)
    sigs = []
    attempts = 0
    used = set()
    while len(sigs) < m and attempts < 200000:
        attempts += 1
        k1 = rng.randint(0, k1_bound - 1)
        k2 = rng.randint(0, k2_bound - 1)
        if distinct:
            if (k1, k2) in used: continue
            used.add((k1, k2))
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

def d_from_k(sig, k1, k2, n, lam):
    """d = B^{-1}(k1 + lam*k2 - A) mod n, from a single signature."""
    try:
        return (k1 + lam * k2 - sig['A']) * modinv(sig['B'], n) % n
    except ValueError:
        return None

# ---------------------------------------------------------------------------
# Variant A — baseline Kannan lattice, dim 2m+2
# ---------------------------------------------------------------------------

def build_A(sigs, n, lam, k1_bound, k2_bound):
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

def planted_A(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_A(rows, sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        last = row[dim - 1]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Variant B — projected lattice (d-column deleted), dim 2m+1
#   columns: 0..m-1 = k1 slots, m..2m-1 = k2 slots, 2m = Kannan
#   NOTE 2m+2 generators of a rank-(2m+1) lattice; LLL emits one zero row.
# ---------------------------------------------------------------------------

def build_B(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = []
    for i in range(m):                                     # c_i rows
        r = [0] * dim; r[i] = n * S_K1; M.append(r)
    r = [0] * dim                                          # d row (no d slot)
    for i in range(m): r[i] = sigs[i]['B'] * S_K1
    M.append(r)
    for i in range(m):                                     # k2_i rows
        r = [0] * dim; r[i] = -lam * S_K1; r[m + i] = S_K2; M.append(r)
    r = [0] * dim                                          # Kannan row
    for i in range(m): r[i] = sigs[i]['A'] * S_K1
    r[2 * m] = S_KAN
    M.append(r)
    return M

def planted_B(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

def recover_B(rows, sigs, n, lam, k1_bound, k2_bound, d_secret):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in rows:
        last = row[2 * m]
        if abs(last) != S_KAN: continue
        sign = 1 if last > 0 else -1
        k1_0, rem1 = divmod(sign * row[0], S_K1)
        if rem1 != 0: continue
        k2_0, rem2 = divmod(sign * row[m], S_K2)
        if rem2 != 0: continue
        d_cand = d_from_k(sigs[0], k1_0, k2_0, n, lam)
        if d_cand is not None and d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Variant C — explicit CVP, dim 2m (no d column, no Kannan column)
#   columns: 0..m-1 = k1 slots, m..2m-1 = k2 slots
#   target  = (-A_i*S_K1 | 0);  closest vector x satisfies
#             x - target = (k1_i*S_K1 | k2_i*S_K2)
# ---------------------------------------------------------------------------

def build_C(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = []
    for i in range(m):
        r = [0] * dim; r[i] = n * S_K1; M.append(r)
    r = [0] * dim
    for i in range(m): r[i] = sigs[i]['B'] * S_K1
    M.append(r)
    for i in range(m):
        r = [0] * dim; r[i] = -lam * S_K1; r[m + i] = S_K2; M.append(r)
    return M

def target_C(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    t = [0] * (2 * m)
    for i in range(m):
        t[i] = -sigs[i]['A'] * S_K1
    return t

def gso_float(basis):
    """Plain float Gram-Schmidt; dims here are <= 30 and entries <= ~2^40."""
    bstar, mu = [], []
    for v in basis:
        w = [float(x) for x in v]
        row = []
        for j, bs in enumerate(bstar):
            nb = sum(y * y for y in bs)
            c = sum(a * b for a, b in zip([float(x) for x in v], bs)) / nb if nb else 0.0
            row.append(c)
            w = [a - c * b for a, b in zip(w, bs)]
        bstar.append(w); mu.append(row)
    return bstar, mu

def babai_nearest_plane(basis, target):
    bstar, _ = gso_float(basis)
    b = [float(x) for x in target]
    coeffs = [0] * len(basis)
    for i in range(len(basis) - 1, -1, -1):
        nb = sum(y * y for y in bstar[i])
        if nb == 0: continue
        c = sum(a * y for a, y in zip(b, bstar[i])) / nb
        ci = int(math.floor(c + 0.5))
        coeffs[i] = ci
        b = [a - ci * float(x) for a, x in zip(b, basis[i])]
    out = [0] * len(target)
    for i, ci in enumerate(coeffs):
        if ci:
            for j in range(len(target)):
                out[j] += ci * basis[i][j]
    return out

def solve_C(sigs, n, lam, k1_bound, k2_bound, d_secret, use_exact_cvp=True,
            babai_only=False):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = build_C(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    basis = [[A[i][j] for j in range(2 * m)] for i in range(A.nrows)]
    basis = [r for r in basis if any(r)]          # strip the rank-deficiency
    t = target_C(sigs, n, k1_bound, k2_bound)

    cands = []
    if use_exact_cvp and not babai_only and _FPCVP is not None:
        try:
            B = IntegerMatrix.from_matrix(basis)
            LLL.reduction(B)
            cands.append(list(_FPCVP.closest_vector(B, tuple(t))))
        except Exception:
            pass
    cands.append(babai_nearest_plane(basis, t))

    for x in cands:
        diff = [xi - ti for xi, ti in zip(x, t)]
        k1_0, rem1 = divmod(diff[0], S_K1)
        if rem1 != 0: continue
        k2_0, rem2 = divmod(diff[m], S_K2)
        if rem2 != 0: continue
        d_cand = d_from_k(sigs[0], k1_0, k2_0, n, lam)
        if d_cand is not None and d_cand == d_secret:
            return d_cand, diff
    return None, (cands[0] if cands else None)

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

VARIANTS = ('kannan', 'proj', 'cvp')
ALL_VARIANTS = ('kannan', 'proj', 'cvp', 'babai')

def run_once(curve, m, d_secret, k1_bound, variant, seed=42,
             use_bkz=False, bkz_beta=20, want_geom=False, distinct=False):
    """Returns dict with 'ok', and if want_geom also 'pv','sv','det','gh'."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed,
                          distinct=distinct)
    if len(sigs) < m:
        return {'ok': False, 'note': 'sig-gen-failed'}

    if variant in ('cvp', 'babai'):
        d_cand, diff = solve_C(sigs, n, lam, k1_bound, k2_bound, d_secret,
                               babai_only=(variant == 'babai'))
        out = {'ok': d_cand is not None}
        if want_geom:
            S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
            pv = [0] * (2 * m)
            for i in range(m):
                pv[i] = sigs[i]['k1'] * S_K1
                pv[m + i] = sigs[i]['k2'] * S_K2
            out['pv'] = norm(pv)
            out['sv'] = norm(diff) if diff is not None else float('nan')
            M = build_C(sigs, n, lam, k1_bound, k2_bound)
            out.update(_lattice_geometry(M, 2 * m))
        return out

    if variant == 'kannan':
        M = build_A(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 2
    else:
        M = build_B(sigs, n, lam, k1_bound, k2_bound)
        dim = 2 * m + 1

    Amat = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(Amat, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(Amat)
    rows = [[Amat[i][j] for j in range(dim)] for i in range(Amat.nrows)]
    nz = [r for r in rows if any(r)]

    if variant == 'kannan':
        d_cand = recover_A(rows, sigs, n, lam, k1_bound, k2_bound, d_secret)
        pv = planted_A(sigs, d_secret, n, k1_bound, k2_bound)
    else:
        d_cand = recover_B(rows, sigs, n, lam, k1_bound, k2_bound, d_secret)
        pv = planted_B(sigs, d_secret, n, k1_bound, k2_bound)

    out = {'ok': d_cand is not None}
    if want_geom:
        out['pv'] = norm(pv)
        out['sv'] = min(norm(r) for r in nz) if nz else float('nan')
        out.update(_lattice_geometry(M, dim))
    return out

def _lattice_geometry(M, dim):
    """log2 det (via LLL'd GS norms) and the Gaussian heuristic."""
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]
    rows = [r for r in rows if any(r)]
    bstar, _ = gso_float(rows)
    rank = 0
    logdet = 0.0
    for bs in bstar:
        nb = math.sqrt(sum(y * y for y in bs))
        if nb > 1e-6:
            rank += 1
            logdet += math.log2(nb)
    gh = math.sqrt(rank / (2 * math.pi * math.e)) * 2 ** (logdet / rank) if rank else float('nan')
    return {'rank': rank, 'log2det': logdet, 'gh': gh}

def rate(curve, m, k1_bound, variant, seeds, use_bkz=False, bkz_beta=20):
    ok = 0
    for s in seeds:
        d_secret = (s * 7919 + 13) % curve[2]
        if d_secret == 0: d_secret = 3
        if run_once(curve, m, d_secret, k1_bound, variant, seed=s,
                    use_bkz=use_bkz, bkz_beta=bkz_beta)['ok']:
            ok += 1
    return ok

def rate_timed(curve, m, k1_bound, variant, seeds, use_bkz=False, bkz_beta=20,
               distinct=False):
    """(#recovered, mean wall-seconds per seed)."""
    import time
    ok, t0 = 0, time.time()
    for s in seeds:
        d_secret = (s * 7919 + 13) % curve[2]
        if d_secret == 0: d_secret = 3
        if run_once(curve, m, d_secret, k1_bound, variant, seed=s,
                    use_bkz=use_bkz, bkz_beta=bkz_beta, distinct=distinct)['ok']:
            ok += 1
    return ok, (time.time() - t0) / len(seeds)

# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]

HIST = [
    # label,             p,    b, n,    lam,  K1, m
    ("8-bit/199",        211,  2, 199,  106,  2,  6),
    ("12-bit/2557",      2557, 2, 2659, 1755, 8,  8),
    ("12-bit/2677 FAIL", 2677, 2, 2647, 185,  8,  10),
]

def load_hist():
    hist = []
    for label, p, b, n, lam, k1, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None, f"no generator for {label}"
        assert (lam * lam + lam + 1) % n == 0, f"bad lam for {label}"
        hist.append((label, (p, b, n, lam, G), k1, m))
    return hist


def main():
  print("=" * 78)
  print("Thread 23 — quotient out the trivial vector n*S_D*e_m  (GLV-HNP Phase 2)")
  print("=" * 78)

  hist = load_hist()

  # ---------------------------------------------------------------------------
  print("\n" + "-" * 78)
  print("EXP E1: geometry per variant.  Does sv/pv rise above 1?")
  print("-" * 78)
  print("sv = shortest nonzero row after LLL (variant cvp: the CVP residual).")
  print("pv = planted vector of that variant.  gh = Gaussian heuristic.")
  print()
  hdr = f"{'curve':<18}{'lam*':>7}{'K1':>4}{'m':>4}  {'variant':<8}{'rank':>5}" \
        f"{'log2det':>10}{'||pv||':>12}{'||sv||':>12}{'sv/pv':>8}{'pv/gh':>8}  rec"
  print(hdr)
  print("-" * len(hdr))
  e1 = {}
  for label, curve, k1, m in hist:
      n = curve[2]; lam = curve[3]
      d_secret = (42 * 7919 + 13) % n
      for var in VARIANTS:
          r = run_once(curve, m, d_secret, k1, var, seed=42, want_geom=True)
          e1[(label, var)] = r
          print(f"{label:<18}{lam_star(lam,n):>7.3f}{k1:>4}{m:>4}  {var:<8}"
                f"{r.get('rank',0):>5}{r.get('log2det',float('nan')):>10.1f}"
                f"{r['pv']:>12.1f}{r['sv']:>12.1f}"
                f"{r['sv']/r['pv']:>8.3f}{r['pv']/r['gh']:>8.3f}"
                f"  {'YES' if r['ok'] else 'no'}")

  # ---------------------------------------------------------------------------
  print("\n" + "-" * 78)
  print("EXP E2: the T4 K1-wall grid, per variant.  Does the wall move?")
  print("-" * 78)
  K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]
  print(f"5 seeds per cell; entries are #recovered/5.  K1 grid = {K1_GRID}")
  print()
  e2 = {}
  for label, curve, _k1, m in hist[1:]:      # the two 12-bit curves (T4 set)
      n = curve[2]; lam = curve[3]
      print(f"\ncurve {label}   n={n}  lam*={lam_star(lam,n):.3f}  m={m}")
      head = f"  {'variant':<8}" + "".join(f"{k:>6}" for k in K1_GRID)
      print(head); print("  " + "-" * (len(head) - 2))
      for var in VARIANTS:
          cells = []
          for k1 in K1_GRID:
              cells.append(rate(curve, m, k1, var, SEEDS))
          e2[(label, var)] = cells
          print(f"  {var:<8}" + "".join(f"{c:>5}/5" if False else f"{c:>6}" for c in cells))
      # wall = largest K1 with >=3/5
      for var in VARIANTS:
          cells = e2[(label, var)]
          wall = max([K1_GRID[i] for i, c in enumerate(cells) if c >= 3] or [0])
          print(f"  wall({var}) = K1 <= {wall}")

  # ---------------------------------------------------------------------------
  print("\n" + "-" * 78)
  print("EXP E3: m-sweep at K1=8 on the lam*=0.07 curve (T4b replication)")
  print("-" * 78)
  label, curve, _k1, _m = hist[2]
  n = curve[2]
  M_GRID = [8, 12, 16, 24, 32]
  head = f"  {'variant':<8}" + "".join(f"{mm:>6}" for mm in M_GRID)
  print(f"curve {label}  n={n}  K1=8   entries #recovered/5")
  print(head); print("  " + "-" * (len(head) - 2))
  e3 = {}
  for var in VARIANTS:
      cells = [rate(curve, mm, 8, var, SEEDS) for mm in M_GRID]
      e3[var] = cells
      print(f"  {var:<8}" + "".join(f"{c:>6}" for c in cells))

  # ---------------------------------------------------------------------------
  print("\n" + "=" * 78)
  print("VERDICT")
  print("=" * 78)
  for label, curve, k1, m in hist:
      for var in VARIANTS:
          r = e1[(label, var)]
          print(f"  E1 {label:<18} {var:<8} sv/pv = {r['sv']/r['pv']:.3f}"
                f"  ({'planted IS lambda_1-scale' if r['sv']/r['pv'] >= 0.999 else 'shorter vector exists'})")
  walls = {}
  for label, _c, _k, _m in hist[1:]:
      for var in VARIANTS:
          cells = e2[(label, var)]
          walls[(label, var)] = max([K1_GRID[i] for i, c in enumerate(cells) if c >= 3] or [0])
  print()
  for (label, var), w in walls.items():
      print(f"  E2 wall {label:<18} {var:<8} K1 <= {w}")


def main_e4():
  """E4 — is the CVP gain the FORMULATION or the exact-CVP enumeration?

  E3 showed 'cvp' reaching 5/5 at m=24 on the lam*=0.07 / K1=8 curve where
  'kannan' and 'proj' stall at 1/5.  Two possible causes:
    (i)  the d-free, Kannan-free CVP formulation is genuinely better, in which
         case Babai nearest-plane alone (polynomial, no enumeration) suffices;
    (ii) fplll's exact CVP enumeration is simply stronger reduction, in which
         case the equivalent extra power given to the Kannan embedding (BKZ)
         should close the gap too.
  E4 separates them: 'babai' = nearest-plane only, exact CVP disabled;
  plus BKZ-20/BKZ-40 controls on the Kannan embedding at the same m.
  Timings are reported because (ii) would be exponential at real scale.
  """
  print("\n" + "=" * 78)
  print("E4 — formulation or enumeration?  Babai-only vs exact CVP vs BKZ-Kannan")
  print("=" * 78)
  hist = load_hist()
  label, curve, _k1, _m = hist[2]          # 12-bit/2677, lam*=0.070
  n = curve[2]
  K1 = 8
  M_GRID = [12, 16, 20, 24, 28, 32]
  print(f"curve {label}  n={n}  K1={K1}  eff=K1*sqrt(n)/n={K1*(math.isqrt(n)+1)/n:.3f}")
  print(f"entries: #recovered/5   (mean seconds per seed in parentheses)\n")

  runs = [
      ('kannan',  False, 0),
      ('kannan',  True, 20),
      ('kannan',  True, 40),
      ('proj',    False, 0),
      ('babai',   False, 0),
      ('cvp',     False, 0),
  ]
  head = f"  {'method':<14}" + "".join(f"{mm:>14}" for mm in M_GRID)
  print(head); print("  " + "-" * (len(head) - 2))
  for var, bkz, beta in runs:
      name = f"{var}-BKZ{beta}" if bkz else f"{var}-LLL"
      cells = []
      for mm in M_GRID:
          ok, t = rate_timed(curve, mm, K1, var, SEEDS, use_bkz=bkz, bkz_beta=beta)
          cells.append(f"{ok}/5 ({t:.2f}s)")
      print(f"  {name:<14}" + "".join(f"{c:>14}" for c in cells))


def main_e5():
  """E5 — the (K1, m) recovery frontier: exact-CVP vs Kannan embedding.

  E4 established that on the lam*=0.070 curve at K1=8 the exact-CVP call
  recovers d for m >= 20 while the Kannan embedding never exceeds 2/5 even
  with BKZ-40.  So "the K1 wall" is really a (K1, m) trade-off and the two
  formulations sit on different frontiers.  E5 maps both frontiers and
  compares them to the information-theoretic bound
        m_thresh = ceil( log n / log(1/eff) ),   eff = K1*K2/n
  which is the point where the planted vector stops being expected-unique.
  """
  print("\n" + "=" * 78)
  print("E5 — (K1, m) recovery frontier: exact CVP vs Kannan embedding")
  print("=" * 78)
  hist = load_hist()
  label, curve, _k1, _m = hist[2]          # 12-bit/2677, lam*=0.070
  n = curve[2]
  k2_bound = math.isqrt(n) + 1
  K1_GRID = [8, 12, 16, 24, 32]
  M_GRID = [6, 8, 12, 16, 20, 24, 32, 40]
  print(f"curve {label}  n={n}  K2={k2_bound}  5 seeds/cell")
  print(f"min-m = smallest m in {M_GRID} reaching >=3/5\n")
  head = (f"  {'K1':>4}{'eff':>8}{'m_thresh':>10}"
          f"{'min-m(kannan)':>15}{'min-m(cvp)':>12}{'t_cvp@min-m':>13}")
  print(head); print("  " + "-" * (len(head) - 2))
  for K1 in K1_GRID:
      eff = K1 * k2_bound / n
      mt = math.ceil(math.log(n) / math.log(1.0 / eff)) if eff < 1 else float('inf')
      res = {}
      for var in ('kannan', 'cvp'):
          minm, tmin = None, None
          for mm in M_GRID:
              ok, t = rate_timed(curve, mm, K1, var, SEEDS)
              if ok >= 3:
                  minm, tmin = mm, t
                  break
          res[var] = (minm, tmin)
      kn = res['kannan'][0]; cv, tcv = res['cvp']
      print(f"  {K1:>4}{eff:>8.3f}{mt:>10}"
            f"{(str(kn) if kn else '>40'):>15}{(str(cv) if cv else '>40'):>12}"
            f"{(f'{tcv:.2f}s' if tcv else '-'):>13}")


def main_e6():
  """E6 — nonce-collision control for the E4/E5 exact-CVP gain.

  Only K1*K2 nonces exist (416 at K1=8, K2=52).  With replacement, m=24
  signatures collide with probability ~1-exp(-C(24,2)/416) = 50%, and a
  repeated ECDSA nonce yields d = (A_j-A_i)/(B_i-B_j) with no lattice at all.
  So the E4 result "cvp reaches 5/5 at m>=20" must be re-run with (k1,k2)
  drawn WITHOUT replacement before it can be called a lattice result.
  """
  print("\n" + "=" * 78)
  print("E6 — distinct-nonce control (draw (k1,k2) without replacement)")
  print("=" * 78)
  hist = load_hist()
  label, curve, _k1, _m = hist[2]
  n = curve[2]
  k2_bound = math.isqrt(n) + 1
  K1 = 8
  M_GRID = [12, 16, 20, 24, 28, 32]
  N_NONCE = K1 * k2_bound
  print(f"curve {label}  n={n}  K1={K1}  K2={k2_bound}  "
        f"nonce space = {N_NONCE}")
  print("P(>=1 repeat) with replacement, per m:")
  print("   " + "  ".join(
      f"m={mm}:{1-math.exp(-mm*(mm-1)/(2*N_NONCE)):.2f}" for mm in M_GRID))
  print("\nentries: #recovered/5\n")
  head = f"  {'method':<22}" + "".join(f"{mm:>7}" for mm in M_GRID)
  print(head); print("  " + "-" * (len(head) - 2))
  for var in ('kannan', 'cvp'):
      for dis in (False, True):
          tag = f"{var}-{'distinct' if dis else 'with-repl'}"
          cells = [rate_timed(curve, mm, K1, var, SEEDS, distinct=dis)[0]
                   for mm in M_GRID]
          print(f"  {tag:<22}" + "".join(f"{c:>7}" for c in cells))


if __name__ == "__main__":
    main()
    main_e4()
    main_e5()
    main_e6()
