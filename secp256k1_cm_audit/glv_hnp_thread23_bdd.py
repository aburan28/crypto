"""
Thread 23 — reformulate the GLV-HNP Phase-2 lattice so the planted vector is
findable, i.e. remove the trivial vector  n*S_D*e_m  as the obstruction.

Background (RESEARCH_AUTOLAB_LOG.md, 2026-07-29, EXP T5):
  The Phase-2 Kannan lattice always contains  t_triv = n*S_D*e_m  (take
  n*row_m - sum_i B_i*row_i).  With S_D = 1 its norm is n, while
  ||v_planted||^2 ~ n^2*(2m/3 + 4/3) > n^2, so the planted vector is NEVER
  lambda_1 and recovery is a BDD/coset condition, not an SVP condition.

  Proposed falsifier (2026-07-29 "Next step proposal"):
    reformulate so the target IS lambda_1 (or drop the embedding and solve
    BDD directly); if the K1 wall on the lam*=0.07 curve (currently K1~4-6)
    moves outward, the reformulation is a real improvement; if the wall stays
    put, the wall is information-theoretic and Phase 2 is at its ceiling.

Arms tested here:
  A0  baseline  : Kannan embedding, S_D = 1, LLL            (2026-07-26 setup)
  A1  S_D sweep : Kannan embedding, S_D in {1,2,4,8,16,32}, LLL
  A2  Babai     : NO embedding; LLL-reduce the (2m+1)-dim lattice and run
                  Babai nearest-plane against the explicit target
  A3  Babai/BKZ : as A2 but BKZ-20 reduction before nearest-plane
  A4  exact CVP : fpylll CVP.closest_vector (enumeration) on a small subset

Key identity used throughout: since t_triv is in the lattice, the planted
class {v_planted + j*t_triv} has a *reduced representative* whose d-coordinate
is d' = d - n when d > n/2.  Only coordinate m changes, so
    ||v_planted^red||^2 = ||v_planted||^2 - (d*S_D)^2 + (d'*S_D)^2 ,
and recovery is unaffected (d' = d mod n).  Predicted crossover where the
planted class beats t_triv:
    n^2*(2m/3 + 1) + (n/2 * S_D)^2 < (n*S_D)^2   <=>   S_D > sqrt((8m+12)/9).

Run: python3 glv_hnp_thread23_bdd.py
"""

import math
import random
import sys

from fpylll import IntegerMatrix, LLL, BKZ, CVP

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (verbatim from glv_hnp_phase2_lambda_threshold.py)
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

def tonelli_shanks(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
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

def lam_star(lam, n):
    return min(lam % n, n - (lam % n)) / n

# ---------------------------------------------------------------------------
# Signatures + scales (verbatim construction, so the comparison to
# 2026-07-26 / 2026-07-29 is exact)
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound, S_D=1):
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

# ---------------------------------------------------------------------------
# A0/A1 : Kannan embedding, dim 2m+2, with an S_D knob
# ---------------------------------------------------------------------------

def build_kannan(sigs, n, lam, k1_bound, k2_bound, S_D=1):
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

def planted_kannan(sigs, d_secret, n, k1_bound, k2_bound, S_D=1, reduced=True):
    """Planted vector; if reduced, use d' = d-n when d > n/2 (same coset)."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    d_use = d_secret - n if (reduced and d_secret > n // 2) else d_secret
    v[m] = d_use * S_D
    v[2 * m + 1] = S_KAN
    return v

def recover_kannan(rows, m, n, S_KAN, S_D, d_secret):
    dim = 2 * m + 2
    for row in rows:
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
# A2/A3/A4 : no embedding, explicit BDD.  dim 2m+1.
# ---------------------------------------------------------------------------

def build_bdd(sigs, n, lam, k1_bound, k2_bound, S_D=1):
    """Lattice rows (dim 2m+1) + target t.  The closest vector w satisfies
       w - t = (k1_i*S_K1, d*S_D, k2_i*S_K2)."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound, S_D)
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m][i] = sigs[i]['B'] * S_K1
    M[m][m] = S_D
    for i in range(m):
        M[m + 1 + i][i] = -lam * S_K1
        M[m + 1 + i][m + 1 + i] = S_K2
    t = [-sigs[i]['A'] * S_K1 for i in range(m)] + [0] * (m + 1)
    return M, t

def error_vector(sigs, d_secret, n, k1_bound, k2_bound, S_D=1, reduced=True):
    m = len(sigs)
    S_K1, S_D, S_K2, _ = scales(n, k1_bound, k2_bound, S_D)
    e = [0] * (2 * m + 1)
    for i in range(m):
        e[i] = sigs[i]['k1'] * S_K1
        e[m + 1 + i] = sigs[i]['k2'] * S_K2
    d_use = d_secret - n if (reduced and d_secret > n // 2) else d_secret
    e[m] = d_use * S_D
    return e

def gram_schmidt(B):
    """Float GS.  Returns (Bstar, mu).  Entries here are <= ~1e7 so f64 is
    exact enough (checked: all recovered d values verified against truth)."""
    k = len(B)
    Bs = []
    mu = [[0.0] * k for _ in range(k)]
    for i in range(k):
        v = [float(x) for x in B[i]]
        for j in range(len(Bs)):
            nj = sum(y * y for y in Bs[j])
            if nj == 0.0:
                mu[i][j] = 0.0
                continue
            c = sum(float(B[i][t]) * Bs[j][t] for t in range(len(v))) / nj
            mu[i][j] = c
            for t in range(len(v)):
                v[t] -= c * Bs[j][t]
        Bs.append(v)
    return Bs, mu

def babai_nearest_plane(B, t):
    """Babai nearest-plane against target t.  B should be LLL/BKZ reduced."""
    Bs, _ = gram_schmidt(B)
    dim = len(t)
    b = [float(x) for x in t]
    w = [0] * dim
    for i in reversed(range(len(B))):
        nj = sum(y * y for y in Bs[i])
        if nj == 0.0:
            continue
        c = sum(b[j] * Bs[i][j] for j in range(dim)) / nj
        c = int(math.floor(c + 0.5))
        if c:
            for j in range(dim):
                b[j] -= c * B[i][j]
                w[j] += c * B[i][j]
    return w

def recover_bdd(w, t, m, n, S_D, d_secret):
    e = [w[j] - t[j] for j in range(len(t))]
    val = e[m]
    for sign in (1, -1):
        v = sign * val
        if v % S_D != 0:
            continue
        d_cand = (v // S_D) % n
        if d_cand == d_secret:
            return d_cand
    return None

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def norm(v):
    return math.sqrt(sum(float(x) * float(x) for x in v))

def gaussian_heuristic(det, dim):
    return det ** (1.0 / dim) * math.sqrt(dim / (2 * math.pi * math.e))

def det_kannan(n, m, S_K1, S_D, S_K2, S_KAN):
    return (n * S_K1) ** m * S_D * S_K2 ** m * S_KAN

def det_bdd(n, m, S_K1, S_D, S_K2):
    return (n * S_K1) ** m * S_D * S_K2 ** m

# ---------------------------------------------------------------------------
# Single trials
# ---------------------------------------------------------------------------

def trial_kannan(curve, m, d_secret, k1_bound, seed, S_D=1, use_bkz=False, beta=20):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M = build_kannan(sigs, n, lam, k1_bound, k2_bound, S_D)
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    S_K1, S_Dv, S_K2, S_KAN = scales(n, k1_bound, k2_bound, S_D)
    ok = recover_kannan(rows, m, n, S_KAN, S_Dv, d_secret) is not None
    pv = planted_kannan(sigs, d_secret, n, k1_bound, k2_bound, S_D, reduced=True)
    pn = norm(pv)
    sn = min(norm(r) for r in rows)
    triv = n * S_Dv
    gh = gaussian_heuristic(float(det_kannan(n, m, S_K1, S_Dv, S_K2, S_KAN)), dim)
    return dict(ok=ok, pn=pn, sn=sn, triv=triv, gh=gh)

def trial_bdd(curve, m, d_secret, k1_bound, seed, S_D=1, use_bkz=False, beta=20,
              exact_cvp=False):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    M, t = build_bdd(sigs, n, lam, k1_bound, k2_bound, S_D)
    dim = 2 * m + 1
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(beta))
    else:
        LLL.reduction(A)
    rows = [[A[i][j] for j in range(dim)] for i in range(dim)]
    if exact_cvp:
        w = list(CVP.closest_vector(A, tuple(t)))
    else:
        w = babai_nearest_plane(rows, t)
    S_K1, S_Dv, S_K2, _ = scales(n, k1_bound, k2_bound, S_D)
    ok = recover_bdd(w, t, m, n, S_Dv, d_secret) is not None
    e_true = error_vector(sigs, d_secret, n, k1_bound, k2_bound, S_D, reduced=True)
    en = norm(e_true)
    Bs, _ = gram_schmidt(rows)
    gs_min = min(math.sqrt(sum(y * y for y in v)) for v in Bs)
    gh = gaussian_heuristic(float(det_bdd(n, m, S_K1, S_Dv, S_K2)), dim)
    return dict(ok=ok, en=en, gs_min=gs_min, gh=gh, dist=norm([w[j] - t[j] for j in range(dim)]))

# ---------------------------------------------------------------------------
# Sweep drivers
# ---------------------------------------------------------------------------

SEEDS = [42, 1234, 9999, 555, 31337]

def rate_kannan(curve, m, k1_bound, seeds, S_D=1, use_bkz=False, beta=20):
    n = curve[2]
    wins, pns, sns, ghs, triv = 0, [], [], [], None
    for s in seeds:
        d = random.Random(s + 7777).randint(1, n - 1)
        r = trial_kannan(curve, m, d, k1_bound, s, S_D, use_bkz, beta)
        if r is None: continue
        wins += bool(r['ok']); pns.append(r['pn']); sns.append(r['sn'])
        ghs.append(r['gh']); triv = r['triv']
    return wins, len(seeds), _mean(pns), _mean(sns), triv, _mean(ghs)

def rate_bdd(curve, m, k1_bound, seeds, S_D=1, use_bkz=False, beta=20, exact_cvp=False):
    n = curve[2]
    wins, ens, gsm, ghs = 0, [], [], []
    for s in seeds:
        d = random.Random(s + 7777).randint(1, n - 1)
        r = trial_bdd(curve, m, d, k1_bound, s, S_D, use_bkz, beta, exact_cvp)
        if r is None: continue
        wins += bool(r['ok']); ens.append(r['en']); gsm.append(r['gs_min']); ghs.append(r['gh'])
    return wins, len(seeds), _mean(ens), _mean(gsm), _mean(ghs)

def _mean(xs):
    return sum(xs) / len(xs) if xs else float('nan')


def main():
    # ===========================================================================
    print("=" * 78)
    print("Thread 23 — is the trivial vector n*S_D*e_m the obstruction?")
    print("=" * 78)

    # Historical Phase-2 curves (RESEARCH_AUTOLAB_LOG.md 2026-06-15 / 07-26 / 07-29)
    HIST = [
        # label,          p,    b, n,    lam,  m
        ("8-bit/199",     211,  2, 199,  106,  6),
        ("12-bit/2557",   2557, 2, 2659, 1755, 8),
        ("12-bit/2677",   2677, 2, 2647, 185,  10),
    ]
    CURVES = []
    for label, p, b, n, lam, m in HIST:
        G = find_generator(p, b, n)
        assert G is not None and (lam * lam + lam + 1) % n == 0, label
        CURVES.append((label, (p, b, n, lam, G), m))
        print(f"  {label:14s} n={n:5d}  lam={lam:5d}  lam*={lam_star(lam,n):.4f}  m={m}")

    K1_GRID = [2, 3, 4, 6, 8, 12, 16, 24]

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("SANITY: the planted vector really is in both lattices")
    print("-" * 78)
    _lbl, _cur, _m = CURVES[1]
    _p, _b, _n, _lam, _G = _cur
    _k2b = math.isqrt(_n) + 1
    _sigs = gen_signatures(_G, 12345 % _n, _m, _n, _lam, _p, 4, _k2b, 42)
    _M, _t = build_bdd(_sigs, _n, _lam, 4, _k2b, 1)
    _e = error_vector(_sigs, 12345 % _n, _n, 4, _k2b, 1, reduced=False)
    _w = [_t[j] + _e[j] for j in range(len(_t))]
    # solve _w as an integer combination by reducing against the row-echelon structure
    _coeff_ok = True
    _d = 12345 % _n
    _comb = [0] * len(_M)
    _comb[_m] = _d
    for i in range(_m):
        _comb[_m + 1 + i] = _sigs[i]['k2']
    _acc = [0] * len(_t)
    for i, c in enumerate(_comb):
        if c:
            for j in range(len(_t)):
                _acc[j] += c * _M[i][j]
    for i in range(_m):
        q = (_acc[i] - _w[i]) // (_n * (_n // 4))
        for j in range(len(_t)):
            _acc[j] -= q * _M[i][j]
    print(f"  BDD lattice contains t+e : {_acc == _w}")
    _MK = build_kannan(_sigs, _n, _lam, 4, _k2b, 1)
    _triv_check = [0] * (2 * _m + 2)
    _triv_check[_m] = _n
    print(f"  trivial vector n*S_D*e_m norm = {_n}, "
          f"||planted^red|| = {norm(planted_kannan(_sigs, _d, _n, 4, _k2b, 1)):.0f}")

    # ===========================================================================
    print("\n" + "=" * 78)
    print("EXP U1 — A0 baseline vs A2 Babai-BDD on the T4 K1 grid")
    print("=" * 78)
    print("A0 = Kannan+LLL (2026-07-26 setup).  A2 = no embedding, LLL + Babai NP.")
    print()
    hdr = "curve            arm   " + "".join(f"K1={k:<4d}" for k in K1_GRID)
    print(hdr)
    print("-" * len(hdr))
    U1 = {}
    for label, curve, m in CURVES:
        for arm in ("A0", "A2"):
            cells = []
            for k1 in K1_GRID:
                if arm == "A0":
                    w, tot, *_ = rate_kannan(curve, m, k1, SEEDS)
                else:
                    w, tot, *_ = rate_bdd(curve, m, k1, SEEDS)
                cells.append(w)
                U1[(label, arm, k1)] = w
            print(f"{label:16s} {arm}    " + "".join(f"{c}/5     " for c in cells))
        print()

    # ===========================================================================
    print("=" * 78)
    print("EXP U2 — S_D sweep (Kannan).  Predicted planted<trivial at S_D>sqrt((8m+12)/9)")
    print("=" * 78)
    SD_GRID = [1, 2, 4, 8, 16, 32]
    for label, curve, m in CURVES:
        thresh = math.sqrt((8 * m + 12) / 9.0)
        print(f"\n{label}  (m={m}, predicted S_D crossover = {thresh:.2f})")
        print("  S_D  |  ||pv^red||/n   triv/n   pv<triv | " +
              " ".join(f"K1={k:<3d}" for k in K1_GRID))
        print("  " + "-" * 74)
        for S_D in SD_GRID:
            cells, pn_rep, triv_rep = [], None, None
            for k1 in K1_GRID:
                w, tot, pn, sn, triv, gh = rate_kannan(curve, m, k1, SEEDS, S_D=S_D)
                cells.append(w)
                if k1 == 8:
                    pn_rep, triv_rep = pn, triv
            n = curve[2]
            flag = "YES" if pn_rep < triv_rep else "no "
            print(f"  {S_D:<4d} |  {pn_rep/n:11.2f}   {triv_rep/n:6.1f}   {flag}    | " +
                  " ".join(f"{c}/5   " for c in cells))

    # ===========================================================================
    print("\n" + "=" * 78)
    print("EXP U3 — Babai variants: LLL vs BKZ-20 vs exact CVP (enumeration)")
    print("=" * 78)
    print("Focused on the historical failure curve 12-bit/2677 (lam*=0.070).")
    label, curve, m = CURVES[2]
    print(f"\n{label}, m={m}")
    print("  arm            " + " ".join(f"K1={k:<3d}" for k in K1_GRID))
    print("  " + "-" * 62)
    for arm, kw in [("A2 LLL+Babai", dict()),
                    ("A3 BKZ20+Bab", dict(use_bkz=True, beta=20)),
                    ("A4 exact CVP", dict(exact_cvp=True))]:
        cells = []
        for k1 in K1_GRID:
            try:
                w, tot, *_ = rate_bdd(curve, m, k1, SEEDS, **kw)
            except Exception as exc:   # enumeration can blow up
                cells.append("ERR")
                continue
            cells.append(f"{w}/5")
        print(f"  {arm:14s} " + " ".join(f"{c:<6s}" for c in cells))

    # ===========================================================================
    print("\n" + "=" * 78)
    print("EXP U4 — why:  Babai succeeds iff ||e|| < min_i ||b*_i|| / 2  (roughly)")
    print("=" * 78)
    print("curve            K1   ||e||/n   min||b*||/n   ratio   pred  actual")
    print("-" * 70)
    for label, curve, m in CURVES:
        n = curve[2]
        for k1 in K1_GRID:
            w, tot, en, gsm, gh = rate_bdd(curve, m, k1, SEEDS)
            ratio = en / (gsm / 2.0) if gsm else float('inf')
            pred = "win" if ratio < 1.0 else "lose"
            print(f"{label:16s} {k1:<4d} {en/n:8.2f}  {gsm/n:11.2f}  {ratio:7.2f}   "
                  f"{pred:5s} {w}/5")
        print()

    print("=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()
