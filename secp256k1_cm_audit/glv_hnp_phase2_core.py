"""
GLV-HNP Phase 2 — shared, import-safe core.

Why this file exists
--------------------
`glv_hnp_phase2_lambda_threshold.py` (Thread 20a) is a *script*: importing it
executes its whole experiment suite (~60 s).  Three downstream scripts
committed on 2026-07-29 --

    glv_hnp_phase2_mu_response.py   (Thread 20b)
    glv_hnp_phase2_nuhat_control.py (Thread 20c)
    glv_hnp_nuhat_vs_c1c2.py        (Thread 20d)

-- `exec_module` it and then bind names that the *committed* 20a does not
define (`glv_eigenvalues`, `mu_of`, `identify_twist`, `rival_sublattice_nu`,
and a `run_experiment` with a different signature).  All three therefore die
with AttributeError, so the headline nu_hat result (AUC 0.935) is not
reproducible from the repository as committed.  See RESEARCH_AUTOLAB_LOG.md,
2026-08-05 entry.

This module is the single source of truth: no top-level side effects, and it
exports both name sets.  The lattice construction is copied verbatim from
`glv_hnp_phase2_lambda_threshold.py:254` so numbers stay comparable.

Run `python3 glv_hnp_phase2_core.py` for a self-test.
"""

import math
import random

# ---------------------------------------------------------------------------
# Minimal EC arithmetic (short Weierstrass, a = 0)
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

# ---------------------------------------------------------------------------
# CM theory for j = 0 curves
# ---------------------------------------------------------------------------

def eisenstein_decompose(p):
    """(a,b) with a^2 - a*b + b^2 = p, a,b >= 0.  O(sqrt p)."""
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
    """The 6 Frobenius traces of the 6 sextic twists of j=0 over F_p."""
    return [2*a - b, -2*a + b, -(a + b), a + b, 2*b - a, a - 2*b]

def glv_roots(n):
    """Both roots of x^2 + x + 1 = 0 mod n.  Requires n = 1 mod 3."""
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
    """Symmetric size of lam: min(lam, n-lam)/n in [0, 0.5]."""
    return min(lam % n, n - (lam % n)) / n

def identify_twist(p, n, tries=400, seed=12345):
    """Sextic twist parameter b with #E(F_p): y^2 = x^3 + b equal to n."""
    rng = random.Random(seed)
    for _ in range(tries):
        b = rng.randint(1, p - 1)
        x = rng.randint(0, p - 1)
        rhs = (pow(x, 3, p) + b) % p
        y = tonelli_shanks(rhs, p)
        if y is None or y == 0:
            continue
        if ec_mul((x, y), n, p) is None:
            return b
    return None

def build_curve(p, n, seed=12345):
    """(p, b, n, None, G) for the twist of j=0/F_p with #E = n."""
    b = identify_twist(p, n, seed=seed)
    if b is None:
        return None
    G = find_generator(p, b, n, seed=seed)
    if G is None:
        return None
    return (p, b, n, None, G)

# ---------------------------------------------------------------------------
# 2D rival sublattice and nu_hat (Thread 20c invariant)
# ---------------------------------------------------------------------------

def gauss_reduce_2d(u, v):
    """Lagrange/Gauss reduction; returns the exact shortest nonzero vector."""
    def nrm2(w): return w[0] * w[0] + w[1] * w[1]
    def dot(w, z): return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v):
        u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0:
            break
        q = (2 * num + den) // (2 * den) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u):
            break
        u, v = v, u
    return u

def lambda_block_mu(n, lam, S_K1, S_K2):
    """Exact |shortest vector| of < (n*S_K1,0), (-lam*S_K1, S_K2) >."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] * w[0] + w[1] * w[1]), w

def rival_sublattice_nu(n, lam, k1_bound, k2_bound):
    """(lambda_1(L2), sqrt(det L2), nu_hat) for the non-planted 2D sublattice
    L2 = < (n*S_K1, 0), (-lam*S_K1, S_K2) >.  det L2 = n*S_K1*S_K2 is
    independent of lam, so nu_hat = lambda_1/sqrt(det) isolates lam."""
    S_K1, _, S_K2, _ = scales(n, k1_bound, k2_bound)
    l1, _ = lambda_block_mu(n, lam, S_K1, S_K2)
    root_det = math.sqrt(n * S_K1 * S_K2)
    return l1, root_det, l1 / root_det

# ---------------------------------------------------------------------------
# Scales, signatures
# ---------------------------------------------------------------------------

def scales(n, k1_bound, k2_bound):
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
# Lattice K  — the ORIGINAL Kannan-embedded Phase-2 lattice, dim 2m+2.
# Verbatim from glv_hnp_phase2_lambda_threshold.py:254.
#
#   cols:  0..m-1  k1-block (scale S_K1)
#          m       d        (scale S_D = 1)
#          m+1..2m k2-block (scale S_K2)
#          2m+1    Kannan   (scale S_KANNAN = n)
#
# Contains the trivial vector n*S_D*e_m (= n*row_d - sum_i B_i*row_i), of norm
# n, which is shorter than the planted vector for every m >= 1.  So the planted
# vector is never lambda_1 (Thread 20a, T5).
# ---------------------------------------------------------------------------

def build_glv_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, S_D, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
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
    return M

def planted_vector(sigs, d_secret, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + 1 + i] = sigs[i]['k2'] * S_K2
    v[m] = d_secret * S_D
    v[2 * m + 1] = S_KAN
    return v

def trivial_vector(n, m, k1_bound, k2_bound):
    """The parasitic lambda_1 of lattice K: n*S_D*e_m."""
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 2)
    v[m] = n * S_D
    return v

# ---------------------------------------------------------------------------
# Lattice P  — Thread 23: K projected along e_m (the trivial direction).
#
# L_K contains Z*(n*S_D*e_m) and that is exactly L_K cap ker(pi), where pi is
# the orthogonal projection deleting column m.  So pi(L_K) is a lattice of
# rank 2m+1 and det = det(L_K)/(n*S_D) = n^(m-1) * S_K1^m * S_K2^m * S_KANNAN.
#
# An explicit basis needs the d-row normalised so its leading entry is 1
# (B'_i = B_i * B_0^-1 mod n, B'_0 = 1); then {B', n*e_1, ..., n*e_(m-1)}
# is a basis of the k1-block sublattice {x : x = t*B' mod n}, of det n^(m-1).
# (n*e_0 = n*B' - sum_(i>=1) B'_i * (n*e_i) is in the span.)
#
#   rows:  0        normalised d-row  (B'_i * S_K1)
#          1..m-1   n*S_K1*e_i
#          m..2m-1  GLV rows  (-lam*S_K1 e_i , S_K2 e_(m+i))
#          2m       Kannan row  (A_i*S_K1 , S_KANNAN)
#
#   cols:  0..m-1   k1-block, m..2m-1 k2-block, 2m Kannan.
#
# d is no longer a coordinate; it is recovered from the k's via
# d = (k1_0 + lam*k2_0 - A_0) * B_0^-1 mod n.  Nothing is lost: pi has kernel
# spanned by e_m and lifts differ by n in the d-coordinate, so d mod n is
# well defined.
# ---------------------------------------------------------------------------

def build_projected_lattice(sigs, n, lam, k1_bound, k2_bound):
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KANNAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    Bp = [sigs[i]['B'] * B0inv % n for i in range(m)]
    assert Bp[0] == 1
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[0][i] = Bp[i] * S_K1
    for i in range(1, m):
        M[i][i] = n * S_K1
    for i in range(m):
        M[m + i][i] = -lam * S_K1
        M[m + i][m + i] = S_K2
    for i in range(m):
        M[2 * m][i] = sigs[i]['A'] * S_K1
    M[2 * m][2 * m] = S_KANNAN
    return M

def planted_vector_proj(sigs, n, k1_bound, k2_bound):
    m = len(sigs)
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    v = [0] * (2 * m + 1)
    for i in range(m):
        v[i] = sigs[i]['k1'] * S_K1
        v[m + i] = sigs[i]['k2'] * S_K2
    v[2 * m] = S_KAN
    return v

# ---------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------

def recover_d(M_reduced, m, n, S_KANNAN, d_secret):
    """ORACLE recovery for lattice K, verbatim from Thread 20a: a trial counts
    as a win iff some row has |last| = S_KANNAN and its d-entry equals the
    secret.  Kept for exact comparability with the 2026-07-29 numbers."""
    dim = 2 * m + 2
    for row in M_reduced:
        last = row[dim - 1]
        if abs(last) != S_KANNAN: continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * row[m]) % n
        if d_cand == 0: continue
        if d_cand == d_secret:
            return d_cand
    return None

def recover_verified_K(M_reduced, sigs, n, lam, k1_bound, k2_bound, G, Q, p):
    """ATTACKER-EXECUTABLE recovery for lattice K.

    d is an explicit coordinate of K, so any lattice vector lying in the right
    coset reveals it -- the k-block does NOT have to be the in-box planted
    decomposition.  The attacker discriminates candidates by checking
    d*G == Q against the public key, which costs one scalar multiplication."""
    m = len(sigs)
    dim = 2 * m + 2
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    for row in M_reduced:
        if abs(row[dim - 1]) != S_KAN: continue
        sgn = 1 if row[dim - 1] > 0 else -1
        d = (sgn * row[m]) % n
        if d == 0: continue
        if ec_mul(G, d, p) == Q:
            return d
    return None

def recover_verified_P(M_reduced, sigs, n, lam, k1_bound, k2_bound, G, Q, p):
    """ATTACKER-EXECUTABLE recovery for lattice P.

    d is no longer a coordinate, so it is reconstructed from signature 0:
    d = (k1_0 + lam*k2_0 - A_0) * B_0^-1 mod n.  This only needs the row to lie
    in the right coset (the congruence is mod n), not to be in-box -- so the
    criterion is the exact analogue of recover_verified_K."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    B0inv = modinv(sigs[0]['B'] % n, n)
    for row in M_reduced:
        if abs(row[dim - 1]) != S_KAN: continue
        sgn = 1 if row[dim - 1] > 0 else -1
        a, c = sgn * row[0], sgn * row[m]
        if a % S_K1 or c % S_K2: continue
        d = (a // S_K1 + lam * (c // S_K2) - sigs[0]['A']) * B0inv % n
        if d == 0: continue
        if ec_mul(G, d, p) == Q:
            return d
    return None

def _check_and_solve(k1, k2, sigs, n, lam):
    """Given candidate k1[], k2[], derive d from signature 0 and verify all m
    HNP congruences.  Returns d or None.  Uses no knowledge of the secret."""
    B0inv = modinv(sigs[0]['B'] % n, n)
    d = (k1[0] + lam * k2[0] - sigs[0]['A']) * B0inv % n
    for i in range(len(sigs)):
        if (sigs[i]['A'] + sigs[i]['B'] * d - k1[i] - lam * k2[i]) % n != 0:
            return None
    return d

def recover_d_free(M_reduced, sigs, n, lam, k1_bound, k2_bound):
    """ORACLE-FREE recovery for lattice K: read k1, k2 off the row, bounds-check
    them, derive d, and verify every congruence."""
    m = len(sigs)
    dim = 2 * m + 2
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in M_reduced:
        if abs(row[dim - 1]) != S_KAN: continue
        sgn = 1 if row[dim - 1] > 0 else -1
        k1, k2, ok = [], [], True
        for i in range(m):
            a, c = sgn * row[i], sgn * row[m + 1 + i]
            if a % S_K1 or c % S_K2: ok = False; break
            u, w = a // S_K1, c // S_K2
            if not (0 <= u < k1_bound and 0 <= w < k2_bound): ok = False; break
            k1.append(u); k2.append(w)
        if not ok: continue
        d = _check_and_solve(k1, k2, sigs, n, lam)
        if d is not None:
            return d
    return None

def recover_d_proj(M_reduced, sigs, n, lam, k1_bound, k2_bound):
    """ORACLE-FREE recovery for lattice P."""
    m = len(sigs)
    dim = 2 * m + 1
    S_K1, _, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    for row in M_reduced:
        if abs(row[dim - 1]) != S_KAN: continue
        sgn = 1 if row[dim - 1] > 0 else -1
        k1, k2, ok = [], [], True
        for i in range(m):
            a, c = sgn * row[i], sgn * row[m + i]
            if a % S_K1 or c % S_K2: ok = False; break
            u, w = a // S_K1, c // S_K2
            if not (0 <= u < k1_bound and 0 <= w < k2_bound): ok = False; break
            k1.append(u); k2.append(w)
        if not ok: continue
        d = _check_and_solve(k1, k2, sigs, n, lam)
        if d is not None:
            return d
    return None

# ---------------------------------------------------------------------------
# Experiment drivers
# ---------------------------------------------------------------------------

def _reduce(M, dim, use_bkz=False, bkz_beta=20):
    from fpylll import IntegerMatrix, LLL, BKZ
    A = IntegerMatrix.from_matrix(M)
    if use_bkz:
        BKZ.reduction(A, BKZ.Param(bkz_beta))
    else:
        LLL.reduction(A)
    return [[A[i][j] for j in range(dim)] for i in range(dim)]

def trial(p, n, lam, G, m, d_secret, k1_bound, k2_bound, seed=42,
          lattice='K', use_bkz=False, bkz_beta=20):
    """One end-to-end attack trial.

    Three success criteria, because they differ and the difference is the whole
    point of Thread 23:

      ok_oracle  lattice K only: Thread-20a's criterion (some row's d-entry
                 equals the secret).  Not attacker-executable, kept only so the
                 2026-07-29 numbers can be reproduced exactly.
      ok_ver     ATTACKER-EXECUTABLE: a candidate d that passes d*G == Q.  This
                 is the criterion to compare lattices on.
      ok_exact   the in-box planted decomposition (k1_i < K1, k2_i < K2) was
                 actually recovered, i.e. lambda_1-style success.

      pn, sn     ||planted||, ||shortest reduced row||
      kan        True iff the shortest reduced row has |last| = S_KANNAN
    """
    sigs = gen_signatures(G, d_secret, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    Q = ec_mul(G, d_secret, p)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    if lattice == 'K':
        dim = 2 * m + 2
        M = build_glv_lattice(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_vector(sigs, d_secret, n, k1_bound, k2_bound)
    elif lattice == 'P':
        dim = 2 * m + 1
        M = build_projected_lattice(sigs, n, lam, k1_bound, k2_bound)
        pv = planted_vector_proj(sigs, n, k1_bound, k2_bound)
    else:
        raise ValueError(lattice)
    red = _reduce(M, dim, use_bkz, bkz_beta)
    nz = [r for r in red if any(r)]
    sv = min(nz, key=norm) if nz else [0] * dim
    if lattice == 'K':
        ok_oracle = recover_d(red, m, n, S_KAN, d_secret) is not None
        d_ver = recover_verified_K(red, sigs, n, lam, k1_bound, k2_bound, G, Q, p)
        d_exact = recover_d_free(red, sigs, n, lam, k1_bound, k2_bound)
    else:
        ok_oracle = None
        d_ver = recover_verified_P(red, sigs, n, lam, k1_bound, k2_bound, G, Q, p)
        d_exact = recover_d_proj(red, sigs, n, lam, k1_bound, k2_bound)
    for got in (d_ver, d_exact):
        if got is not None and got != d_secret:
            raise AssertionError("recovery returned a wrong d")
    return {'ok_oracle': ok_oracle, 'ok_ver': d_ver is not None,
            'ok_exact': d_exact is not None,
            'pn': norm(pv), 'sn': norm(sv),
            'kan': abs(sv[dim - 1]) == S_KAN}

def success(p, n, lam, G, m, k1_bound, k2_bound, seeds, lattice='K',
            use_bkz=False, bkz_beta=20):
    """(wins_oracle, wins_verified, wins_exact, total, mean sv/pv,
    frac shortest-is-Kannan)."""
    wo = wv = we = tot = 0
    ratios, kans = [], []
    for s in seeds:
        d = random.Random(s + 7777).randint(1, n - 1)
        r = trial(p, n, lam, G, m, d, k1_bound, k2_bound, s, lattice,
                  use_bkz, bkz_beta)
        if r is None:
            continue
        tot += 1
        wo += bool(r['ok_oracle'])
        wv += bool(r['ok_ver'])
        we += bool(r['ok_exact'])
        ratios.append(r['sn'] / r['pn'] if r['pn'] else float('nan'))
        kans.append(bool(r['kan']))
    mean = lambda xs: sum(xs) / len(xs) if xs else float('nan')
    return wo, wv, we, tot, mean(ratios), mean([1.0 if k else 0.0 for k in kans])

# ---------------------------------------------------------------------------
# Compatibility aliases for the Thread 20b/20c/20d scripts
# ---------------------------------------------------------------------------

glv_eigenvalues = glv_roots
mu_of = lam_star

def run_experiment(p, n, lam, G, m, d_secret, k1_bound, k2_bound, seed=42,
                   lattice='K'):
    """Signature expected by mu_response / nuhat_control / nuhat_vs_c1c2.
    Returns bool under the Thread-20a oracle criterion (lattice K)."""
    r = trial(p, n, lam, G, m, d_secret, k1_bound, k2_bound, seed, lattice)
    if r is None:
        return False
    return bool(r['ok_oracle'] if lattice == 'K' else r['ok_ver'])


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print("glv_hnp_phase2_core self-test")
    p, b, n, lam = 211, 2, 199, 106
    G = find_generator(p, b, n)
    assert G is not None and ec_mul(G, n, p) is None
    assert (lam * lam + lam + 1) % n == 0
    k2 = math.isqrt(n) + 1
    for lat in ('K', 'P'):
        wo, wv, we, tot, r, kf = success(p, n, lam, G, 6, 2, k2,
                                         [42, 1234, 9999], lattice=lat)
        print(f"  lattice {lat}: oracle={wo}/{tot} verified={wv}/{tot} "
              f"exact={we}/{tot} sv/pv={r:.3f} kannan-shortest={kf:.2f}")
    l1, rd, nu = rival_sublattice_nu(n, lam, 2, k2)
    print(f"  nu_hat(199, lam=106) = {nu:.4f}")
    print("OK")
