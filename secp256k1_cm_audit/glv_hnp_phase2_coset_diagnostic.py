"""
Thread 23b — what does the Phase-2 attack actually recover?

Thread 23 (this same day) projected the trivial vector n*e_m out of the
Phase-2 Kannan lattice and found the K1 wall did not move.  Diagnosing why
turned up a discrepancy between two success criteria:

  HIST  (glv_hnp_phase2_20bit.py:recover_d, used for every Phase-2 number in
        the log since 2026-06-15) scans reduced rows with |last| = S_KANNAN
        and accepts  d_cand = sign*row[m] mod n.  In a real attack this is
        checked by testing d_cand*G == Q against the public key -- cheap and
        legitimate, ~dim scalar multiplications.

  STRICT additionally demands that the k1/k2 blocks of that same row be the
        true nonce shares (k1_i in [0,K1), k2_i in [0,K2)).

On the 12-bit/2557 curve at m=12 these disagree sharply:

  K1      2    3    4    6    8   12   16   24
  HIST   5/5  5/5  5/5  5/5  5/5  4/5  1/5  0/5      (= the 2026-07-29 T4 row)
  STRICT 5/5  5/5  2/5  0/5  0/5  0/5  0/5  0/5

So from K1 >= 4 the attack recovers d off rows that are NOT the planted
vector.  This script characterises those rows.

Structure being tested.  Write a reduced row w with last coordinate S_KANNAN
as w = (Kannan row) + delta*(d-row) + sum_i gamma_i*(k2-row_i) - sum_i c_i*(n*e_i).
Then, since S_D = 1,

    w[m]     = delta                      -> HIST succeeds iff delta = d (mod n)
    w[m+1+i] = gamma_i * S_K2
    w[i]     = (k1_i + lam*(k2_i - gamma_i) - c_i*n) * S_K1

gamma_i = k2_i gives the planted vector.  But ANY gamma_i with
lam*(k2_i - gamma_i) mod n small also gives a short row with the SAME
delta = d.  Those shifts live in the 2-D lattice <(n,0), (-lam,1)>, whose
minimum is the mu of the (falsified-as-a-ratio) 2026-07-29 T2 experiment.

Prediction to check: at K1 = 8 the successful rows have delta = d exactly,
gamma != k2 (nonces NOT recovered), and gamma - k2 a short vector of the 2-D
lam-lattice.

Run: python3 glv_hnp_phase2_coset_diagnostic.py
"""

import math
import random

from fpylll import IntegerMatrix, LLL

# --- shared machinery, imported from the Thread 23 script -------------------
_src = open(__file__.replace('coset_diagnostic', 'projected')).read()
_head = _src.split('# ' + '=' * 75)[0]
_G = {}
exec(compile(_head, 'glv_hnp_phase2_projected.py', 'exec'), _G)

find_generator = _G['find_generator']
gen_signatures = _G['gen_signatures']
scales = _G['scales']
build_base_lattice = _G['build_base_lattice']
lam_star = _G['lam_star']
modinv = _G['modinv']
norm = _G['norm']


def gauss_reduce_2d(u, v):
    """Lagrange-Gauss reduction; returns the shorter basis vector."""
    def n2(w):
        return w[0] * w[0] + w[1] * w[1]
    def dot(w, z):
        return w[0] * z[0] + w[1] * z[1]
    while True:
        if n2(u) > n2(v):
            u, v = v, u
        if n2(u) == 0:
            return v
        q = round(dot(u, v) / n2(u))
        w = [v[0] - q * u[0], v[1] - q * u[1]]
        if n2(w) >= n2(v):
            return u
        v = w


def analyse(curve, m, k1_bound, seed):
    """Return a record describing the row that HIST recovers d from."""
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    d = random.Random(seed + 7777).randint(1, n - 1)
    sigs = gen_signatures(G, d, m, n, lam, p, k1_bound, k2_bound, seed)
    if len(sigs) < m:
        return None
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound, k2_bound)
    M = build_base_lattice(sigs, n, lam, k1_bound, k2_bound)
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    dim = 2 * m + 2
    red = [[A[i][j] for j in range(dim)] for i in range(A.nrows)]

    k2_true = [s['k2'] for s in sigs]
    k1_true = [s['k1'] for s in sigs]

    n_kannan = 0
    hit = None
    for row in red:
        if abs(row[dim - 1]) != S_KAN:
            continue
        n_kannan += 1
        r = row if row[dim - 1] > 0 else [-x for x in row]
        if r[m] % n != d:
            continue
        gamma = [r[m + 1 + i] // S_K2 for i in range(m)]
        shift = [gamma[i] - k2_true[i] for i in range(m)]
        # is the k1-block consistent with the coset formula?
        consistent = all(
            (r[i] // S_K1 - k1_true[i] - lam * (k2_true[i] - gamma[i])) % n == 0
            for i in range(m))
        hit = {
            'delta_eq_d': True,
            'gamma_is_k2': shift == [0] * m,
            'max_abs_shift': max(abs(x) for x in shift),
            'nz_shift': sum(1 for x in shift if x),
            'k1_in_range': all(0 <= r[i] // S_K1 < k1_bound for i in range(m)),
            'coset_formula_holds': consistent,
            'row_norm_over_planted': None,
            'shift': shift,
        }
        # planted-vector norm for reference
        pv = [0] * dim
        for i in range(m):
            pv[i] = k1_true[i] * S_K1
            pv[m + 1 + i] = k2_true[i] * S_K2
        pv[m] = d * S_D
        pv[dim - 1] = S_KAN
        hit['row_norm_over_planted'] = norm(r) / norm(pv)
        break
    return {'d': d, 'n_kannan_rows': n_kannan, 'hit': hit, 'lam': lam, 'n': n}


print("=" * 78)
print("Thread 23b — the Phase-2 attack recovers d WITHOUT recovering the nonces")
print("=" * 78)

SEEDS = [42, 1234, 9999, 555, 31337]
HIST = [
    ("12-bit/2557", 2557, 2, 2659, 1755),
    ("12-bit/2677", 2677, 2, 2647, 185),
]

for label, p, b, n, lam in HIST:
    G = find_generator(p, b, n)
    curve = (p, b, n, lam, G)
    k2_bound = math.isqrt(n) + 1
    S_K1, _S_D, S_K2, _S = scales(n, 8, k2_bound)
    mu = gauss_reduce_2d([n, 0], [-lam, 1])
    print(f"\n{label}: n={n}, lam={lam}, lam*={lam_star(lam,n):.4f}, K2={k2_bound}")
    print(f"  2-D lam-lattice <(n,0),(-lam,1)> shortest vector: "
          f"({mu[0]}, {mu[1]}), norm {math.hypot(*mu):.1f}")
    print(f"  {'K1':>3} {'seed':>6} {'#Kan':>5} {'d found':>8} {'gamma==k2':>10} "
          f"{'max|shift|':>11} {'nz':>3} {'k1 in range':>12} {'coset ok':>9} {'|w|/|pv|':>9}")
    for k1_bound in (2, 4, 8, 12):
        for seed in SEEDS[:3]:
            rec = analyse(curve, 12, k1_bound, seed)
            h = rec['hit']
            if h is None:
                print(f"  {k1_bound:>3} {seed:>6} {rec['n_kannan_rows']:>5} "
                      f"{'no':>8}")
                continue
            print(f"  {k1_bound:>3} {seed:>6} {rec['n_kannan_rows']:>5} "
                  f"{'yes':>8} {str(h['gamma_is_k2']):>10} "
                  f"{h['max_abs_shift']:>11} {h['nz_shift']:>3} "
                  f"{str(h['k1_in_range']):>12} {str(h['coset_formula_holds']):>9} "
                  f"{h['row_norm_over_planted']:>9.4f}")

print("\n" + "=" * 78)
print("done")
print("=" * 78)
