"""
Thread 23 follow-up: does the centring defect found in the GLV Phase-2 lattice
also cost one bit in the STANDARD (non-GLV) Boneh-Venkatesan lattice used by
`src/cryptanalysis/hnp_ecdsa.rs` (and exercised by
`tests/curve_audit.rs::hnp_recovers_p256_key_via_lll_phase15`)?

The Rust basis (hnp_ecdsa.rs:216-229, doc-comment at hnp_ecdsa.rs:30-40) is

    row i (i<m) : n^2 * e_i
    row m       : (n*a_0, ..., n*a_{m-1},  2^kb,  0        )
    row m+1     : (n*t_0, ..., n*t_{m-1},  0,     n*2^kb   )

with target w = (n*k_0, ..., n*k_{m-1}, d*2^kb, n*2^kb) and k_i in [0, 2^kb).
The k_i are ONE-SIDED, so every one of the m informative coordinates carries a
mean offset of n*2^kb/2:  E[coord^2] = (n*2^kb)^2/3 instead of /12.

Centring fix (one line): put  n*(t_i - 2^(kb-1))  in row m+1 instead of n*t_i.
Then the target coordinate is n*(k_i - 2^(kb-1)), zero-mean.

The HNP instance depends only on (a_i, t_i) mod n, so no EC arithmetic is
needed: draw a_i uniform, d uniform, k_i uniform in [0,2^kb), set
t_i = (k_i - a_i*d) mod n.  This is exactly the distribution the Rust code sees.

Run: python3 hnp_centring_256bit.py
"""

import random

from fpylll import IntegerMatrix, LLL, BKZ

# secp256k1 group order (matches src/ecc curve params)
N_SECP256K1 = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
# NIST P-256 group order (the curve the Rust regression test uses)
N_P256 = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551


def make_instance(n, m, kb, seed):
    rng = random.Random(seed)
    d = rng.randrange(1, n)
    a, t, k = [], [], []
    for _ in range(m):
        ai = rng.randrange(1, n)
        ki = rng.randrange(0, 1 << kb)
        a.append(ai)
        t.append((ki - ai * d) % n)
        k.append(ki)
    return d, a, t, k


def build(n, a, t, kb, centred):
    m = len(a)
    dim = m + 2
    two_l = 1 << kb
    half = (1 << (kb - 1)) if centred else 0
    M = [[0] * dim for _ in range(dim)]
    for i in range(m):
        M[i][i] = n * n
    for i in range(m):
        M[m][i] = n * a[i]
    M[m][m] = two_l
    for i in range(m):
        M[m + 1][i] = n * ((t[i] - half) % n)
    M[m + 1][m + 1] = n * two_l
    return M


def recover(n, m, kb, M, d_secret):
    dim = m + 2
    A = IntegerMatrix.from_matrix(M)
    LLL.reduction(A)
    two_l = 1 << kb
    for i in range(dim):
        v = A[i]
        if v[m] % two_l:
            continue
        for sgn in (1, -1):
            cand = (sgn * v[m] // two_l) % n
            if cand == d_secret:
                return True
    return False


def sweep(n, label, m_list, kb_bits_list, seeds):
    print(f"\n{'='*74}")
    print(f"  {label}   (n is {n.bit_length()} bits)")
    print(f"{'='*74}")
    print(f"{'m':>4} {'bias bits':>10} {'m*bias/nbits':>13} "
          f"{'uncentred':>10} {'centred':>9}")
    nb = n.bit_length()
    for m in m_list:
        for bias in kb_bits_list:
            kb = nb - bias
            wu = wc = 0
            for seed in seeds:
                d, a, t, k = make_instance(n, m, kb, seed)
                if recover(n, m, kb, build(n, a, t, kb, False), d):
                    wu += 1
                if recover(n, m, kb, build(n, a, t, kb, True), d):
                    wc += 1
            print(f"{m:>4} {bias:>10} {m * bias / nb:>13.3f} "
                  f"{str(wu) + '/' + str(len(seeds)):>10} "
                  f"{str(wc) + '/' + str(len(seeds)):>9}")


if __name__ == '__main__':
    SEEDS = [1, 2, 3, 4, 5]
    print("Thread 23 follow-up — centring the standard Boneh-Venkatesan lattice")
    print("Rust reference: src/cryptanalysis/hnp_ecdsa.rs:216-229 (uncentred)")

    # The interesting region is the LLL threshold m*bias ~ 1.5*n_bits
    # (hnp_ecdsa.rs:174).  Walk bias down until each variant breaks.
    sweep(N_P256, "NIST P-256", [4, 6, 8, 12, 20],
          [80, 64, 48, 40, 32, 24, 20, 16, 12, 10, 8], SEEDS)
    sweep(N_SECP256K1, "secp256k1", [8, 20], [32, 20, 16, 12, 10, 8], SEEDS)
