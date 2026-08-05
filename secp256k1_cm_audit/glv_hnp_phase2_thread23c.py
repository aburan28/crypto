"""
GLV-HNP Thread 23, part C — list-CVP and the rank of the true solution.

Parts A/B established:
  * every lattice point u in L2 encodes a candidate secret d'(u);
  * exact CVP recovers d iff argmin_u ||u-w|| has d'(u) = d;
  * "true error is the closest point" was sufficient in 89/89 trials.

But U4 also showed the baseline V0 (Kannan, dim 2m+2) sometimes BEATS exact
CVP at large K1 — because `recover_V0` scans all 2m+2 rows of the reduced
basis, i.e. it is already a crude list-decoder, while exact CVP returns one
point.  The principled version of that is list-CVP: enumerate the K closest
lattice points to w and test each candidate d'.  Testing is free in a real
attack (check d'*G == Q), so list size K is the honest cost parameter.

V3  list-CVP: fpylll Enumeration with BEST_N_SOLUTIONS, radius seeded at the
    Babai distance (a valid upper bound on dist(w,L2), so the K best points
    inside it really are the K closest).

Headline metric: RANK = 1-based index of the first enumerated point whose
d'(u) equals the true d.  rank=1 means exact CVP already wins; rank<=K means
list size K wins; rank=None means d is not in the K closest at all.

Run: python3 glv_hnp_phase2_thread23c.py
"""

import math
import random
import statistics

from fpylll import IntegerMatrix, LLL, GSO, Enumeration, EvaluatorStrategy, CVP

from glv_hnp_phase2_thread23 import (
    find_generator, modinv, lam_star, norm, scales, gen_signatures,
    build_V0, planted_V0, recover_V0,
    build_V2, target_V2, true_error_V2, d_from_error,
)

MAXK = 512


def expected_err2(m, k1_bound, k2_bound, S_K1, S_K2):
    """E[||e_true||^2] from PUBLIC data only (k1~U[0,K1), k2~U[0,K2)).
    Used as an enumeration radius; uses no knowledge of the secret."""
    return m * (k1_bound * S_K1) ** 2 / 3.0 + m * (k2_bound * S_K2) ** 2 / 3.0


def list_cvp_rank(sigs, n, lam, k1_bound, k2_bound, d_secret, maxk=MAXK,
                  radius_mode='public', slack=1.35):
    """Enumerate the <=maxk closest lattice points to w; return
    (rank_of_true_d, n_solutions, dist_true, dist_best).

    radius_mode:
      'babai'  -- radius = Babai distance (only finds points at least as close
                  as Babai's; can EXCLUDE the true point).
      'public' -- radius = max(Babai, slack^2 * E[||e_true||^2]) with the
                  expectation computed from the public bounds K1, K2 only.
    """
    m = len(sigs)
    S_K1, _S_D, S_K2, _S_KAN = scales(n, k1_bound, k2_bound)
    M2 = build_V2(sigs, n, lam, k1_bound, k2_bound)
    w = target_V2(sigs, n, k1_bound, k2_bound)
    e_true = true_error_V2(sigs, n, k1_bound, k2_bound)
    d_true_dist = norm(e_true)

    A = IntegerMatrix.from_matrix(M2)
    LLL.reduction(A)
    G = GSO.Mat(A)
    G.update_gso()

    # Babai gives a real lattice point => a valid enumeration radius.
    u_bab = CVP.babai(A, tuple(w))
    r0 = sum((u_bab[i] - w[i]) ** 2 for i in range(2 * m))
    r0 = float(r0) * 1.0001 + 1.0
    if radius_mode == 'public':
        r0 = max(r0, slack ** 2 * expected_err2(m, k1_bound, k2_bound, S_K1, S_K2))

    tc = G.from_canonical(tuple(w))
    enum = Enumeration(G, nr_solutions=maxk,
                       strategy=EvaluatorStrategy.BEST_N_SOLUTIONS)
    try:
        sols = enum.enumerate(0, G.d, r0, 0, target=tc)
    except Exception:
        return None, 0, d_true_dist, float('nan')

    sols = sorted(sols, key=lambda s: s[0])
    rank = None
    best = math.sqrt(sols[0][0]) if sols else float('nan')
    for idx, (dist, coeffs) in enumerate(sols, start=1):
        v = [0] * (2 * m)
        for i, ci in enumerate(coeffs):
            ci = int(round(ci))
            if not ci:
                continue
            for j in range(2 * m):
                v[j] += ci * A[i][j]
        err = [v[j] - w[j] for j in range(2 * m)]
        dc = d_from_error(err, sigs, m, n, lam, S_K1, S_K2)
        if dc == d_secret:
            rank = idx
            break
    return rank, len(sols), d_true_dist, best


def run_grid(label, curve, m, k1_grid, seeds, list_sizes=(1, 8, 64, 512)):
    p, b, n, lam, G = curve
    k2_bound = math.isqrt(n) + 1
    S_K1, S_D, S_K2, S_KAN = scales(n, k1_bound=1, k2_bound=k2_bound)  # placeholder
    print(f"\n{label}: n={n} ({n.bit_length()}b), lam*={lam_star(lam,n):.4f}, m={m}")
    print(f"  {'solver':<14}" + "".join(f"{('K1=%d' % k):>8}" for k in k1_grid))
    v0_row, rank_rows = [], {K: [] for K in list_sizes}
    med_rank = []
    for k1b in k1_grid:
        ranks, v0w, tot = [], 0, 0
        kb2 = math.isqrt(n) + 1
        for s in seeds:
            d_secret = random.Random(s).randrange(1, n)
            sigs = gen_signatures(G, d_secret, m, n, lam, p, k1b, kb2, s)
            if len(sigs) < m:
                continue
            tot += 1
            # V0 baseline
            M0 = build_V0(sigs, n, lam, k1b, kb2)
            A0 = IntegerMatrix.from_matrix(M0)
            LLL.reduction(A0)
            rows0 = [[A0[i][j] for j in range(A0.ncols)] for i in range(A0.nrows)]
            _sk1, _sd, _sk2, skan = scales(n, k1b, kb2)
            if recover_V0(rows0, m, n, skan, d_secret) is not None:
                v0w += 1
            r, _ns, _dt, _db = list_cvp_rank(sigs, n, lam, k1b, kb2, d_secret)
            ranks.append(r)
        v0_row.append(f"{v0w}/{tot}")
        for K in list_sizes:
            hit = sum(1 for r in ranks if r is not None and r <= K)
            rank_rows[K].append(f"{hit}/{tot}")
        got = [r for r in ranks if r is not None]
        med_rank.append(f"{int(statistics.median(got))}" if got else "-")
    print(f"  {'V0 (Kannan)':<14}" + "".join(f"{v:>8}" for v in v0_row))
    for K in list_sizes:
        nm = f"V3 list K={K}"
        print(f"  {nm:<14}" + "".join(f"{v:>8}" for v in rank_rows[K]))
    print(f"  {'med rank(hit)':<14}" + "".join(f"{v:>8}" for v in med_rank))


# ===========================================================================
if __name__ == '__main__':
    print("=" * 78)
    print("Thread 23c — list-CVP: how far past the wall does list-decoding reach?")
    print("=" * 78)

    SEEDS = [42, 1234, 9999, 555, 31337, 7, 20260805, 88, 4242, 123456]
    FINE = [4, 5, 6, 7, 8, 9, 10, 12, 16]

    HIST = [
        ("12-bit/2557", 2557, 2, 2659, 1755, 8),
        ("12-bit/2677", 2677, 2, 2647, 185, 10),
    ]
    for label, p, b, n, lam, m in HIST:
        G = find_generator(p, b, n)
        assert (lam * lam + lam + 1) % n == 0
        run_grid(label, (p, b, n, lam, G), m, FINE, SEEDS)

    print("\n" + "-" * 78)
    print("EXP U8: rank distribution just past the V0 wall (hard curve, K1=6..10)")
    print("-" * 78)
    p, b, n, lam, m = 2677, 2, 2647, 185, 10
    G = find_generator(p, b, n)
    kb2 = math.isqrt(n) + 1
    print(f"  {'K1':>4}{'ranks (None = not in top %d)' % MAXK:>60}")
    for k1b in (6, 7, 8, 9, 10):
        out = []
        for s in SEEDS:
            d_secret = random.Random(s).randrange(1, n)
            sigs = gen_signatures(G, d_secret, m, n, lam, p, k1b, kb2, s)
            if len(sigs) < m:
                continue
            r, ns, dt, db = list_cvp_rank(sigs, n, lam, k1b, kb2, d_secret)
            out.append(str(r) if r is not None else f"-({ns})")
        print(f"  {k1b:>4}  " + " ".join(f"{v:>6}" for v in out))

    print("\n" + "-" * 78)
    print("EXP U9: enumeration radius matters — Babai-capped vs public-data radius")
    print("        (public radius uses only K1, K2, m; no secret knowledge)")
    print("-" * 78)
    for label, p, b, n, lam, m in HIST:
        G = find_generator(p, b, n)
        kb2 = math.isqrt(n) + 1
        print(f"\n{label}: n={n}, lam*={lam_star(lam,n):.4f}, m={m}")
        print(f"  {'solver':<16}" + "".join(f"{('K1=%d' % k):>8}" for k in FINE))
        keys = ('V0', 'bab1', 'bab64', 'pub1', 'pub8', 'pub64', 'medr')
        rows = {k: [] for k in keys}
        for k1b in FINE:
            cnt = {k: 0 for k in keys[:-1]}
            ranks, tot = [], 0
            for s in SEEDS:
                d_secret = random.Random(s).randrange(1, n)
                sigs = gen_signatures(G, d_secret, m, n, lam, p, k1b, kb2, s)
                if len(sigs) < m:
                    continue
                tot += 1
                M0 = build_V0(sigs, n, lam, k1b, kb2)
                A0 = IntegerMatrix.from_matrix(M0)
                LLL.reduction(A0)
                r0 = [[A0[i][j] for j in range(A0.ncols)] for i in range(A0.nrows)]
                _s1, _sd, _s2, skan = scales(n, k1b, kb2)
                if recover_V0(r0, m, n, skan, d_secret) is not None:
                    cnt['V0'] += 1
                rb, _, _, _ = list_cvp_rank(sigs, n, lam, k1b, kb2, d_secret,
                                            radius_mode='babai')
                if rb is not None:
                    if rb <= 1: cnt['bab1'] += 1
                    if rb <= 64: cnt['bab64'] += 1
                rp, _, _, _ = list_cvp_rank(sigs, n, lam, k1b, kb2, d_secret,
                                            radius_mode='public')
                ranks.append(rp)
                if rp is not None:
                    if rp <= 1: cnt['pub1'] += 1
                    if rp <= 8: cnt['pub8'] += 1
                    if rp <= 64: cnt['pub64'] += 1
            for k in cnt:
                rows[k].append(f"{cnt[k]}/{tot}")
            got = [r for r in ranks if r is not None]
            rows['medr'].append(str(int(statistics.median(got))) if got else "-")
        names = {'V0': 'V0 (Kannan)', 'bab1': 'babai-rad K=1',
                 'bab64': 'babai-rad K=64', 'pub1': 'public-rad K=1',
                 'pub8': 'public-rad K=8', 'pub64': 'public-rad K=64',
                 'medr': 'med rank(hit)'}
        for k in keys:
            print(f"  {names[k]:<16}" + "".join(f"{v:>8}" for v in rows[k]))

    print("\n" + "=" * 78)
    print("done")
