"""
Thread 23, addendum: where does the K1 wall sit on the Gaussian-heuristic axis,
and does nu_hat (Thread 20c) explain the residual spread?

E3/E5 of glv_hnp_phase2_thread23.py showed the d-eliminated reformulation does
not move the wall.  This addendum asks the follow-up question: what DOES pin the
wall?  For each curve it reports

    tau  = ||v_planted|| / GH(det, dim)        at the measured wall K1
    nu_hat = lambda_1(L2) / sqrt(det L2),  L2 = <(n*S_K1,0), (-lam*S_K1, S_K2)>
             (the Thread 20c lattice-geometric invariant, log 2026-07-29)

tau is blind to lam (it depends only on n, K1, K2, m), so if the wall were purely
a tau phenomenon every curve would wall at the same tau.  Any spread in
tau-at-wall must be explained by a lam-sensitive quantity, and nu_hat is the only
one that survived Thread 20b/20c.

Run: python3 glv_hnp_phase2_thread23_tau_wall.py
"""

import math

from glv_hnp_phase2_thread23 import (
    find_generator, build_curve, glv_roots, lam_star, scales,
    det_log2_new, expected_target_new, gh_ratio, run_old, rate,
)

def gauss_reduce_2d(u, v):
    def nrm2(w): return w[0] * w[0] + w[1] * w[1]
    def dot(w, z): return w[0] * z[0] + w[1] * z[1]
    if nrm2(u) > nrm2(v): u, v = v, u
    while True:
        num, den = dot(v, u), nrm2(u)
        if den == 0: break
        q = ((2 * num + den) // (2 * den)) if num >= 0 else -((-2 * num + den) // (2 * den))
        v = (v[0] - q * u[0], v[1] - q * u[1])
        if nrm2(v) >= nrm2(u): break
        u, v = v, u
    return u

def nu_hat(n, lam, S_K1, S_K2):
    """lambda_1(L2) / sqrt(det L2),  det L2 = n*S_K1*S_K2."""
    w = gauss_reduce_2d((n * S_K1, 0), (-(lam % n) * S_K1, S_K2))
    return math.sqrt(w[0] ** 2 + w[1] ** 2) / math.sqrt(n * S_K1 * S_K2)

SEEDS = [42, 1234, 9999, 555, 31337]
M = 12
GRID = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]

CURVES = []
for label, p, b, n, lam in [
    ("12-bit/2557", 2557, 2, 2659, 1755),
    ("12-bit/2677", 2677, 2, 2647, 185),
]:
    G = find_generator(p, b, n)
    CURVES.append((label, (p, b, n, lam, G)))
for p, n in [(131113, 130981), (131149, 130483), (131203, 131707)]:
    cur = build_curve(p, n)
    lam = glv_roots(n)[0]
    CURVES.append((f"17-bit/{n}", (p, cur[1], n, lam, cur[4])))

print("=" * 84)
print("Thread 23 addendum — tau at the measured K1 wall, and nu_hat")
print("=" * 84)
print(f"m = {M}, K2 = isqrt(n)+1, seeds = {SEEDS}, old lattice, wall = largest K1 >= 4/5")
print()
print(f"{'curve':<14} {'n':>8} {'lam*':>7} {'nu_hat':>7} {'wall K1':>8} "
      f"{'tau@wall':>9} {'tau@wall+1':>11}")

for label, curve in CURVES:
    p, b, n, lam, G = curve
    K2 = math.isqrt(n) + 1
    taus, wins = {}, {}
    for K1 in GRID:
        S_K1, _, S_K2, S_KAN = scales(n, K1, K2)
        taus[K1] = gh_ratio(det_log2_new(n, M, S_K1, S_K2, S_KAN), 2 * M + 1,
                            expected_target_new(M, n, K1, K2, S_K1, S_K2, S_KAN))
        w, _t, _ = rate(run_old, curve, M, K1, SEEDS)
        wins[K1] = w
    ok = [K1 for K1 in GRID if wins[K1] >= 4]
    wall = max(ok) if ok else 0
    nxt = [K1 for K1 in GRID if K1 > wall]
    tau_w = taus[wall] if wall else float('nan')
    tau_n = taus[nxt[0]] if nxt else float('nan')
    S_K1, _, S_K2, _ = scales(n, max(wall, 2), K2)
    nh = nu_hat(n, lam, S_K1, S_K2)
    print(f"{label:<14} {n:>8} {lam_star(lam, n):>7.4f} {nh:>7.3f} {wall:>8} "
          f"{tau_w:>9.3f} {tau_n:>11.3f}")
    print(f"{'':<14} wins: " + "  ".join(f"K1={k}:{wins[k]}/5" for k in GRID))

print()
print("Reading: the wall is bracketed by tau in [tau@wall, tau@wall+1].  If every")
print("curve's bracket contains a common tau, the wall is a Gaussian-heuristic")
print("phenomenon and no det/dim-preserving reformulation can move it.")


# ---------------------------------------------------------------------------
# Addendum part 2: the large-m limit of tau, and what it says about the ceiling.
#
#   log2 det_new = (m-1)L + m*log2(S_K1) + m*log2(S_K2) + L,  L = log2 n
#                = m*(3L - log2(K1*K2))              (S_Kx = n/Kx exactly)
#   dim          = 2m+1
#   ||target||^2 = m*(K1*S_K1)^2/3 + m*(K2*S_K2)^2/3 + n^2  ~ n^2*(2m/3 + 1)
# so as m -> infinity
#   tau -> sqrt(2*pi*e/3) * sqrt(K1*K2/n) = 2.3900 * sqrt(eff),  eff = K1*K2/n.
# The wall tau = 1 therefore sits at eff = 3/(2*pi*e) = 0.17547, independent of
# n and m.  That is exactly the "eff" variable T3 (2026-07-29) found empirically.
# ---------------------------------------------------------------------------

print()
print("=" * 84)
print("Addendum part 2 — large-m limit:  tau -> sqrt(2*pi*e/3)*sqrt(eff)")
print("=" * 84)
C = math.sqrt(2 * math.pi * math.e / 3)
print(f"constant sqrt(2*pi*e/3) = {C:.4f};  tau = 1  <=>  eff = 3/(2*pi*e) = "
      f"{3 / (2 * math.pi * math.e):.5f}")

print(f"\nconvergence of exact tau to the limit (n = 2^256-ish, K2 = isqrt(n)):")
n_big = (1 << 256) - 432420386565659656852420866394968145599
K2b = math.isqrt(n_big)
print(f"{'eff':>7} {'K1/2^128':>10} | " + " ".join(f"m={m:<5}".rjust(8)
      for m in [4, 12, 48, 200, 1000]) + f" {'limit':>8}")
for eff in [0.05, 0.10, 0.17547, 0.25, 0.50, 1.00]:
    K1b = max(1, int(eff * n_big / K2b))
    cells = []
    for m in [4, 12, 48, 200, 1000]:
        S_K1, _, S_K2, S_KAN = scales(n_big, K1b, K2b)
        t = gh_ratio(det_log2_new(n_big, m, S_K1, S_K2, S_KAN), 2 * m + 1,
                     expected_target_new(m, n_big, K1b, K2b, S_K1, S_K2, S_KAN))
        cells.append(f"{t:.4f}".rjust(8))
    print(f"{eff:>7.5f} {K1b / 2**128:>10.4f} | " + " ".join(cells)
          + f" {C * math.sqrt(eff):>8.4f}")

print("\nCross-check against T3 (2026-07-29, 20 fresh 17-bit curves, m=12).")
print("The limit is only approached for m in the hundreds, so the exact finite-m")
print("tau is the one to compare against; it runs ~1.5x the limit at m=12.")
n17 = 131101
K2_17 = math.isqrt(n17) + 1
print(f"  {'eff':>6} {'tau_limit':>10} {'tau_exact(m=12)':>16} {'T3 result':>20}")
for eff, res in [(0.05, "19/20 recover 5/5"), (0.15, "3/20 recover"),
                 (0.25, "0/20 recover")]:
    K1_17 = max(2, int(round(eff * n17 / K2_17)))
    S_K1, _, S_K2, S_KAN = scales(n17, K1_17, K2_17)
    t_ex = gh_ratio(det_log2_new(n17, 12, S_K1, S_K2, S_KAN), 25,
                    expected_target_new(12, n17, K1_17, K2_17, S_K1, S_K2, S_KAN))
    print(f"  {eff:>6.2f} {C * math.sqrt(eff):>10.3f} {t_ex:>16.3f} {res:>20}")
print("  => the T3 transition brackets tau_exact = 1, same wall as E3/E5.")

print("\nCEILING for secp256k1 (n = order of secp256k1, K2 = sqrt(n) = 2^128):")
print(f"  tau < 1 requires eff = K1*K2/n < {3/(2*math.pi*math.e):.5f},")
print(f"  i.e. K1 < {3/(2*math.pi*math.e):.5f} * 2^128 = 2^{128 + math.log2(3/(2*math.pi*math.e)):.2f}.")
print(f"  Unbiased GLV gives K1 = K2 = 2^128, eff = 1.  The attack therefore needs")
print(f"  log2(1/{3/(2*math.pi*math.e):.5f}) = {-math.log2(3/(2*math.pi*math.e)):.2f} bits of combined bias in (k1,k2)")
print(f"  before the planted vector is even expected to be lambda_1 -- and LLL/BKZ")
print(f"  need more than that.  This is Phase 2's ceiling, in bits.")
