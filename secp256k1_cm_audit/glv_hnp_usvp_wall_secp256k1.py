"""
GLV-HNP — Thread 23, part 4: what the Psi law says at secp256k1 size.

Parts 1-3 fitted, on 14/17/20/22-bit j=0 GLV curves with m=12 and plain LLL,

    Psi = nu * nu_hat^gamma,   gamma ~ 0.65 (flat over [0.5, 1.0])

    nu     = E||v_planted|| / GH(L)            [Kannan lattice, dim 2m+2]
    nu_hat = lambda_1(L2) / sqrt(det L2),  L2 = <(n*S_K1,0), (-lam*S_K1,S_K2)>

with a sharp transition: Psi <= 0.80 -> 100% recovery (576 trials), Psi ~ 1.1
-> 51%, Psi >= 1.46 -> 2%.  Best cut Psi <= 0.960, AUC 0.957, and the four bit
strata collapse onto one curve.

The m-independent floor from part 1 carries over:

    nu >= nu_inf = sqrt(2*pi*e/3) * sqrt(r),  r = K1*K2/n
    ==> Psi >= Psi_floor = 2.386 * sqrt(r) * nu_hat^gamma

Psi_floor does not depend on m, so it bounds what ANY number of signatures can
achieve.  This script evaluates it on secp256k1 itself.

READ THIS BEFORE QUOTING ANY NUMBER BELOW
-----------------------------------------
This is NOT an attack on ECDSA.  It is conditional on a non-standard nonce
generator k = k1 + lam*k2 (mod n) with k1 drawn from [0, K1) and k2 from
[0, sqrt(n)) -- the "GLV-side-channel" model of RESEARCH_GLV_HNP_PHASE2.md.
Real secp256k1 signers draw k uniformly from [1, n).  Additionally:
  * Psi_crit = 0.96 was fit at 14-22 bits with plain LLL.  At 256 bits an
    attacker would use BKZ, which moves Psi_crit up by an unknown amount, so
    the K1 bound below is a LOWER bound on what is exploitable, not a security
    margin.
  * The extrapolation of a 22-bit fit to 256 bits is heuristic.

Run: python3 glv_hnp_usvp_wall_secp256k1.py
"""

import math
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = open(os.path.join(_HERE, "glv_hnp_usvp_wall_residual.py")).read()
_CUT = _SRC.index('\nprint("=" * 78)\nprint("Thread 23 part 2')
_lib = {"__name__": "_t23p2b", "__file__": os.path.join(_HERE, "glv_hnp_usvp_wall_residual.py")}
exec(compile(_SRC[:_CUT], "glv_hnp_usvp_wall_residual.py(prefix)", "exec"), _lib)

nu_of = _lib["nu_of"]
nu_inf = _lib["nu_inf"]
nu_hat = _lib["nu_hat"]
lam_star = _lib["lam_star"]
scales = _lib["scales"]

GAMMA = 0.65
PSI_CRIT = 0.9603
NU_INF_CONST = math.sqrt(2 * math.pi * math.e / 3.0)

CURVES = {
    "secp256k1": (
        0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141,
        0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72,
    ),
    "secp192k1": (
        0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE26F2FC170F69466A74DEFD8D,
        None,
    ),
}

print("=" * 78)
print("Thread 23 part 4 — the Psi floor at cryptographic size")
print("=" * 78)
print(f"gamma = {GAMMA}, Psi_crit = {PSI_CRIT} "
      f"(fit at 14-22 bits, m=12, plain LLL)")


def glv_lambda(n):
    """Smaller root of x^2+x+1 = 0 mod n, via sqrt(-3)."""
    ts = _lib["_t20"]["tonelli_shanks"]
    sq = ts((n - 3) % n, n)
    if sq is None:
        return None
    inv2 = pow(2, -1, n)
    r1 = (n - 1 + sq) * inv2 % n
    if (r1 * r1 + r1 + 1) % n != 0:
        r1 = (n - 1 + (n - sq)) * inv2 % n
    assert (r1 * r1 + r1 + 1) % n == 0
    return min(r1, n - 1 - r1)


for name, (n, lam_known) in CURVES.items():
    lam = glv_lambda(n)
    if lam is None:
        print(f"\n{name}: n has no GLV eigenvalue (n !≡ 1 mod 3) — skipped")
        continue
    if lam_known is not None:
        assert lam in (lam_known, n - 1 - lam_known), \
            f"{name}: computed lambda disagrees with the standard one"
    k2 = math.isqrt(n) + 1
    print(f"\n{'=' * 78}\n{name}")
    print(f"  n      = {n}  ({n.bit_length()} bits)")
    print(f"  lambda = {lam}   lam* = {lam_star(lam, n):.6f}")
    print(f"  K2     = isqrt(n)+1 = 2^{math.log2(k2):.2f}  (generic GLV half)")

    def m_needed(k1, nh, m_hi=10 ** 7):
        """Smallest m with nu(m)*nu_hat^gamma <= PSI_CRIT (None if the floor
        itself is above PSI_CRIT).  Doubling search then bisection."""
        if nu_inf(k1, k2, n) * (nh ** GAMMA) > PSI_CRIT:
            return None
        m = 4
        while m <= m_hi and nu_of(m, n, k1, k2) * (nh ** GAMMA) > PSI_CRIT:
            m *= 2
        if m > m_hi:
            return None
        lo, hi = m // 2, m
        while lo < hi:
            mid = (lo + hi) // 2
            if nu_of(max(2, mid), n, k1, k2) * (nh ** GAMMA) <= PSI_CRIT:
                hi = mid
            else:
                lo = mid + 1
        return lo

    print(f"\n  {'log2 K1':>8} {'eff=r':>10} {'nu_inf':>8} {'nu_hat':>8} "
          f"{'Psi_floor':>10} {'floor':>9} {'m needed':>10}")
    lo_ok = None
    for lg in range(128, 100, -1):
        k1 = 1 << lg
        r = k1 * k2 / n
        ninf = nu_inf(k1, k2, n)
        nh = nu_hat(n, lam, k1, k2)
        psi = ninf * (nh ** GAMMA)
        ok = psi <= PSI_CRIT
        if ok and lo_ok is None:
            lo_ok = lg
        mn = m_needed(k1, nh)
        print(f"  {lg:>8} {r:>10.5f} {ninf:>8.4f} {nh:>8.4f} {psi:>10.4f} "
              f"{('clears' if ok else 'blocked'):>9} "
              f"{(str(mn) if mn else '-- never --'):>10}")
        if lo_ok is not None and lg < lo_ok - 5:
            break

    if lo_ok is None:
        print("\n  no K1 in the scanned range clears Psi_crit")
    else:
        bits_short = math.log2(k2) - lo_ok
        print(f"\n  => the floor first clears Psi_crit at K1 ~ 2^{lo_ok}")
        print(f"     i.e. k1 must be at least {bits_short:.2f} bits shorter "
              f"than the generic GLV half sqrt(n) ~ 2^{math.log2(k2):.1f}")
        print(f"     and this is INDEPENDENT of the number of signatures.")

print(f"""
{'=' * 78}
Reading

The floor is shallow: for secp256k1 the requirement is a bias of a couple of
bits in k1 relative to sqrt(n), not a large one.  That is the same qualitative
place the 2026-07-29 run left it (secp256k1 nu_hat = 0.664 at eff = 0.0993,
"just above the C2 ceiling"), now with an m-independent reason attached: no
quantity of signatures moves Psi below 2.386*sqrt(r)*nu_hat^{GAMMA}, so the
attack is bounded by the nonce generator's bias alone.

This says nothing about ECDSA as deployed -- see the module docstring.
""")
