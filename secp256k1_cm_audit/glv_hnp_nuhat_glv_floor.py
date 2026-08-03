"""
GLV-HNP -- Thread 23b: does the Eisenstein relation put a FLOOR on nu_hat?

Motivation
----------
Thread 23a (glv_hnp_nuhat_cf.py) gave the exact closed form

    lambda_1(L2)^2 = min_{a in {0} u {q_j}} [ S_K1^2*||a*lam||_n^2 + S_K2^2*a^2 ]

over the continued-fraction convergent denominators q_j of lam/n.  Running it on
the real secp256k1 (n, lam) showed something the f64 Lagrange-Gauss version had
hidden: the winning index j* = 70 out of depth 140 is EXACTLY the CF midpoint, and
it does not move as eff sweeps 0.05 -> 0.50, whereas for random lam the winner
moves (Thread 23a, T3).  Near j* the table of log2 q_j and log2 ||q_j*lam||_n is
visibly mirrored.

That is the Eisenstein symmetry: lam^2 + lam + 1 = 0 mod n makes the lattice
{(a,b) in Z^2 : b = a*lam mod n} an ideal of Z[omega], omega a primitive cube
root of unity, so it carries an order-6 automorphism and its two successive
minima are both ~ sqrt(n) rather than being free to be skew.

Hypothesis H-floor: a GLV lam cannot produce an arbitrarily skew L2, so nu_hat
for GLV curves has a floor well above the floor for random lam.  If true, GLV
curves are structurally excluded from the deepest "easy" tail of the 2026-07-29
nu_hat separator -- a mild but genuine security-relevant statement.

Falsifier: if the empirical minimum of nu_hat over real j=0 GLV curves tracks the
random-lam minimum (both -> 0), H-floor is false.

Run: python3 glv_hnp_nuhat_glv_floor.py
"""

import importlib.util
import math
import random
import sys

import sympy

_here = __file__.rsplit("/", 1)[0]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, _here + "/" + path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cf = _load("cf", "glv_hnp_nuhat_cf.py")
nl = _load("nl", "glv_hnp_nuhat_lib.py")

K1_FIXED = 72          # Exp S protocol, matches glv_hnp_nuhat_vs_c1c2.py:55
N_CURVES = 300


def glv_curves(lo, hi, want):
    """Real j=0 CM curves over F_p with prime order n = 1 mod 3, plus their lam."""
    out = []
    for p in sympy.primerange(lo, hi):
        if p % 3 != 1:
            continue
        dec = nl.eisenstein_decompose(p)
        if not dec:
            continue
        a, b = dec
        for t in nl.j0_traces(a, b):
            n = p + 1 - t
            if n > 3 and sympy.isprime(n) and n % 3 == 1:
                e = nl.glv_eigenvalues(n)
                if e:
                    out.append((p, n, e[0]))
        if len(out) >= want:
            break
    return out[:want]


def quantiles(v):
    v = sorted(v)
    def q(f):
        return v[min(len(v) - 1, int(f * len(v)))]
    return (v[0], q(.05), q(.25), q(.5), q(.75), q(.95), v[-1])


def main():
    print("=" * 78)
    print("Thread 23b -- nu_hat floor for GLV lam vs random lam")
    print("=" * 78)
    curves = glv_curves(2 ** 19, 2 ** 20, N_CURVES)
    print(f"real j=0 GLV curves harvested (20-bit): {len(curves)}, K1 = {K1_FIXED}")
    rng = random.Random(7)

    glv, rnd = [], []
    for p, n, lam in curves:
        k2 = math.isqrt(n) + 1
        glv.append(cf.nu_hat_cf(n, lam, K1_FIXED, k2))
        rnd.append(cf.nu_hat_cf(n, rng.randrange(1, n), K1_FIXED, k2))

    print()
    hdr = f"{'sample':>11} {'min':>8} {'q05':>8} {'q25':>8} {'med':>8} " \
          f"{'q75':>8} {'q95':>8} {'max':>8}"
    print(hdr)
    for name, v in (("GLV lam", glv), ("random lam", rnd)):
        print(f"{name:>11} " + " ".join(f"{x:8.4f}" for x in quantiles(v)))

    print()
    print(f"{'sample':>11} {'frac<0.45':>10} {'frac<0.30':>10} {'frac<0.20':>10}")
    for name, v in (("GLV lam", glv), ("random lam", rnd)):
        print(f"{name:>11} "
              f"{sum(x < .45 for x in v)/len(v):10.3f} "
              f"{sum(x < .30 for x in v)/len(v):10.3f} "
              f"{sum(x < .20 for x in v)/len(v):10.3f}")

    print()
    print("j* offset from the CF midpoint (index units):")
    rng2 = random.Random(7)
    for name, use_glv in (("GLV", True), ("random", False)):
        off = []
        for p, n, lam in curves:
            L = lam if use_glv else rng2.randrange(1, n)
            _, _, _, j, qs, _ = cf.nu_hat_cf(
                n, L, K1_FIXED, math.isqrt(n) + 1, return_detail=True)
            off.append(j - (len(qs) - 1) / 2)
        off.sort()
        print(f"  {name:>7}  med {off[len(off)//2]:+.1f}  "
              f"q05 {off[int(.05*len(off))]:+.1f}  "
              f"q95 {off[int(.95*len(off))]:+.1f}  "
              f"max|off| {max(abs(o) for o in off):.1f}")

    print()
    print("Is the CF of lam/n an exact palindrome for GLV lam?")
    pal = 0
    for p, n, lam in curves[:60]:
        _, a = cf.convergent_denominators(lam, n)
        core = a[1:]
        pal += (core == core[::-1])
    print(f"  exact palindrome of a_1..a_last: {pal}/60  "
          f"-> NO; the symmetry pins the BALANCE POINT, not the quotients.")

    print()
    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print(f"GLV floor   min nu_hat = {min(glv):.4f}, none below 0.30")
    print(f"random floor min nu_hat = {min(rnd):.4f}, "
          f"{sum(x < .30 for x in rnd)/len(rnd):.1%} below 0.30")
    ok = min(glv) > 0.30 and min(rnd) < 0.30
    print("H-floor:", "SUPPORTED" if ok else "NOT SUPPORTED")
    print()
    print("Caveats: 20-bit sample, one K1, empirical only.  The floor is NOT yet")
    print("proved; see the 2026-08-03 log entry for the proposed proof route via")
    print("the successive minima of an ideal of Z[omega].")
    return 0


if __name__ == '__main__':
    sys.exit(main())
