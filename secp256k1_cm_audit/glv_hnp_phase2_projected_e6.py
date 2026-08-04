"""
GLV-HNP Phase 2, Thread 23 / E6: can the spectral (Bleichenbacher) route reach
into the gap that E5 opened?

E5 showed that for the historical failure curve (n=2647, m=10, K2=52) the
private key d is UNIQUELY determined for K1 = 8, 12, 16 -- exactly one
survivor out of 2646 candidates -- while LLL and BKZ-40 recover it 0/5 times.
The information is present; the lattice cannot extract it.  That gap spans
eff = K1*K2/n in [0.157, 0.47], about a 3x range in K1.

The repo already carries `src/cryptanalysis/bleichenbacher.rs`, an FFT bias
detector that "works below the lattice threshold".  Its bias function is

    Z(d) = sum_i exp(2*pi*I * omega * (A_i + B_i*d) / n)

For the classic HNP the nonce set is an interval [0, 2^b) and the bias of the
correct d is sinc-like and large.  Here the nonce set is the GLV set

    S = { k1 + lam*k2 mod n : 0 <= k1 < K1, 0 <= k2 < K2 }

which is NOT an interval -- it is a union of K2 short intervals of length K1
scattered by lam.  Its bias factorises exactly:

    B(omega) = |G(omega, K1)| / K1  *  |G(omega*lam, K2)| / K2
    G(a, K)  = sum_{j<K} exp(2*pi*I*a*j/n)   (closed-form geometric sum)

The two factors fight each other: omega = 1 makes the k1 factor ~1 but leaves
lam generic in the k2 factor; omega = lam^{-1} does the reverse.  Whether ANY
omega makes the product large is a 2646-term search -- cheap, and it decides
whether the spectral route is viable before anyone writes the attack.

Sample complexity: Bleichenbacher needs roughly m ~ 1/B^2 signatures to lift
the peak above the sqrt(m) noise floor of the wrong candidates.

  B large  (m_needed ~ tens)   -> spectral route is viable; write the attack.
  B ~ 1/sqrt(n)                -> no better than exhaustive search; route dead.

Run: python3 glv_hnp_phase2_projected_e6.py
"""

import cmath
import math

from glv_hnp_phase2_projected import find_generator, lam_star, HIST


def geom(a, K, n):
    """|sum_{j=0}^{K-1} exp(2*pi*I*a*j/n)|, exact closed form."""
    a %= n
    if a == 0:
        return float(K)
    num = 1 - cmath.exp(2j * math.pi * a * K / n)
    den = 1 - cmath.exp(2j * math.pi * a / n)
    return abs(num / den)


def glv_bias(omega, n, lam, K1, K2):
    return (geom(omega, K1, n) / K1) * (geom(omega * lam, K2, n) / K2)


def interval_bias(omega, n, K):
    """Bias of the plain-interval nonce set [0,K), for reference."""
    return geom(omega, K, n) / K


print("=" * 78)
print("Thread 23 / E6 - spectral (Bleichenbacher) bias of the GLV nonce set")
print("=" * 78)

CURVES = [(lb, p, b, n, lam, m) for (lb, p, b, n, lam, k1, m) in HIST
          if '12-bit' in lb]

for label, p, b, n, lam, m in CURVES:
    assert find_generator(p, b, n) is not None
    K2 = math.isqrt(n) + 1
    print()
    print(f"  {label}  n={n}  lam={lam}  lam*={lam_star(lam, n):.4f}  K2={K2}")
    print(f"  {'K1':>3} {'eff':>6} | {'B(1)':>9} {'B(lam^-1)':>10} "
          f"{'max_w B':>9} {'argmax w':>9} | {'m~1/B^2':>9} {'vs interval':>12}")
    print("  " + "-" * 76)
    for K1 in (4, 8, 12, 16, 24):
        eff = K1 * K2 / n
        lam_inv = pow(lam, -1, n)
        b1 = glv_bias(1, n, lam, K1, K2)
        bl = glv_bias(lam_inv, n, lam, K1, K2)
        best_w, best_b = max(
            ((w, glv_bias(w, n, lam, K1, K2)) for w in range(1, n)),
            key=lambda t: t[1])
        m_need = (1.0 / best_b ** 2) if best_b > 0 else float('inf')
        # reference: an interval nonce set of the same cardinality K1*K2
        ib = interval_bias(1, n, min(n, K1 * K2))
        print(f"  {K1:>3} {eff:>6.3f} | {b1:>9.5f} {bl:>10.5f} "
              f"{best_b:>9.5f} {best_w:>9} | {m_need:>9.1f} {ib:>12.5f}")

print()
print("  Reference scales:")
for label, p, b, n, lam, m in CURVES:
    print(f"    {label}: 1/sqrt(n) = {1/math.sqrt(n):.5f}   "
          f"(bias of a uniform set; anything at this level is worthless)")


# ===========================================================================
# E6b - is the optimal spectral multiplier the lambda-block short vector?
# ===========================================================================
# H23b: the omega maximising B(omega) is exactly the omega-component of the
# Gauss-reduced shortest vector of
#       W = < (K1, lam*K2), (0, n*K2) >     [ = (omega*K1, (omega*lam mod n)*K2) ]
# i.e. the SAME 2-D lambda-block whose vector mu blocks LLL in E4a.
# If so, the obstruction and the alternative attack's key are the same object.

def gauss2(u, v):
    def n2(w): return w[0] * w[0] + w[1] * w[1]
    if n2(u) > n2(v): u, v = v, u
    while True:
        q = int(math.floor((u[0]*v[0] + u[1]*v[1]) / n2(u) + 0.5))
        v = (v[0] - q*u[0], v[1] - q*u[1])
        if n2(v) >= n2(u): return u
        u, v = v, u

def centered(a, n):
    a %= n
    return a - n if a > n // 2 else a

print()
print("=" * 78)
print("E6b  H23b: optimal spectral multiplier == lambda-block short vector?")
print("=" * 78)
print(f"  {'curve':<16} {'K1':>3} | {'w_pred':>7} {'w_argmax':>9} {'match':>6} | "
      f"{'w~':>6} {'(w*lam)~':>9} | {'B(w_pred)':>10} {'B(max)':>8}")
print("  " + "-" * 76)

n_match = n_tot = 0
for label, p, b, n, lam, m in CURVES:
    K2 = math.isqrt(n) + 1
    for K1 in (4, 8, 12, 16, 24):
        sv = gauss2((K1, (lam * K2) % (n * K2)), (0, n * K2))
        w_pred = abs(sv[0]) // K1
        if w_pred == 0:
            continue
        best_w, best_b = max(((w, glv_bias(w, n, lam, K1, K2))
                              for w in range(1, n)), key=lambda t: t[1])
        bp = glv_bias(w_pred, n, lam, K1, K2)
        # omega and n-omega give the same |bias|; compare up to that symmetry
        ok = (w_pred % n in (best_w % n, (n - best_w) % n)) or abs(bp - best_b) < 1e-9
        n_match += ok; n_tot += 1
        print(f"  {label:<16} {K1:>3} | {w_pred:>7} {best_w:>9} {str(ok):>6} | "
              f"{centered(w_pred, n):>6} {centered(w_pred*lam, n):>9} | "
              f"{bp:>10.5f} {best_b:>8.5f}")
print(f"  H23b confirmed in {n_match}/{n_tot} cells")

print()
print("=" * 78)
print("done")
print("=" * 78)
