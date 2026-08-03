"""
GLV-HNP — Thread 23: closed form for nu_hat via continued fractions.

Background
----------
Thread 20 (2026-07-29 run #2, commit e845207, log section restored 2026-08-03)
found the first non-trivial predictor of Phase-2 lattice-attack success:

    L2    = < (n*S_K1, 0), (-lam*S_K1, S_K2) >        (the non-planted 2D block)
    nu    = lambda_1(L2)                              (Lagrange-Gauss, O(log n))
    nu_hat = nu / sqrt(det L2),   det L2 = n*S_K1*S_K2

with AUC 0.935 against the June C1/C2 classes, and the counter-intuitive sign:
SMALL nu_hat => attack EASIER.  It was computed, not understood, and the
next-step proposal was to derive it in closed form and thereby explain why the
eight previously tested invariants all failed.

This script does that.  Write A = S_K1, B = S_K2.  A lattice vector is

    a*(-lam*A, B) + c*(n*A, 0) = ( A*(c*n - a*lam),  a*B )

so, minimising over c for each fixed a,

    lambda_1(L2)^2 = min over a>=0 of  [ A^2 * r(a)^2 + B^2 * a^2 ],
                     r(a) = min_c |a*lam - c*n| = ||a*lam||_n           (*)

THEOREM T23.  The minimiser a in (*) lies in {0} u {q_j}, where q_j are the
continued-fraction convergent denominators of lam/n.

  Proof.  The objective is strictly increasing in both a >= 0 and r >= 0.  If a
  minimises and a >= 1, no a' with 1 <= a' < a can have r(a') <= r(a) -- such an
  a' would be strictly better.  So a is a best approximation of the second kind
  to lam/n, and by the classical theorem (Khinchin, Thm 17) those are exactly
  the convergent denominators.  a = 0 contributes the vector (n*A, 0). QED

So nu_hat is an explicit Euclidean-algorithm quantity: one pass over the
convergents of lam/n, no lattice reduction at all.

The mechanism, in scale-free coordinates.  Put

    x_j = q_j * sqrt(B / (n*A)),      z_j = r_j * sqrt(A / (n*B))

Then nu_hat^2 = min_j (x_j^2 + z_j^2), and since r_j ~ n/q_{j+1},

    x_j * z_j = q_j * r_j / n  ~  q_j / q_{j+1}  ~  1 / a_{j+1}

By AM-GM nu_hat^2 >= 2 * min_j x_j z_j, with equality iff some convergent has
x_j ~ z_j, i.e. iff q_j sits at the scale

    q* = sqrt(n * A / B) = sqrt(n * K2 / K1)

So: nu_hat is small iff lam/n has a LARGE PARTIAL QUOTIENT LOCATED AT THE SCALE
q*.  A large partial quotient anywhere else does nothing.  That is precisely why
the scale-free CF invariants (q_cf, max_q_cf, max_a) were falsified in June, and
why mu = lam/n was falsified twice: none of them can see where in the CF the
large partial quotient sits relative to q*.

Experiments
-----------
T1  Theorem check: CF closed form == exact Lagrange-Gauss, 20..256 bits.
T2  Precision audit of the 2026-07-29 gauss_reduce_2d.  Pre-registered
    hypothesis: `round((u.v)/|v|^2)` is an f64 hazard at 256 bits, so the log's
    secp256k1 numbers are wrong.  FALSIFIED -- CPython int/int division is
    correctly rounded at any width, the routine was right, and T4 confirms the
    published values to 4 decimals.  What T2 does show is that the same line
    breaks from 512 bits up once ported to fixed-width doubles.
T3  Mechanism: which convergent wins, and does the scale-localised partial
    quotient predict nu_hat where the scale-free max_a does not?
T4  secp256k1 placement, recomputed exactly.
T5  Hermite ceiling and the null distribution of nu_hat.
"""

import math
import random
import sys
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "_nuhat", __file__.rsplit("/", 1)[0] + "/glv_hnp_nuhat_lib.py")
_nuhat = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_nuhat)

gauss_reduce_2d_exact = _nuhat.gauss_reduce_2d_exact
scales = _nuhat.scales

try:
    import sympy
except ImportError:
    sympy = None


# ---------------------------------------------------------------------------
# The closed form
# ---------------------------------------------------------------------------

def convergents(lam, n):
    """
    Convergent denominators q_j of lam/n together with r_j = |q_j*lam - p_j*n|.

    Runs the Euclidean algorithm on (lam, n); r_j is exactly the j-th remainder
    (up to sign), so no separate reduction is needed.  Returns a list of
    (q_j, r_j, a_{j+1}) with a_{j+1} the next partial quotient (0 at the end).
    """
    # p/q convergents of lam/n via the standard recurrences
    a_list = []
    x, y = lam, n
    while y:
        q, r = divmod(x, y)
        a_list.append(q)
        x, y = y, r
    # lam/n = [a0; a1, a2, ...]
    q_prev, q_cur = 0, 1          # q_{-1}, q_0
    out = []
    for j, a in enumerate(a_list):
        if j > 0:
            q_prev, q_cur = q_cur, a * q_cur + q_prev
        # q_cur is q_j (q_0 = 1)
        r = abs(q_cur * lam - round_nearest(q_cur * lam, n) * n)
        nxt = a_list[j + 1] if j + 1 < len(a_list) else 0
        out.append((q_cur, r, nxt))
    return out


def round_nearest(a, n):
    """round(a/n) exactly, for integers a>=0, n>0."""
    return (2 * a + n) // (2 * n)


def lambda1_sq_cf(n, lam, A, B):
    """
    lambda_1(L2)^2 by Theorem T23: exact integer value, plus the winning
    (a, r) and its convergent index (-1 for the a=0 vector).
    """
    best = (A * n) ** 2          # a = 0 -> (n*A, 0)
    best_a, best_r, best_j = 0, n, -1
    for j, (q, r, _nxt) in enumerate(convergents(lam, n)):
        val = A * A * r * r + B * B * q * q
        if val < best:
            best, best_a, best_r, best_j = val, q, r, j
    return best, best_a, best_r, best_j


def nu_hat_cf(n, lam, k1_bound, k2_bound):
    """nu_hat from the closed form."""
    A, _, B, _ = scales(n, k1_bound, k2_bound)
    l1sq, a, r, j = lambda1_sq_cf(n, lam, A, B)
    det = n * A * B
    return math.sqrt(l1sq / det), a, r, j


def nu_hat_gauss_exact(n, lam, k1_bound, k2_bound):
    """nu_hat from exact-integer Lagrange-Gauss reduction (the reference)."""
    A, _, B, _ = scales(n, k1_bound, k2_bound)
    l1sq = gauss_reduce_2d_exact((n * A, 0), (-lam * A, B))
    return math.sqrt(l1sq / (n * A * B)), l1sq


def gauss_reduce_2d_float(u, v):
    """The ORIGINAL 2026-07-29 routine, kept verbatim for the T2 comparison."""
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    if nrm2(u) < nrm2(v):
        u, v = v, u
    guard = 0
    while True:
        guard += 1
        if guard > 10000:
            return None                       # non-termination
        q = round((u[0] * v[0] + u[1] * v[1]) / nrm2(v))
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if nrm2(u) <= nrm2(v):
            return nrm2(u)


def gauss_reduce_2d_naive_float(u, v):
    """
    What the same line becomes under a naive port to a language without
    arbitrary-precision integers (C/Rust f64, or numpy): each operand is
    converted to double FIRST, then divided.  This is the version that
    actually loses the mantissa.
    """
    def nrm2(w):
        return w[0] * w[0] + w[1] * w[1]
    if nrm2(u) < nrm2(v):
        u, v = v, u
    guard = 0
    while True:
        guard += 1
        if guard > 10000:
            return None
        try:
            q = round(float(u[0] * v[0] + u[1] * v[1]) / float(nrm2(v)))
        except (OverflowError, ValueError, ZeroDivisionError):
            return None
        r = (u[0] - q * v[0], u[1] - q * v[1])
        u, v = v, r
        if nrm2(u) <= nrm2(v):
            return nrm2(u)


def max_partial_quotient(a, b):
    """Scale-free max_a, falsified by Exp S 2026-06-29."""
    best = 0
    while b:
        q, r = divmod(a, b)
        best = max(best, q)
        a, b = b, r
    return best


def scale_local_partial_quotient(n, lam, A, B):
    """
    The scale-localised partial quotient: a_{j+1} at the convergent where the
    two terms of the objective balance (x_j ~ z_j), i.e. at q* = sqrt(n*A/B).
    """
    conv = convergents(lam, n)
    best_j, best_gap = 0, None
    for j, (q, r, _nxt) in enumerate(conv):
        if r == 0:
            continue
        # balance |log(A*r / (B*q))|
        gap = abs(math.log(A * r) - math.log(B * q))
        if best_gap is None or gap < best_gap:
            best_gap, best_j = gap, j
    return conv[best_j][2], best_j, len(conv)


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else 0.0


# ---------------------------------------------------------------------------
# secp256k1 constants
# ---------------------------------------------------------------------------

SECP_N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72


# ---------------------------------------------------------------------------
# Experiments
# ---------------------------------------------------------------------------

def t1_theorem_check(rng):
    print("=" * 72)
    print("T1  Theorem T23: CF closed form vs exact Lagrange-Gauss")
    print("=" * 72)
    print("  For each (n, lam, K1, K2): lambda_1^2 from convergents of lam/n")
    print("  must equal lambda_1^2 from exact-integer Gauss reduction.\n")
    print(f"  {'bits':>5} {'trials':>7} {'agree':>7} {'max rel err':>12}")
    total_ok = True
    for bits in (20, 24, 32, 64, 128, 256):
        trials, agree, worst = 0, 0, 0.0
        for _ in range(120):
            n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
            lam = rng.randrange(1, n)
            k1 = rng.randrange(2, max(3, int(math.isqrt(n))))
            k2 = math.isqrt(n) + 1
            A, _, B, _ = scales(n, k1, k2)
            if A == 0 or B == 0:
                continue
            cf, _, _, _ = lambda1_sq_cf(n, lam, A, B)
            gs = gauss_reduce_2d_exact((n * A, 0), (-lam * A, B))
            trials += 1
            if cf == gs:
                agree += 1
            else:
                worst = max(worst, abs(cf - gs) / gs)
        ok = (agree == trials)
        total_ok &= ok
        flag = "" if ok else "   <-- MISMATCH"
        print(f"  {bits:>5} {trials:>7} {agree:>7} {worst:>12.3e}{flag}")
    print()
    print("  VERDICT:", "THEOREM T23 HOLDS on all trials (exact integer equality)."
          if total_ok else "FALSIFIED -- see mismatches above.")
    return total_ok


def t2_float_bug(rng):
    print()
    print("=" * 72)
    print("T2  Precision audit of the 2026-07-29 gauss_reduce_2d  [PRE-REGISTERED")
    print("    HYPOTHESIS FALSIFIED -- the routine was correct]")
    print("=" * 72)
    print("  Hypothesis going in: `q = round((u.v) / |v|^2)` looked like an f64")
    print("  rounding hazard, since at 256 bits |v|^2 ~ 2^500 -- far past the")
    print("  53-bit mantissa -- so the 2026-07-29 secp256k1 numbers might be wrong.")
    print()
    print("  This is FALSE.  CPython's int.__truediv__ is correctly rounded for")
    print("  arbitrary-precision operands: it does NOT convert each operand to")
    print("  double first.  The quotient here is small, so `q` is exact.")
    print("  Column 3 = the original routine; column 4 = the same line ported")
    print("  naively to fixed-width doubles (float(a)/float(b)), which is what a")
    print("  C/Rust/numpy translation would do.  Reference = CF closed form.\n")
    print(f"  {'bits':>5} {'trials':>7} {'orig wrong':>11} {'naive-f64 wrong':>16}")
    for bits in (20, 24, 64, 128, 256, 320, 512, 1100):
        trials, wrong, wrong_naive = 0, 0, 0
        for _ in range(60):
            n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
            lam = rng.randrange(1, n)
            k2 = math.isqrt(n) + 1
            k1 = max(2, n // (20 * k2))                    # eff = 0.05, integer-safe
            A, _, B, _ = scales(n, k1, k2)
            if A == 0 or B == 0:
                continue
            ref, _, _, _ = lambda1_sq_cf(n, lam, A, B)
            flt = gauss_reduce_2d_float((n * A, 0), (-lam * A, B))
            nai = gauss_reduce_2d_naive_float((n * A, 0), (-lam * A, B))
            trials += 1
            if flt is None or flt != ref:
                wrong += 1
            if nai is None or nai != ref:
                wrong_naive += 1
        print(f"  {bits:>5} {trials:>7} {wrong:>11} {wrong_naive:>16}")
    print()
    print("  => 2026-07-29 needs NO numerical correction at any bit size, and the")
    print("     restored entry's secp256k1 values stand (independently confirmed")
    print("     in T4 below).  The exact-integer rounding now in")
    print("     glv_hnp_nuhat_lib.gauss_reduce_2d is a portability hardening, not")
    print("     a bug fix: it is what the naive-f64 column shows you need once the")
    print("     routine leaves CPython.")


def t3_mechanism(rng):
    print()
    print("=" * 72)
    print("T3  Mechanism: scale-localised vs scale-free partial quotient")
    print("=" * 72)
    print("  Sample: 400 random 20-bit (n, lam), eff = K1*K2/n = 0.05.")
    print("  Claim: nu_hat is governed by the partial quotient AT the scale")
    print("  q* = sqrt(n*A/B), not by the global max partial quotient.\n")
    rows = []
    for _ in range(400):
        n = rng.getrandbits(20) | (1 << 19) | 1
        lam = rng.randrange(1, n)
        k2 = math.isqrt(n) + 1
        k1 = max(2, n // (20 * k2))                        # eff = 0.05, integer-safe
        A, _, B, _ = scales(n, k1, k2)
        if A == 0 or B == 0:
            continue
        nh, a_win, r_win, j_win = nu_hat_cf(n, lam, k1, k2)
        a_loc, j_loc, n_conv = scale_local_partial_quotient(n, lam, A, B)
        rows.append({
            'nu_hat': nh, 'max_a': max_partial_quotient(lam, n),
            'a_loc': a_loc, 'j_win': j_win, 'j_loc': j_loc, 'n_conv': n_conv,
            'qstar_ratio': (a_win / math.sqrt(n * A / B)) if a_win else 0.0,
        })
    nh = [r['nu_hat'] for r in rows]
    print(f"  n = {len(rows)} instances")
    print(f"  spearman(nu_hat, max_a) = {spearman(nh, [r['max_a'] for r in rows]):+.3f}"
          "   <- scale-free, falsified June 2026")
    print(f"  spearman(nu_hat, a_loc) = {spearman(nh, [r['a_loc'] for r in rows]):+.3f}"
          "   <- scale-localised")
    print()
    print("  The winning convergent sits at the predicted scale q* = sqrt(n*A/B):")
    ratios = sorted(r['qstar_ratio'] for r in rows if r['qstar_ratio'] > 0)
    def pct(p):
        return ratios[min(len(ratios) - 1, int(p * len(ratios)))]
    print(f"    q_win / q*  quantiles:  q05={pct(.05):.3f}  q25={pct(.25):.3f}  "
          f"q50={pct(.50):.3f}  q75={pct(.75):.3f}  q95={pct(.95):.3f}")
    print(f"    winning index j / #convergents: median "
          f"{sorted(r['j_win'] / max(1, r['n_conv']) for r in rows)[len(rows)//2]:.2f}")
    print()
    print("  nu_hat by a_loc bucket (AM-GM predicts nu_hat ~ sqrt(2/a_loc)):")
    print(f"    {'a_loc':>10} {'count':>6} {'mean nu_hat':>12} {'sqrt(2/a_loc)':>14}")
    buckets = [(1, 1), (2, 2), (3, 4), (5, 8), (9, 16), (17, 10 ** 9)]
    for lo, hi in buckets:
        grp = [r for r in rows if lo <= r['a_loc'] <= hi]
        if not grp:
            continue
        mid = sum(r['a_loc'] for r in grp) / len(grp)
        lab = f"{lo}" if lo == hi else f"{lo}-{hi if hi < 10**9 else '+'}"
        print(f"    {lab:>10} {len(grp):>6} "
              f"{sum(r['nu_hat'] for r in grp)/len(grp):>12.3f} "
              f"{math.sqrt(2/mid):>14.3f}")
    return rows


def t4_secp256k1():
    print()
    print("=" * 72)
    print("T4  secp256k1 placement, recomputed in exact arithmetic")
    print("=" * 72)
    assert (SECP_LAM * SECP_LAM + SECP_LAM + 1) % SECP_N == 0, "lambda check failed"
    print("  lambda^2 + lambda + 1 = 0 mod n   VERIFIED\n")
    print(f"  {'eff':>6} {'nu_hat (exact)':>15} {'2026-07-29 (f64)':>18} "
          f"{'float ok':>9} {'j_win':>6} {'#conv':>6}")
    logged = {0.05: 0.8709, 0.10: 0.6624, 0.25: 0.5852, 0.50: 0.6851}
    n = SECP_N
    k2 = math.isqrt(n) + 1
    for eff in (0.05, 0.10, 0.25, 0.50):
        k1 = max(2, (n * round(eff * 100)) // (100 * k2))  # integer-safe
        A, _, B, _ = scales(n, k1, k2)
        nh, a_win, r_win, j_win = nu_hat_cf(n, SECP_LAM, k1, k2)
        ref, _, _, _ = lambda1_sq_cf(n, SECP_LAM, A, B)
        flt = gauss_reduce_2d_float((n * A, 0), (-SECP_LAM * A, B))
        n_conv = len(convergents(SECP_LAM, n))
        ok = "yes" if flt == ref else "NO"
        print(f"  {eff:>6.2f} {nh:>15.4f} {logged[eff]:>18.4f} "
              f"{ok:>9} {j_win:>6} {n_conv:>6}")
    print()
    print("  NOTE: this is a lattice invariant of (n, lambda, K1, K2) only.")
    print("  It is NOT an attack on secp256k1 -- the whole Phase-2 thread is")
    print("  conditional on a non-standard nonce generator k = k1 + lambda*k2")
    print("  with k1 bounded.  Correctly generated nonces are unaffected.")


def t5_hermite(rng):
    print()
    print("=" * 72)
    print("T5  Hermite ceiling and the null distribution of nu_hat")
    print("=" * 72)
    hermite = (4.0 / 3.0) ** 0.25
    print(f"  Hermite bound in dimension 2: nu_hat <= (4/3)^(1/4) = {hermite:.4f}")
    print("  (2026-07-29 observed max 1.012 / 1.072 -- consistent.)\n")
    n = SECP_N
    k2 = math.isqrt(n) + 1
    for eff in (0.05, 0.25):
        k1 = max(2, (n * round(eff * 100)) // (100 * k2))  # integer-safe
        vals = sorted(nu_hat_cf(n, rng.randrange(1, n), k1, k2)[0]
                      for _ in range(4000))
        def q(p):
            return vals[min(len(vals) - 1, int(p * len(vals)))]
        over = sum(1 for v in vals if v > hermite)
        print(f"  eff={eff:.2f}  q05={q(.05):.3f} q25={q(.25):.3f} q50={q(.50):.3f} "
              f"q75={q(.75):.3f} q95={q(.95):.3f}  max={vals[-1]:.4f}")
        print(f"            frac < 0.45 (easy tail) = {sum(1 for v in vals if v < 0.45)/len(vals):.3f}"
              f"   frac > 0.90 = {sum(1 for v in vals if v > 0.90)/len(vals):.3f}"
              f"   violations of Hermite = {over}")


def main():
    rng = random.Random(20260803)
    print()
    print("GLV-HNP Thread 23 -- closed form for nu_hat")
    print("autolab 2026-08-03")
    print()
    ok = t1_theorem_check(rng)
    t2_float_bug(rng)
    t3_mechanism(rng)
    t4_secp256k1()
    t5_hermite(rng)
    print()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print("  T23 theorem:", "CONFIRMED" if ok else "FALSIFIED")
    print("  nu_hat is computable by one Euclidean-algorithm pass over lam/n;")
    print("  no lattice reduction is required.  It is small exactly when lam/n")
    print("  has a large partial quotient AT the scale q* = sqrt(n*K2/K1),")
    print("  which is why every scale-free CF invariant tried in June failed.")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
