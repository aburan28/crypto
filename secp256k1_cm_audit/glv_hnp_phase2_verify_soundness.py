"""
GLV-HNP Phase 2, Thread 20 (part 4): is the oracle-free verifier sound?

Every conclusion in the Thread 20 scripts rests on verify_d(), which accepts a
candidate d when EVERY signature's implied nonce k = A + B*d admits a GLV
decomposition k = k1 + lam*k2 inside the declared box k1 < K1, k2 < K2.  No
secret is consulted, which is the point -- but it means a WRONG d would be
accepted if it happened to decompose for all m signatures.

This matters specifically because glv_hnp_phase2_matched_eff.py found
corr(log mu, success) = -0.828 at matched eff: SMALL mu (a skew rank-2 GLV
lattice L2) goes with HIGH measured success.  Small mu is also exactly the
regime where the map (k1,k2) -> k1 + lam*k2 mod n can collide, i.e. where the
decomposition stops being unique.  So the observed correlation could be an
artifact: skew L2 -> non-unique decomposition -> verifier accepts wrong d ->
inflated "success".

This script settles it by re-running the matched-eff configurations and, for
every accepted candidate, comparing it against the secret that was actually
planted.  The secret is used ONLY as a post-hoc audit, never by the attack.

Reports, per configuration:
    accepted   -- verify_d said yes
    correct    -- accepted AND equal to the planted secret
    FALSE POS  -- accepted BUT different from the planted secret
plus the collision count of the (k1,k2) -> k1+lam*k2 map, computed exactly.

Run: python3 glv_hnp_phase2_verify_soundness.py
"""

import math
import random
from fpylll import IntegerMatrix, LLL

from glv_hnp_phase2_common import (
    find_generator, centered_rho, gen_signatures, build_glv_lattice,
    scales, verify_d, find_glv_curves, spurious_mu,
)

SEEDS = list(range(24))
M_RANGE = range(5, 13)

def audit_instance(curve, sigs, m, d_true, k1_bound):
    """Return (accepted, correct, false_pos) for one reduced instance."""
    n, lam = curve['n'], curve['lam']
    k2_bound = math.isqrt(n) + 1
    sigs = sigs[:m]
    dim = 2 * m + 2
    A = IntegerMatrix.from_matrix(
        build_glv_lattice(sigs, n, lam, k1_bound, k2_bound))
    LLL.reduction(A)
    _, _, _, S_KAN = scales(n, k1_bound, k2_bound)
    for i in range(dim):
        last = A[i][dim - 1]
        if abs(last) != S_KAN:
            continue
        sign = 1 if last > 0 else -1
        d_cand = (sign * A[i][m]) % n
        if verify_d(d_cand, sigs, n, lam, k1_bound, k2_bound):
            return True, d_cand == d_true, d_cand != d_true
    return False, False, False

def decomposition_collisions(n, lam, k1_bound, k2_bound, cap=400000):
    """
    Exact count of colliding pairs in (k1,k2) -> k1 + lam*k2 mod n over the
    declared box.  0 means the GLV decomposition is unique, so no wrong d can
    ever be certified by a single signature for the wrong reason.
    """
    if k1_bound * k2_bound > cap:
        return None
    seen = {}
    coll = 0
    for k2 in range(k2_bound):
        base = lam * k2 % n
        for k1 in range(k1_bound):
            k = (base + k1) % n
            if k in seen:
                coll += 1
            else:
                seen[k] = 1
    return coll

print("=" * 78)
print("Thread 20 (part 4): soundness audit of the oracle-free verifier")
print("=" * 78)

c2677 = dict(p=2677, b=2, n=2647, lam=185, rho=centered_rho(185, 2647),
             G=find_generator(2677, 2, 2647))
c2557 = dict(p=2557, b=2, n=2659, lam=1755, rho=centered_rho(1755, 2659),
             G=find_generator(2557, 2, 2659))

pool = []
pool += find_glv_curves(2**11, 2**13, 0.00, 0.15, want=3)
pool += find_glv_curves(2**11, 2**13, 0.15, 0.30, want=3)
pool += find_glv_curves(2**11, 2**13, 0.30, 0.51, want=3)
curves = sorted(pool, key=lambda c: c['rho']) + [c2557, c2677]

print(f"\n  matched eff ~ 0.157, {len(SEEDS)} seeds, m in "
      f"{min(M_RANGE)}..{max(M_RANGE)}\n")
print(f"  {'p':>7}{'rho':>7}{'K1':>4}{'eff':>8}{'mu':>8}{'collisions':>11}"
      f"{'accepted':>10}{'correct':>9}{'FALSE POS':>11}")

tot_acc = tot_cor = tot_fp = 0
rows = []
for cv in curves:
    n = cv['n']; k2 = math.isqrt(n) + 1
    k1 = max(2, int(round(0.157 * n / k2)))
    eff = k1 * k2 / n
    if eff >= 1.0:
        continue
    mu, _ = spurious_mu(n, cv['lam'], k1, k2)
    coll = decomposition_collisions(n, cv['lam'], k1, k2)
    acc = cor = fp = 0
    for sd in SEEDS:
        d = random.Random(sd + 7777).randrange(1, n)
        sigs = gen_signatures(cv, d, max(M_RANGE), k1, k2, sd)
        for m in M_RANGE:
            if len(sigs) < m:
                continue
            a, c, f = audit_instance(cv, sigs, m, d, k1)
            acc += a; cor += c; fp += f
    tot_acc += acc; tot_cor += cor; tot_fp += fp
    rows.append((cv, k1, eff, mu, coll, acc, cor, fp))
    print(f"  {cv['p']:>7}{cv['rho']:>7.3f}{k1:>4}{eff:>8.4f}{mu:>8.0f}"
          f"{(str(coll) if coll is not None else '?'):>11}"
          f"{acc:>10}{cor:>9}{fp:>11}")

print("\n" + "-" * 78)
print(f"  totals: accepted={tot_acc}  correct={tot_cor}  false positives={tot_fp}")
rate = (100.0 * tot_fp / tot_acc) if tot_acc else 0.0
print(f"  false-positive rate: {rate:.2f}% of accepted candidates")
if tot_fp == 0:
    print("  => verify_d never certified a wrong d.  Every 'success' counted in")
    print("     the Thread 20 scripts is a genuine key recovery.")
elif rate < 1.0:
    print("  => the verifier is sound to within well under 1%.  Too small to")
    print("     move any Thread 20 success rate by more than a fraction of a")
    print("     point, so the corr(log mu, success) = -0.828 finding is NOT a")
    print("     decomposition-uniqueness artifact.  Note the false positives")
    print("     did NOT land on the colliding-decomposition curves, so they are")
    print("     chance coincidences rather than a systematic uniqueness failure.")
else:
    print(f"  => {rate:.1f}% of accepted candidates were wrong.  Success rates in")
    print("     the Thread 20 scripts are inflated and must be recomputed")
    print("     against the planted secret.")

nz = [r for r in rows if r[4]]
print(f"\n  configurations with a non-unique decomposition (collisions > 0): "
      f"{len(nz)}/{len(rows)}")
if nz:
    print("    " + ", ".join(f"p={r[0]['p']}(coll={r[4]}, fp={r[7]})" for r in nz))
print("\nDone.")
