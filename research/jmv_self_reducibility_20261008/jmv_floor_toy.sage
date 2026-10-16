# Toy test of the JMV self-reducibility across endomorphism-ring levels.
# Surface curves E0 come from the Hilbert class polynomial of D_K, so
# End(E0) = O_K by construction.  For each prime ell exactly dividing the
# Frobenius conductor f_pi, E1 is one descending ell-step down (the floor).
#   T2  a DLP instance on E1 transfers to E0 through the ascending
#       ell-isogeny and the answer pulls back correctly;
#   T3  E0[ell] is rational exactly over F_{p^r}, r = ord(t/2 mod ell);
#       the vertical step is computed from F_{p^r} torsion (cost model);
#   T4  on the floor, the unique rational ell-subgroup of E1 is killed by
#       pi - lambda, i.e. the ascending kernel is E1[ell] ∩ ker(pi - lambda).
import json, sys, time
from sage.all import *

def frob_scalar_order(p, t, ell):
    lam = (t * inverse_mod(2, ell)) % ell
    return lam, Mod(lam, ell).multiplicative_order()

def order_over_ext(p, t, r):
    s = [2, t]
    for k in range(2, r + 1): s.append(t*s[-1] - p*s[-2])
    return p**r + 1 - s[r], s

def rational_ell_subgroups(E, p, t, ell, tries=600):
    """F_p-rational cyclic ell-subgroups of E (ell | t^2 - 4p), each as
    (kernel polynomial over F_p, generator over F_{p^r}, r)."""
    lam, r = frob_scalar_order(p, t, ell)
    Fr = GF(p**r); Er = E.change_ring(Fr)
    nr, _ = order_over_ext(p, t, r)
    k = 0
    while nr % ell**(k+1) == 0: k += 1
    cof = nr // ell**k
    pts = []
    for _ in range(tries):
        R = cof * Er.random_point()
        while not R.is_zero() and (ell*R).order() != 1 and R.order() != ell:
            R = ell*R
        if R.is_zero() or R.order() != ell: continue
        if not pts:
            pts.append(R)
        elif len(pts) == 1 and all(R != i*pts[0] for i in range(ell)):
            pts.append(R); break
    if not pts: raise RuntimeError("no ell-torsion found over F_{p^r}")
    if len(pts) == 1:
        gens = [pts[0]]
    else:
        P, Q = pts
        gens = [P + i*Q for i in range(ell)] + [Q]
    x = polygen(GF(p))
    out = []
    for G in gens:
        xs = [(i*G)[0] for i in range(1, (ell+1)//2)]
        hr = prod(polygen(Fr) - xi for xi in xs)
        try:
            h = hr.change_ring(GF(p))
        except (TypeError, ValueError):
            continue
        out.append((h, G, r))
    return out

def instances(ells, pmax_mult=40):
    out = []
    for ell in ells:
        found = 0
        for p in prime_range(max(50, ell*ell), pmax_mult*ell*ell + 2000):
            if found >= 2: break
            B = isqrt(4*p)
            for t in range(-B, B+1):
                D = t*t - 4*p
                if D >= 0: continue
                DK = fundamental_discriminant(D)
                f2 = D // DK
                if not is_square(f2): continue
                f = isqrt(f2)
                if f % ell != 0 or (f // ell) % ell == 0: continue
                H = hilbert_class_polynomial(DK)
                roots = H.change_ring(GF(p)).roots(multiplicities=False)
                if not roots: continue
                j = roots[0]
                if j == 0 or j == 1728: continue
                E = EllipticCurve_from_j(j)
                if E.trace_of_frobenius() != t: E = E.quadratic_twist()
                assert E.trace_of_frobenius() == t
                out.append((ell, p, E, t, D, DK, f)); found += 1
                break
    return out

def run(ell, p, E0, t, D, DK, f):
    n = E0.order()
    subs0 = rational_ell_subgroups(E0, p, t, ell)
    assert len(subs0) == ell + 1, (ell, len(subs0))        # surface: pi scalar on E0[ell]
    isos = [E0.isogeny(h) for h, _, _ in subs0]
    tv0 = time.time()
    down = [phi for phi in isos if len(rational_ell_subgroups(phi.codomain(), p, t, ell)) == 1]
    tv = time.time() - tv0
    horiz = ell + 1 - len(down)
    assert horiz == 1 + kronecker(DK, ell), (horiz, DK, ell)
    phi = down[0]; E1 = phi.codomain(); assert E1.order() == n
    # T2
    qs = [pf for pf, _ in factor(n) if pf != ell]
    if not qs: raise ValueError("no prime factor of n other than ell")
    q = max(qs); h = n // q
    P1 = h*E1.random_point()
    while P1.is_zero(): P1 = h*E1.random_point()
    x = randint(1, q-1); Q1 = x*P1
    subs1 = rational_ell_subgroups(E1, p, t, ell); assert len(subs1) == 1
    up = E1.isogeny(subs1[0][0]); to_E0 = up.codomain().isomorphism_to(E0)
    lift = lambda P: to_E0(up(P))
    x0 = discrete_log(lift(Q1), lift(P1), q, operation='+')
    t2 = (x0 == x)
    # T3
    lam, r = frob_scalar_order(p, t, ell)
    nr, s = order_over_ext(p, t, r)
    t3_a = (nr % (ell*ell) == 0)
    t3_b = all((p**k + 1 - s[k]) % ell != 0 for k in range(1, r))
    # T4
    Fr = GF(p**r); E1r = E1.change_ring(Fr)
    xs = subs1[0][0].change_ring(Fr).roots(multiplicities=False)
    t4 = all(E1r(T[0]**p, T[1]**p) == lam*T for T in (E1r.lift_x(xr) for xr in xs))
    return dict(ell=int(ell), p=int(p), t=int(t), D=int(D), DK=int(DK), f=int(f), n=int(n), order_q=int(q),
                horizontal=int(horiz), descending=len(down), T2_dlp_transfers=bool(t2),
                T3_r=int(r), T3_full_torsion_over_Fpr=bool(t3_a), T3_no_torsion_below_r=bool(t3_b),
                T3_descend_classify_seconds=float("%.3f" % tv),
                T4_ascending_kernel_is_pi_minus_lambda=bool(t4), T4_kernel_roots_checked=len(xs))

set_random_seed(20261008)
ells = [3, 5, 7, 11, 13, 17, 19]
res = []
for inst in instances(ells):
    try:
        r = run(*inst); res.append(r); print(json.dumps(r)); sys.stdout.flush()
    except Exception as e:
        print('ERROR', int(inst[0]), int(inst[1]), repr(e)[:160]); sys.stdout.flush()
import os; os.makedirs('results', exist_ok=True)
with open('results/pilot_results.jsonl', 'w') as fh:
    for r in res: fh.write(json.dumps(r) + '\n')
print('summary', json.dumps(dict(instances=len(res),
      T2=all(r['T2_dlp_transfers'] for r in res),
      T3=all(r['T3_full_torsion_over_Fpr'] and r['T3_no_torsion_below_r'] for r in res),
      T4=all(r['T4_ascending_kernel_is_pi_minus_lambda'] for r in res))))
