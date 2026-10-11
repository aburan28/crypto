"""
Eigen-Seta toy prototype (correctness only; toy parameters, no security).

Scheme (RESEARCH_ENDRING_PKE_STRATEGY.md, section 4, with the parity fix
found while implementing: decrypt with the trace-zero element
theta' = 2*theta - tr(theta), which has odd norm Delta' = 4w - t^2).

  pk = (E_A, P, Q, Pd, Qd):   (P, Q) eigenbasis of E_A[2^e] for a secret
       endomorphism theta of E_A (eigenvalues mu, nu),  (Pd, Qd) basis of E_A[3^f].
  ct = (E_C, S, T) = (E_A/<Pd + s Qd>, [a] phi(P), [b] phi(Q)),  a, b random odd.
  dec: psi' = phi theta' phi_dual acts on S, T by +d*s0, -d*s0 (s0 = mu - nu);
       the masks cancel.  alpha = psi' + c,  F = [[alpha, c'], [-c', alpha_dual]]
       with deg(alpha) + c'^2 = d^2 Delta' + c^2 + c'^2 = 2^e is a (2^e,2^e)-isogeny
       of E_C x E_C with kernel <(c'S, -(d s0 + c) S), (c'T, -(c - d s0) T)>.
       Evaluating F on E_C[3^f] gives alpha, hence ker(phi_dual) = ker(alpha - c),
       hence s.

Key variants implemented here:
  cm    : E_A = curve with CM by Z[(1+sqrt(-7))/2], theta = omega of degree 2
          (trace +-1, norm 2, Delta' = 7).  Insecure as a key (End public) but
          exercises the complete mechanism with a tiny Kani chain.
  walk  : E_A = eps(E_cm) for a secret 5^k-isogeny eps, theta = eps omega eps_dual
          (trace +-5^k, norm 2*25^k, Delta' = 7*25^k).  The "hidden special
          curve" variant (Castryck-Vercauteren section 5.3 backdoor, used by
          its owner).  Chain condition 2^e > 7 * 25^k * 9^f.
Both use a single (2^e,2^e) Richelot chain, i.e. the N > deg regime.  The
N^2 > deg two-halves evaluation (needed for random curves) is not implemented.
"""
import sys, time, random
from sage.all import (ZZ, GF, EllipticCurve, EllipticCurve_from_j, is_prime, two_squares,
                      discrete_log, inverse_mod, randint, set_random_seed)
from kani_chain import chain_evaluate

# ----------------------------------------------------------------------------
def find_prime(e, f, k):
    """p = 2^(e+1) 3^f 5^k c - 1 prime with p = 3,5,6 mod 7 (so that the CM
    curve with discriminant -7 is supersingular)."""
    base = ZZ(2) ** (e + 1) * ZZ(3) ** f * ZZ(5) ** k
    c = 1
    while True:
        p = base * c - 1
        if is_prime(p) and p % 7 in (3, 5, 6):
            return p, c
        c += 1

def random_point_of_order(E, n, cofactor):
    primes = n.prime_factors()
    while True:
        R = cofactor * E.random_point()
        if (n * R).is_zero() and all(not ((n // l) * R).is_zero() for l in primes):
            return R

def torsion_basis(E, n, cofactor):
    l = n.prime_factors()[0]
    P = random_point_of_order(E, n, cofactor)
    while True:
        Q = random_point_of_order(E, n, cofactor)
        if P.weil_pairing(Q, n) ** (n // l) != 1:
            return P, Q

def bidlp(R, P, Q, n):
    ePQ = P.weil_pairing(Q, n)
    x = ZZ(discrete_log(R.weil_pairing(Q, n), ePQ, ord=n))
    y = ZZ(discrete_log(P.weil_pairing(R, n), ePQ, ord=n))
    assert R == x * P + y * Q
    return x, y

def sqrt_mod_2power(a, k):
    a = ZZ(a) % 2 ** k
    assert a % 8 == 1
    s = ZZ(1)
    for j in range(3, k):
        if (s * s - a) % 2 ** (j + 1) != 0:
            s += 2 ** (j - 1)
    assert (s * s - a) % 2 ** k == 0
    return s

def full_order_curve(F, j, p):
    """F_{p^2}-model of the supersingular curve with invariant j having
    (p+1)^2 points (Frobenius = -p)."""
    E = EllipticCurve_from_j(F(j))
    R = E.random_point()
    if not ((p + 1) * R).is_zero():
        E = E.quadratic_twist()
        R = E.random_point()
        assert ((p + 1) * R).is_zero()
    E.set_order((p + 1) ** 2)
    return E

# ----------------------------------------------------------------------------
def cm7_endomorphism(E, p):
    """omega with omega^2 - omega + 2 = 0 (or its negative) on the CM(-7)
    curve E, as (evaluator, trace, norm).  omega = iso o (2-isogeny)."""
    for T in E(0).division_points(2):
        if T.is_zero():
            continue
        phi2 = E.isogeny(T)
        E2 = phi2.codomain()
        if not E2.is_isomorphic(E):
            continue
        for iso in [u * E2.isomorphism_to(E) for u in E.automorphisms()]:
            def th(R, _iso=iso, _phi2=phi2):
                return _iso(_phi2(R))
            R = E.random_point()
            if (th(th(R)) - th(R) + 2 * R).is_zero():
                return th, ZZ(1), ZZ(2)
            if (th(th(R)) + th(R) + 2 * R).is_zero():
                return th, ZZ(-1), ZZ(2)
    raise RuntimeError("no degree-2 endomorphism found (is p inert in Q(sqrt(-7))?)")

def admissible_kani_scalars(e, d, Dp):
    """largest e' <= e with 2^e' - d^2 Delta' = c^2 + c'^2, c even, c' odd"""
    for ep in range(e, 0, -1):
        Mp = ZZ(2) ** ep - d * d * Dp
        if Mp <= 0:
            continue
        try:
            cc, cp = two_squares(Mp)
        except ValueError:
            continue
        if cc % 2 == 1:
            cc, cp = cp, cc
        if cp % 2 == 1:
            return ep, ZZ(cc), ZZ(cp)
    raise RuntimeError("no admissible (c, c') up to 2^e; increase e")

def keygen(p, e, f, k, variant="walk", verbose=True):
    F = GF(p ** 2, name="z2", modulus=[1, 0, 1])
    d = ZZ(3) ** f
    cof = lambda n: (p + 1) // n
    Ecm = full_order_curve(F, -3375, p)
    omega, t_cm, w_cm = cm7_endomorphism(Ecm, p)
    if variant == "cm" or k == 0:
        EA, theta_eval, t, w = Ecm, omega, t_cm, w_cm
    else:
        DA = ZZ(5) ** k
        PD, QD = torsion_basis(Ecm, DA, cof(DA))
        eps = Ecm.isogeny(PD + randint(0, DA - 1) * QD, algorithm="factored")
        EA, eps_dual = eps.codomain(), None
        eps_dual = eps.dual()
        def theta_eval(R):          # theta = eps o omega o eps_dual
            return eps(omega(eps_dual(R)))
        t, w = DA * t_cm, DA * DA * w_cm
    Dp = 4 * w - t * t
    assert Dp % 8 == 7
    ep, cc, cp = admissible_kani_scalars(e, d, Dp)
    N2 = ZZ(2) ** ep
    assert cc * cc + cp * cp + d * d * Dp == N2 and cc % 2 == 0 and cp % 2 == 1

    # eigenbasis of theta on E_A[2^e']
    P1, Q1 = torsion_basis(EA, N2, cof(N2))
    a11, a21 = bidlp(theta_eval(P1), P1, Q1, N2)
    a12, a22 = bidlp(theta_eval(Q1), P1, Q1, N2)
    assert (a11 + a22 - t) % N2 == 0 and (a11 * a22 - a12 * a21 - w) % N2 == 0
    s0 = sqrt_mod_2power(-Dp, ep + 2)
    mu = ZZ(((t + s0) // 2) % N2); nu = ZZ(((t - s0) // 2) % N2)
    def eigvec(lam):
        for v in ((-a12, a11 - lam), (a22 - lam, -a21)):
            if v[0] % 2 == 1 or v[1] % 2 == 1:
                return v
        raise RuntimeError("no eigenvector of full order")
    v = eigvec(mu); P = v[0] * P1 + v[1] * Q1
    v = eigvec(nu); Q = v[0] * P1 + v[1] * Q1
    assert theta_eval(P) == mu * P and theta_eval(Q) == nu * Q
    assert P.weil_pairing(Q, N2) ** (N2 // 2) != 1
    Pd, Qd = torsion_basis(EA, d, cof(d))

    pk = dict(EA=EA, P=P, Q=Q, Pd=Pd, Qd=Qd)
    sk = dict(t=t, w=w, Dp=Dp, s0=s0 % N2, c=cc, cp=cp, mu=mu, nu=nu, theta_eval=theta_eval)
    par = dict(F=F, N2=N2, e=ep, d=d, f=f, p=p)
    if verbose:
        print(f"  keygen[{variant}]: theta trace {t}, norm {w}, Delta' {Dp}; using 2^{ep} torsion, "
              f"c={cc}, c'={cp}; d^2 Delta'/2^e' = {float(d*d*Dp)/float(N2):.3f}")
    return pk, sk, par

# ----------------------------------------------------------------------------
def encaps(pk, par, s=None, a=None, b=None):
    N2, d = par["N2"], par["d"]
    s = randint(0, d - 1) if s is None else s
    a = 2 * randint(0, N2 // 2 - 1) + 1 if a is None else a
    b = 2 * randint(0, N2 // 2 - 1) + 1 if b is None else b
    phi = pk["EA"].isogeny(pk["Pd"] + s * pk["Qd"], algorithm="factored")
    ct = dict(EC=phi.codomain(), S=a * phi(pk["P"]), T=b * phi(pk["Q"]))
    return ct, dict(s=s, a=a, b=b, phi=phi)

def lemma_check(pk, sk, par, ct, wit):
    """independent check of the decoding lemma: phi theta' phi_dual (S) = d s0 S"""
    phi, d, N2 = wit["phi"], par["d"], par["N2"]
    th = sk["theta_eval"]; t = sk["t"]
    def theta_prime(R):
        return 2 * th(R) - t * R
    phid = phi.dual()
    for X, sign in ((ct["S"], 1), (ct["T"], -1)):
        lhs = phi(theta_prime(phid(X)))
        assert lhs == (sign * d * sk["s0"]) * X, "lollipop eigenvalue mismatch"
    return True

def decaps(pk, sk, par, ct, verbose=True):
    N2, d, e = par["N2"], par["d"], par["e"]
    EC, S, T = ct["EC"], ct["S"], ct["T"]
    c, cp, s0 = sk["c"], sk["cp"], sk["s0"]
    lamS = (d * s0 + c) % N2; lamT = (c - d * s0) % N2        # alpha on S, T
    assert lamS % 2 == 1 and lamT % 2 == 1 and cp % 2 == 1
    cofd = (par["p"] + 1) // d
    # The Mumford-coordinate formulas of the splitting step have rare degenerate
    # inputs (divisors through a Weierstrass point etc.).  They depend on the
    # payload points, so on failure we retry with a fresh basis of E_C[3^f].
    last_err = None
    for attempt in range(6):
        try:
            t0 = time.time()
            X1, X2 = torsion_basis(EC, d, cofd)
            O = EC(0)
            K1 = (cp * S, (-lamS) * S)
            K2 = (cp * T, (-lamT) * T)
            payload = [(X1, O), (O, X1), (X2, O)]
            (E1, E2), imgs, log = chain_evaluate(EC, EC, K1, K2, e, payload)
            if verbose:
                from collections import Counter
                print(f"  decaps: (2^{e},2^{e})-chain in {time.time()-t0:.2f}s; steps: {dict(Counter(log))}; "
                      f"first steps {log[:4]}")
            def iso_candidates(E):
                base = E.isomorphism_to(EC)
                return [u * base for u in EC.automorphisms()]
            # F(X,0) = (alpha X, -c' X),  F(0,Y) = (c' Y, alpha_dual Y) fix the component
            # isomorphisms; the chain may return the two factors swapped.
            (A1, A2), (Z1, Z2), (B1, B2) = imgs
            found = None
            for swap in (False, True):
                a1, a2, z1, z2, b1, b2 = (A1, A2, Z1, Z2, B1, B2) if not swap else (A2, A1, Z2, Z1, B2, B1)
                Ea, Eb = (E1, E2) if not swap else (E2, E1)
                for iso2 in iso_candidates(Eb):
                    if iso2(a2) != -cp * X1 or iso2(b2) != -cp * X2:
                        continue
                    for iso1 in iso_candidates(Ea):
                        if iso1(z1) == cp * X1:
                            found = (iso1(a1), iso1(b1))
                            break
                    if found:
                        break
                if found:
                    break
            if found is None:
                raise RuntimeError("could not match component isomorphism")
            alphaX1, alphaX2 = found
            break
        except (TypeError, ZeroDivisionError, AssertionError, RuntimeError) as ex:
            last_err = ex
            if verbose:
                print(f"  decaps: degenerate chain evaluation, retrying ({type(ex).__name__}: {str(ex)[:50]})")
    else:
        raise RuntimeError(f"chain evaluation failed repeatedly: {last_err}")

    psiX1 = alphaX1 - c * X1
    psiX2 = alphaX2 - c * X2
    m11, m21 = bidlp(psiX1, X1, X2, d)
    m12, m22 = bidlp(psiX2, X1, X2, d)
    # ker(phi_dual) is a cyclic subgroup of order d inside ker(psi') cap E_C[d].
    # Generically ker(psi') cap E_C[d] is exactly ker(phi_dual); when theta'
    # stabilises part of ker(phi) the kernel is larger (CV section 3.2) and we
    # try each cyclic order-d subgroup (toy: brute force over (Z/d)^2).
    kernel_vecs = [(x, y) for x in range(d) for y in range(d)
                   if (m11 * x + m12 * y) % d == 0 and (m21 * x + m22 * y) % d == 0]
    assert len(kernel_vecs) >= d
    seen, candidates = set(), []
    for (x, y) in kernel_vecs:
        if x % 3 == 0 and y % 3 == 0:
            continue                                    # order < d
        u = x if x % 3 else y                            # normalise the generator
        inv = inverse_mod(u, d)
        key = ((x * inv) % d, (y * inv) % d)
        if key not in seen:
            seen.add(key); candidates.append(key)
    if verbose and len(candidates) > 1:
        print(f"  decaps: ker(psi') on E_C[3^f] has order {len(kernel_vecs)}; "
              f"{len(candidates)} cyclic candidates for ker(phi_dual)")
    phi_dual = None
    for (x, y) in candidates:
        Khat = x * X1 + y * X2
        cand = EC.isogeny(Khat, algorithm="factored")
        if cand.codomain().is_isomorphic(pk["EA"]):
            phi_dual = cand
            break
    if phi_dual is None:
        raise RuntimeError("no kernel candidate lands on E_A")
    EAp = phi_dual.codomain()
    iso = EAp.isomorphism_to(pk["EA"])
    Krec = None
    for Xk in (X1, X2, X1 + X2):                 # need phi_dual(Xk) of order exactly d
        img = phi_dual(Xk)
        if not ((d // 3) * img).is_zero():
            Krec = iso(img)
            break
    assert Krec is not None
    x, y = bidlp(Krec, pk["Pd"], pk["Qd"], d)
    assert x % 3 != 0
    return ZZ((y * inverse_mod(x, d)) % d)

# ----------------------------------------------------------------------------
def main(e=24, f=4, k=2, trials=3, variant="walk", seed=1):
    set_random_seed(seed); random.seed(seed)
    p, c = find_prime(e, f, k)
    print(f"variant={variant}: p = 2^{e+1} 3^{f} 5^{k} * {c} - 1  ({p.nbits()} bits)")
    t0 = time.time()
    pk, sk, par = keygen(p, e, f, k, variant=variant)
    print(f"  keygen time {time.time()-t0:.2f}s")
    ok = 0
    for i in range(trials):
        ct, wit = encaps(pk, par)
        lemma_check(pk, sk, par, ct, wit)
        t1 = time.time()
        s_rec = decaps(pk, sk, par, ct)
        print(f"  trial {i}: s = {wit['s']}, recovered {s_rec}, masks a={wit['a']} b={wit['b']} "
              f"({time.time()-t1:.2f}s)")
        ok += int(s_rec == wit["s"])
    print(f"RESULT[{variant}]: {ok}/{trials} decapsulations correct")
    return ok == trials

if __name__ == "__main__":
    e = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    f = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    trials = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    variant = sys.argv[5] if len(sys.argv) > 5 else "walk"
    sys.exit(0 if main(e, f, k, trials, variant) else 1)
