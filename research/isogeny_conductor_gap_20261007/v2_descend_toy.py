#!/usr/bin/env sage -python
"""V2 toy: descend one volcano level from the Koblitz curve K_0 at a toy rung.

Run through the checked Sage launcher:
    /Volumes/SSD990/cryptanalysis/sage -python v2_descend_toy.py --n 23 --ell 967 --seed 1

What it does (all exact, no DLP):
  1. K_0: y^2 + xy = x^3 + 1 over F_q, q = 2^n.  #K_0 = q + 1 - t from the Lucas
     recurrence (checked against Sage's point count), conductor f with
     t^2 - 4q = -7 f^2, and the requested prime ell | f.
  2. d = ord(t/2 mod ell): the Frobenius pi acts on E[ell] as the scalar t/2
     (K_0 is on the surface), so E[ell] is rational over F_{q^d}.
  3. Find a point P of order ell in E(F_{q^d}); reject the two tau-eigenlines
     (horizontal edges); the remaining ell - 1 kernels are descending.
  4. Vélu isogeny phi: E -> E' of degree ell over F_{q^d}; j(E') must lie in F_q
     (the kernel subgroup is F_q-rational); descend: build E_1 over F_q with
     j(E_1) = j(E') and the twist whose order equals #K_0.
  5. Transfer the generator: phi(G) has the same order as G; record it.
Outputs a JSON record on stdout; exits non-zero if any check fails.
"""
import argparse, json, sys, time
from sage.all import GF, EllipticCurve, PolynomialRing, ZZ, Integer, randint, set_random_seed


def lucas(n):
    V = [2, -1]; U = [0, 1]
    for k in range(1, n):
        V.append(-V[k] - 2 * V[k - 1]); U.append(-U[k] - 2 * U[k - 1])
    return V[n], abs(U[n])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=23)
    ap.add_argument("--ell", type=int, default=967)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-tries", type=int, default=50)
    a = ap.parse_args()
    set_random_seed(a.seed)
    t0 = time.time()
    n, ell = a.n, a.ell
    q = ZZ(2) ** n
    t, f = lucas(n)
    t, f = ZZ(t), ZZ(f)
    assert t * t - 4 * q == -7 * f * f, "conductor identity failed"
    assert f % ell == 0, "ell does not divide the conductor"
    N = q + 1 - t
    chi = Integer(-7).kronecker(ell)
    a_scalar = (t * Integer(2).inverse_mod(ell)) % ell
    # d = multiplicative order of a mod ell, computed by hand
    x, dd = a_scalar, 1
    while x != 1:
        x = (x * a_scalar) % ell; dd += 1
        if dd > ell: raise SystemExit("order computation failed")
    d = dd
    F = GF(q, "a")
    E = EllipticCurve(F, [1, 0, 0, 0, 1])          # y^2 + xy = x^3 + 1
    N_sage = E.cardinality()
    assert N_sage == N, f"point count mismatch: sage {N_sage} vs lucas {N}"
    # extension field containing E[ell]
    K = GF(q ** d, "b")
    EK = EllipticCurve(K, [1, 0, 0, 0, 1])
    # #E(F_{q^d}) from the Frobenius eigenvalues: alpha^d + beta^d via Lucas on (alpha,beta) with alpha+beta=t, alpha*beta=q
    Vd = [2, t]
    for k in range(1, d):
        Vd.append(t * Vd[k] - q * Vd[k - 1])
    NK = q ** d + 1 - Vd[d]
    v = 0; m = NK
    while m % ell == 0:
        m //= ell; v += 1
    assert v >= 2, f"E[ell] not fully rational over F_q^d: v_ell(#E(K)) = {v}"
    cof = NK // (ell ** v)
    # find a point of exact order ell whose line is not a tau-eigenline
    def frob2(P):  # 2-power Frobenius = tau on points
        return EK(P[0] ** 2, P[1] ** 2)
    tries = 0; P = None
    while tries < a.max_tries:
        tries += 1
        R = EK.random_point()
        S = cof * R
        # reduce to exact order ell
        while S != EK(0) and (ell * S) != EK(0):
            S = ell * S
        if S == EK(0):
            continue
        # eigenline test: tau(S) in <S>  <=> tau(S) and S linearly dependent;
        # check via Weil pairing e(S, tau(S)) == 1
        if S.weil_pairing(frob2(S), ell) == 1:
            continue  # horizontal edge; skip
        P = S; break
    if P is None:
        raise SystemExit("no descending kernel found")
    phi = EK.isogeny(P)             # Vélu, degree ell, over K
    Eprime = phi.codomain()
    jp = Eprime.j_invariant()
    assert jp ** q == jp, "j(E') not in F_q: kernel was not F_q-rational"
    assert jp ** 2 != jp, "j(E') in F_2: that is a Koblitz curve, not a descent"
    # descend: minimal polynomial of j' over F_2 has a root in F_q
    mp = jp.minimal_polynomial()
    assert mp.degree() <= n
    roots = mp.roots(F, multiplicities=False)
    assert roots, "no root of minpoly(j') in F_q"
    j1 = roots[0]
    # E_1: y^2 + xy = x^3 + a2 x^2 + a6 with j = 1/a6 ; choose twist by trace of a2
    a6 = 1 / j1
    cands = []
    for a2 in [F(0), F(1)]:
        E1 = EllipticCurve(F, [1, a2, 0, 0, a6])
        cands.append((E1.cardinality(), a2))
    match = [a2 for (c, a2) in cands if c == N]
    assert match, f"neither twist has order {N}: {[c for c, _ in cands]}"
    E1 = EllipticCurve(F, [1, match[0], 0, 0, a6])
    # transfer a generator: order of phi(G) equals order of G
    G = E.random_point()
    GK = EK(G)
    G1 = phi(GK)
    oG = G.order()
    # order check without factoring #E(K): oG*G1 = O and (oG/p)*G1 != O for every prime p | oG
    assert oG * G1 == EK(0), "phi(G) not killed by ord(G)"
    for (p_, _) in ZZ(oG).factor():
        assert (oG // p_) * G1 != EK(0), "phi(G) has smaller order than G"
    oG1 = oG
    # Galois orbit size of j1 under the 2-Frobenius
    orb = 1; y = j1 ** 2
    while y != j1:
        y = y ** 2; orb += 1
    rec = {
        "n": n, "ell": ell, "q_bits": n, "trace": int(t), "conductor_f": int(f),
        "chi_minus7_ell": int(chi), "pi_mod_ell": int(a_scalar), "kernel_field_degree_d": d,
        "order_K0": int(N), "v_ell_of_order_over_extension": v,
        "kernel_search_tries": tries, "isogeny_degree": int(phi.degree()),
        "j_E1_in_Fq_not_F2": True, "galois_orbit_size_of_j1": orb,
        "E1_a2_twist": str(match[0]), "E1_order_equals_K0_order": True,
        "generator_order": int(oG), "phi_G_order": int(oG1),
        "wall_s": round(time.time() - t0, 2), "status": "PASS",
    }
    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
