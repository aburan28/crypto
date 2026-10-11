"""
Generic (2^e, 2^e)-isogeny chain between products of elliptic curves whose
intermediate steps may be of any of the four kinds:

  product  -> product   (kernel = graph of an isomorphism on E1[2], or a
                         product kernel <(a,0),(0,b)>)
  product  -> Jacobian  (gluing, Castryck-Decru / Oudompheng formula)
  Jacobian -> Jacobian  (Richelot correspondence)
  Jacobian -> product   (splitting, Smith's formula)

The existing FESTA / Castryck-Decru code assumes glue-once / split-once.  The
Eigen-Seta decryption needs the general case because the embedded endomorphism
is diagonal on the eigenbasis, which makes the first kernel 2-torsion
degenerate (the diagonal of E_C[2] x E_C[2]).

Jacobian-side formulas are taken from richelot_aux.py (MIT, Pope-Oudompheng-
Panny); product-side steps use Sage's elliptic isogenies/isomorphisms.
Points are Sage points (the zero point for "nothing"), divisors are Mumford
pairs (u, v) of polynomials over F_{p^2}.
"""
from sage.all import Matrix, HyperellipticCurve, PolynomialRing
from richelot_aux import (FromProdToJac, FromJacToProd, RichelotCorr,
                          jacobian_iter_double)


def _order2_pair(K, m):
    return (m * K[0], m * K[1])


def _is_zero_pt(P):
    return P.is_zero()


def _product_step(E1, E2, A, B):
    """A, B order-2 pairs generating a maximal isotropic subgroup of
    (E1 x E2)[2].  Returns ('prod', map, E1', E2') or None if the kernel is a
    gluing kernel."""
    elems = [A, B, (A[0] + B[0], A[1] + B[1])]
    # product kernel <(a,0),(0,b)> ?
    a = [z[0] for z in elems if _is_zero_pt(z[1]) and not _is_zero_pt(z[0])]
    b = [z[1] for z in elems if _is_zero_pt(z[0]) and not _is_zero_pt(z[1])]
    if a and b:
        fa = E1.isogeny(a[0])
        fb = E2.isogeny(b[0])
        def m(pair):
            return (fa(pair[0]), fb(pair[1]))
        return ("prod", m, fa.codomain(), fb.codomain())
    # graph of an isomorphism sigma : E1 -> E2 on E1[2] ?
    if any(_is_zero_pt(z[0]) or _is_zero_pt(z[1]) for z in elems):
        return None
    if E1.is_isomorphic(E2):
        for sigma in E1.isomorphisms(E2):
            if sigma(A[0]) == A[1] and sigma(B[0]) == B[1]:
                sinv = sigma.dual() if hasattr(sigma, "dual") else ~sigma
                def m(pair, sinv=sinv):
                    y = sinv(pair[1])
                    return (pair[0] + y, pair[0] - y)
                return ("prod", m, E1, E1)
    return None


def chain_evaluate(E1, E2, K1, K2, e, payload, verbose=False):
    """Compute the (2^e,2^e)-isogeny of E1 x E2 with kernel <K1, K2> (pairs of
    points of order 2^e, maximal isotropic) and push `payload` (list of pairs
    of points) through it.  Returns (codomain (E1', E2'), images, log)."""
    state = ("prod", E1, E2)
    kers = [K1, K2]
    pay = list(payload)
    log = []
    for j in range(e, 0, -1):
        m = 2 ** (j - 1)
        if state[0] == "prod":
            _, C1, C2 = state
            A = _order2_pair(kers[0], m)
            B = _order2_pair(kers[1], m)
            ps = _product_step(C1, C2, A, B)
            if ps is not None:
                _, f, D1, D2 = ps
                kers = [f(K) for K in kers]
                pay = [f(P) for P in pay]
                state = ("prod", D1, D2)
                log.append("prod")
                continue
            # gluing into a Jacobian
            h, u1, v1, u2, v2, f = FromProdToJac(C1, C2, kers[0][0], kers[1][0],
                                                 kers[0][1], kers[1][1], j)
            kers = [(u1, v1), (u2, v2)]
            def fpair(P, f=f):
                D = f((P[0] if not P[0].is_zero() else None,
                       P[1] if not P[1].is_zero() else None))
                return (D[0], D[1])
            pay = [fpair(P) for P in pay]
            state = ("jac", h)
            log.append("glue")
            continue
        # Jacobian: find the order-2 part of the kernel
        _, h = state
        G1, _ = jacobian_iter_double(h, kers[0][0], kers[0][1], j - 1)
        G2, _ = jacobian_iter_double(h, kers[1][0], kers[1][1], j - 1)
        G3, r3 = h.quo_rem(G1 * G2)
        assert r3 == 0
        delta = Matrix(G.padded_list(3) for G in (G1, G2, G3))
        if delta.determinant() == 0:
            f, (D1, D2) = FromJacToProd(G1, G2, G3)
            kers = [f(K) for K in kers] if j > 1 else []
            pay = [f(P) for P in pay]
            state = ("prod", D1, D2)
            log.append("split")
            continue
        x = h.parent().gen()
        dinv = delta.inverse()
        H1 = -dinv[0][0] * x ** 2 + 2 * dinv[1][0] * x - dinv[2][0]
        H2 = -dinv[0][1] * x ** 2 + 2 * dinv[1][1] * x - dinv[2][1]
        H3 = -dinv[0][2] * x ** 2 + 2 * dinv[1][2] * x - dinv[2][2]
        hnew = H1 * H2 * H3
        R = RichelotCorr(G1, G2, H1, H2, hnew)
        kers = [R.map(K) for K in kers] if j > 1 else []
        pay = [R.map(P) for P in pay]
        state = ("jac", hnew)
        log.append("rich")
    if state[0] != "prod":
        raise RuntimeError("chain did not end on a product: " + ",".join(log))
    return (state[1], state[2]), pay, log
