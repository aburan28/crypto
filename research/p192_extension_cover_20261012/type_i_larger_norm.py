"""Point-level Type I quotient and three-step norm checks for P-224/P-256.

Uses the checked Sage launcher. Search is deterministic in m=1..32 and
selects the first m for which the two cover fibers above mG and 2mG split
over F_(p^6). Exact genus-three class conclusions use degree-two push-pull
and the Klein-four quotient identity as recorded in BRANCH_CONTINUATION.md.
"""

import hashlib
import itertools
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, PolynomialRing, matrix, prod


HERE = Path(__file__).resolve().parent
CURVES = {
    "P-224": {
        "p": "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001",
        "b": "B4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4",
        "gx": "B70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21",
        "gy": "BD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34",
        "r": "FFFFFFFFFFFFFFFFFFFFFFFFFFFF16A2E0B8F03E13DD29455C5C2A3D",
    },
    "P-256": {
        "p": "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF",
        "b": "5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B",
        "gx": "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296",
        "gy": "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5",
        "r": "FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551",
    },
}


def run(name, params):
    p, b, gx, gy, r = [Integer(params[k], 16) for k in ("p", "b", "gx", "gy", "r")]
    a = p - 3
    q = p**2
    F = GF(p)
    R = PolynomialRing(F, "x")
    x = R.gen()
    f = x**3 + F(a) * x + F(b)
    assert f.is_irreducible()
    # Native absolute finite-field arithmetic is essential here: a quadratic
    # polynomial quotient over F_(p^3) makes the large Frobenius powers slow.
    L = GF(p**6, "u")
    assert L.cardinality() == p**6
    two_division_roots = f.change_ring(L).roots(multiplicities=False)
    assert len(two_division_roots) == 3
    z = two_division_roots[0]
    eps = z
    eps_q = eps**q
    lam = (z**(p**2) - z) / (z**p - z)
    mu = (lam * eps - eps_q) / (lam - 1)
    nu = (lam * eps_q - eps) / (lam - 1)
    A = matrix(L, [[nu, -eps * eps_q], [1, -mu]])
    frob = lambda M, power: matrix(L, 2, 2, [v**power for v in M.list()])
    B = frob(A, q**2) * frob(A, q) * A
    assert B[1, 0] != 0
    disc = (B[1, 1] - B[0, 0])**2 + 4 * B[1, 0] * B[0, 1]
    assert disc.is_square()
    root = disc.sqrt()
    betas = []
    for sgn in (1, -1):
        beta = (B[0, 0] - B[1, 1] + sgn * root) / (2 * B[1, 0])
        den = A[1, 0] * beta + A[1, 1]
        if den == 0 or beta**q != (A[0, 0] * beta + A[0, 1]) / den:
            continue
        if len({eps, eps_q, beta, beta**q}) != 4:
            continue
        replay = ((beta - eps_q) / (beta - eps)) * ((beta**q - eps) / (beta**q - eps_q))
        if replay == lam:
            betas.append(beta)
    assert betas
    beta = betas[0]
    quartic_roots = [eps, eps_q, beta, beta**q]
    elliptic_roots = [z, z**p, z**(p**2)]
    P = PolynomialRing(L, "X")
    X = P.gen()
    Q = prod(X - v for v in quartic_roots)
    assert Q.is_squarefree()

    models = []
    for rho in quartic_roots:
        finite = [v for v in quartic_roots if v != rho]
        for targets in itertools.permutations(elliptic_roots):
            reciprocal = lambda v: 1 / (v - rho)
            D = (targets[0] - targets[1]) / (reciprocal(finite[0]) - reciprocal(finite[1]))
            C = targets[0] - D * reciprocal(finite[0])
            if C + D * reciprocal(finite[2]) != targets[2]:
                continue
            linear = C * (X - rho) + D
            cleared = linear**3 * (X - rho) + L(a) * linear * (X - rho)**3 + L(b) * (X - rho)**4
            scale = cleared.leading_coefficient()
            assert cleared == scale * Q
            if scale.is_square():
                models.append((rho, C, D, scale.sqrt()))
    assert models
    rho, C, D, kappa = models[0]
    E = EllipticCurve(L, [L(a), L(b)])
    G = E(L(gx), L(gy))
    O = E(0)
    assert G != O and r.is_prime() and r * G == O

    def quotient(Pt, i):
        x0, ys = Pt
        yi = ys[i]
        rhoi, Ci, Di, ki = [v**(q**i) for v in (rho, C, D, kappa)]
        assert x0 != rhoi
        return E(Ci + Di / (x0 - rhoi), ki * yi / (x0 - rhoi)**2)

    def fiber(point):
        if point == O or point[0] == C:
            return None
        x0 = rho + D / (point[0] - C)
        y0 = point[1] * (x0 - rho)**2 / kappa
        rhs = [prod(x0 - v**(q**i) for v in quartic_roots) for i in range(3)]
        if y0 == 0 or y0**2 != rhs[0] or rhs[1] == 0 or not rhs[1].is_square():
            return None
        y1 = rhs[1].sqrt()
        norm_u = prod((x0 - eps**(q**i)) * (x0 - beta**(q**i)) for i in range(3))
        y2 = norm_u / (y0 * y1)
        assert y2**2 == rhs[2]
        plus = (x0, (y0, y1, y2))
        minus = (x0, (y0, -y1, -y2))
        assert quotient(plus, 0) == point and quotient(minus, 0) == point
        assert quotient(plus, 1) == -quotient(minus, 1)
        assert quotient(plus, 2) == -quotient(minus, 2)
        return [plus, minus]

    def sigma(Pt):
        x0, (y0, y1, y2) = Pt
        # These quotient labels are q^i, unlike the p^i labels in the
        # P-192 coordinate verifier. The q-action sends f_i to f_(i+1).
        result = (x0**q, (y2**q, y0**q, y1**q))
        assert prod(result[1]) == prod(Pt[1])**q
        return result

    attempted = []
    chosen = None
    for m in range(1, 33):
        positive = fiber(m * G)
        negative = fiber(2 * m * G)
        attempted.append({"m": m, "fiber_mG_split": positive is not None, "fiber_2mG_split": negative is not None})
        if positive is not None and negative is not None:
            chosen = (m, positive, negative)
            break
    assert chosen, f"no paired split fibers through m=32 for {name}"
    m, positive, negative = chosen
    target = -2 * m * G
    assert target != O and r * target == O
    vectors = []
    for j in range(3):
        pos = positive
        neg = negative
        for _ in range(j):
            pos = [sigma(Pt) for Pt in pos]
            neg = [sigma(Pt) for Pt in neg]
        row = [sum((quotient(Pt, i) for Pt in pos), O) -
               sum((quotient(Pt, i) for Pt in neg), O) for i in range(3)]
        expected = [O, O, O]
        expected[j % 3] = target
        assert row == expected, (name, j, row, expected)
        vectors.append(row)
    norm = [sum((vectors[j][i] for j in range(3)), O) for i in range(3)]
    assert norm == [target, target, target]
    assert all(Pt != O and r * Pt == O for Pt in norm)
    assert [sigma(sigma(sigma(Pt))) for Pt in positive + negative] == positive + negative

    def encode_point(Pt):
        return None if Pt == O else [str(Pt[0]), str(Pt[1])]

    def encode_element(v):
        return [str(c) for c in v.polynomial().list()]

    def encode_cover_point(Pt):
        return [encode_element(Pt[0]), [encode_element(y) for y in Pt[1]]]

    result = {
        "curve": name,
        "base_prime_decimal": str(p),
        "subgroup_order_decimal": str(r),
        "field_bits": int(p).bit_length(),
        "beta_degree_over_Fp": 3 if beta**(p**3) == beta else 6,
        "beta_Fp6_coefficients": encode_element(beta),
        "elliptic_quotient_mobius_models": len(models),
        "mobius_rho_C_D_kappa_Fp6_coefficients": [encode_element(v) for v in (rho, C, D, kappa)],
        "search_cap_m": 32,
        "attempted_split_fibers": attempted,
        "selected_m": m,
        "selected_mG": encode_point(m * G),
        "selected_2mG": encode_point(2 * m * G),
        "fiber_mG_Fp6_coefficients": [encode_cover_point(Pt) for Pt in positive],
        "fiber_2mG_Fp6_coefficients": [encode_cover_point(Pt) for Pt in negative],
        "expected_component_minus_2mG": encode_point(target),
        "three_quotient_vectors": [[encode_point(Pt) for Pt in row] for row in vectors],
        "norm_quotient_vector": [encode_point(Pt) for Pt in norm],
        "three_step_norm_fixed_by_q_frobenius": True,
        "three_step_norm_nonzero_exact_order_r": True,
        "exact_cover_fibers_replayed": True,
        "quotient_map_polynomial_identity_replayed": True,
        "geometric_reduction": "degree-two push-pull and Klein-four three-quotient identity on odd torsion",
    }
    print(f"{name}: m={m}, beta degree={result['beta_degree_over_Fp']}, norm diagonal exact order r", flush=True)
    return result


rows = [run(name, params) for name, params in CURVES.items()]
out = {
    "status": "verified_larger_type_i_three_step_norm_quotients",
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "rows": rows,
}
(HERE / "type_i_larger_norm.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps({"status": out["status"], "selected": {v["curve"]: v["selected_m"] for v in rows}}, indent=2))
