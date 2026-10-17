"""Independent replay of frozen P-224/P-256 Type I norm point coordinates.

The verifier reads exact F_(p^6) coefficients from the producer receipt. It
reconstructs the quartic, cover fibers, three quotient maps, q-Frobenius
cycle, first-split search, and subgroup arithmetic without importing the
producer or using its beta solver.
"""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, PolynomialRing, prod


HERE = Path(__file__).resolve().parent
producer_path = HERE / "type_i_larger_norm.json"
receipt = json.loads(producer_path.read_text())
assert receipt["status"] == "verified_larger_type_i_three_step_norm_quotients"
assert receipt["source_sha256"] == hashlib.sha256((HERE / "type_i_larger_norm.py").read_bytes()).hexdigest()


def verify(row):
    p = Integer(row["base_prime_decimal"])
    r = Integer(row["subgroup_order_decimal"])
    q = p**2
    F = GF(p)
    L = GF(p**6, "u")
    zL = L.gen()
    R = PolynomialRing(F, "x")
    x = R.gen()
    a = p - 3
    b = {
        "P-224": Integer("B4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4", 16),
        "P-256": Integer("5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B", 16),
    }[row["curve"]]
    generator = {
        "P-224": ("B70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21", "BD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34"),
        "P-256": ("6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5"),
    }[row["curve"]]
    f = x**3 + F(a) * x + F(b)
    assert f.is_irreducible()
    roots = f.change_ring(L).roots(multiplicities=False)
    assert len(roots) == 3

    def dec(c):
        return sum((L(Integer(v)) * zL**i for i, v in enumerate(c)), L(0))

    beta = dec(row["beta_Fp6_coefficients"])
    valid_epsilon = []
    for eps in roots:
        lam = (eps**(p**2) - eps) / (eps**p - eps)
        eps_q = eps**q
        if len({eps, eps_q, beta, beta**q}) != 4:
            continue
        replay = ((beta - eps_q) / (beta - eps)) * ((beta**q - eps) / (beta**q - eps_q))
        if replay == lam:
            valid_epsilon.append(eps)
    assert valid_epsilon
    eps = valid_epsilon[0]
    quartic_roots = [eps, eps**q, beta, beta**q]
    P = PolynomialRing(L, "X")
    X = P.gen()
    Q = prod(X - v for v in quartic_roots)
    assert Q.is_squarefree()
    rho, C, D, kappa = [dec(c) for c in row["mobius_rho_C_D_kappa_Fp6_coefficients"]]
    assert rho in quartic_roots and kappa != 0 and D != 0
    linear = C * (X - rho) + D
    cleared = linear**3 * (X - rho) + L(a) * linear * (X - rho)**3 + L(b) * (X - rho)**4
    assert cleared == kappa**2 * Q

    E = EllipticCurve(L, [L(a), L(b)])
    gx, gy = [Integer(v, 16) for v in generator]
    G = E(L(gx), L(gy))
    O = E(0)
    assert r.is_prime() and r * G == O and G != O
    m = Integer(row["selected_m"])
    assert 1 <= m <= row["search_cap_m"]
    target = -2 * m * G
    assert target != O and r * target == O

    def decode_cover(v):
        return (dec(v[0]), tuple(dec(c) for c in v[1]))

    positive = [decode_cover(v) for v in row["fiber_mG_Fp6_coefficients"]]
    negative = [decode_cover(v) for v in row["fiber_2mG_Fp6_coefficients"]]
    assert len(positive) == len(negative) == 2

    def check_cover(pt):
        x0, ys = pt
        for i in range(3):
            assert ys[i]**2 == prod(x0 - v**(q**i) for v in quartic_roots)
        expected_norm = prod((x0 - eps**(q**i)) * (x0 - beta**(q**i)) for i in range(3))
        assert prod(ys) == expected_norm

    def project(pt, i):
        x0, ys = pt
        rhoi, Ci, Di, ki = [v**(q**i) for v in (rho, C, D, kappa)]
        assert x0 != rhoi
        return E(Ci + Di / (x0 - rhoi), ki * ys[i] / (x0 - rhoi)**2)

    for pt in positive + negative:
        check_cover(pt)
    assert [project(pt, 0) for pt in positive] == [m * G, m * G]
    assert [project(pt, 0) for pt in negative] == [2 * m * G, 2 * m * G]
    for pair in (positive, negative):
        assert pair[0][0] == pair[1][0]
        assert pair[0][1][0] == pair[1][1][0]
        assert all(pair[0][1][i] == -pair[1][1][i] for i in (1, 2))

    # Independently replay every candidate preceding the selected multiple.
    def splits(Pt):
        if Pt == O or Pt[0] == C:
            return False
        x0 = rho + D / (Pt[0] - C)
        rhs1 = prod(x0 - v**q for v in quartic_roots)
        return rhs1 != 0 and rhs1.is_square()

    for item in row["attempted_split_fibers"]:
        i = Integer(item["m"])
        assert item["fiber_mG_split"] == splits(i * G)
        assert item["fiber_2mG_split"] == splits(2 * i * G)
        if i < m:
            assert not (item["fiber_mG_split"] and item["fiber_2mG_split"])
    assert row["attempted_split_fibers"][-1]["m"] == m

    def sigma(pt):
        x0, (y0, y1, y2) = pt
        return (x0**q, (y2**q, y0**q, y1**q))

    def enc(Pt):
        return None if Pt == O else [str(Pt[0]), str(Pt[1])]

    vectors = []
    for j in range(3):
        pos = positive
        neg = negative
        for _ in range(j):
            pos = [sigma(pt) for pt in pos]
            neg = [sigma(pt) for pt in neg]
        for pt in pos + neg:
            check_cover(pt)
        vector = [sum((project(pt, i) for pt in pos), O) -
                  sum((project(pt, i) for pt in neg), O) for i in range(3)]
        assert vector[j] == target
        assert all(vector[i] == O for i in range(3) if i != j)
        assert [enc(pt) for pt in vector] == row["three_quotient_vectors"][j]
        vectors.append(vector)
    norm = [sum((vectors[j][i] for j in range(3)), O) for i in range(3)]
    assert norm == [target, target, target]
    assert [enc(pt) for pt in norm] == row["norm_quotient_vector"]
    assert all(pt != O and r * pt == O for pt in norm)
    assert [sigma(sigma(sigma(pt))) for pt in positive + negative] == positive + negative
    return {"curve": row["curve"], "selected_m": int(m), "all_nine_quotient_components_replayed": True,
            "first_split_policy_replayed": True, "norm_diagonal_nonzero_exact_order_r": True}


results = [verify(row) for row in receipt["rows"]]
out = {"status": "independent_larger_type_i_norm_replay_passed",
       "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       "producer_receipt_sha256": hashlib.sha256(producer_path.read_bytes()).hexdigest(),
       "results": results}
(HERE / "type_i_larger_norm_verification.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps(out, indent=2, sort_keys=True))
