"""Bounded Type I cover screen for prime-field elliptic curves.

Run with /Volumes/SSD990/cryptanalysis/sage -python. The sample law is
lexicographic in (a,b), restricted to nonsingular irreducible 2-division
cubics, with at most four curves per small prime. Standard-curve parameters
are fixed below. This tests the Momose--Chao discriminant and constructs
beta; it does not estimate a relation-collection cost.
"""

import hashlib
import itertools
import json
from pathlib import Path

from sage.all import GF, Integer, PolynomialRing, matrix, prod


HERE = Path(__file__).resolve().parent
STANDARD = {
    "P-192": (
        "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF",
        "64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1",
    ),
    "P-224": (
        "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001",
        "B4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4",
    ),
    "P-256": (
        "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF",
        "5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B",
    ),
}


def nonsquare(F):
    for n in range(2, 128):
        if not F(n).is_square():
            return F(n)
    raise AssertionError("nonsquare search exhausted")


def field_coefficients(v):
    return [str(c) for c in v.polynomial().list()]


def screen(name, p, a, b, verify_beta):
    F = GF(p)
    R = PolynomialRing(F, "x")
    x = R.gen()
    f = x**3 + F(a) * x + F(b)
    assert -16 * (4 * F(a)**3 + 27 * F(b)**2) != 0
    assert f.is_irreducible()
    K3 = F.extension(f, "z")
    z = K3.gen()
    q = p**2
    roots = [z, z**p, z**(p**2)]
    assert len(set(roots)) == 3
    lam = (roots[2] - roots[0]) / (roots[1] - roots[0])
    assert lam not in (0, 1)
    j_legendre = 256 * (lam**2 - lam + 1)**3 / (lam**2 * (lam - 1)**2)
    j_cubic = F(1728 * 4 * F(a)**3 / (4 * F(a)**3 + 27 * F(b)**2))
    assert j_legendre == K3(j_cubic)
    norm = lam * lam**q * lam**(q**2)
    trace = lam + lam**q + lam**(q**2)
    assert norm == -1 and trace**p == trace
    delta = (trace - 1)**2 - 4 * norm
    assert delta**p == delta
    delta_fp = F(delta.polynomial()[0])
    assert K3(delta_fp) == delta
    d = nonsquare(F)
    W = PolynomialRing(F, "w")
    F2 = F.extension(W.gen()**2 - d, "u")
    assert F2(delta_fp).is_square()

    result = {
        "name": name,
        "p": str(p),
        "field_bits": int(p).bit_length(),
        "a": str(a),
        "b": str(b),
        "irreducible_2_division_cubic": True,
        "j_invariant_replayed": True,
        "lambda_coefficients_Fp3": field_coefficients(lam),
        "lambda_norm_Fp6_over_Fp2": str(norm),
        "lambda_trace_Fp6_over_Fp2": str(trace),
        "discriminant_factor_Fp": str(delta_fp),
        "discriminant_factor_square_Fp": bool(delta_fp.is_square()),
        "discriminant_factor_square_Fp2": True,
        "discriminant_factor_zero": bool(delta_fp == 0),
        "beta_verified": False,
    }
    if not verify_beta:
        return result

    # Paper equations (20)--(25), independently checked after solving B's
    # fixed-point quadratic. If the quadratic splits only over K6, construct
    # that exact quadratic extension rather than treating Fp3 as sufficient.
    eps = z
    eps_q = eps**q
    mu = (lam * eps - eps_q) / (lam - 1)
    nu = (lam * eps_q - eps) / (lam - 1)
    A = matrix(K3, [[nu, -eps * eps_q], [1, -mu]])
    frob = lambda M, power: matrix(K3, 2, 2, [v**power for v in M.list()])
    B = frob(A, q**2) * frob(A, q) * A
    assert B[1, 0] != 0
    disc = (B[1, 1] - B[0, 0])**2 + 4 * B[1, 0] * B[0, 1]
    if disc.is_square():
        K = K3
        root = disc.sqrt()
    else:
        T = PolynomialRing(K3, "t")
        K = K3.extension(T.gen()**2 - disc, "v")
        root = K.gen()
    assert root**2 == K(disc)
    valid = []
    for sgn in (1, -1):
        beta = (K(B[0, 0]) - K(B[1, 1]) + sgn * root) / (2 * K(B[1, 0]))
        denominator = K(A[1, 0]) * beta + K(A[1, 1])
        if denominator == 0:
            continue
        if beta**q != (K(A[0, 0]) * beta + K(A[0, 1])) / denominator:
            continue
        if len({K(eps), K(eps_q), beta, beta**q}) != 4:
            continue
        replay = ((beta - K(eps_q)) / (beta - K(eps))) * (
            (beta**q - K(eps)) / (beta**q - K(eps_q))
        )
        if replay == K(lam):
            valid.append(beta)
    assert valid
    result["beta_verified"] = True
    result["beta_count"] = len(valid)
    result["beta_field_degree_over_Fp"] = 3 if K == K3 else 6
    result["beta_field_element"] = str(valid[0])
    result["fixed_point_discriminant_square_Fp3"] = bool(disc.is_square())
    if name in STANDARD:
        beta = valid[0]
        quartic_roots = [K(eps), K(eps)**q, beta, beta**q]
        elliptic_roots = [K(eps), K(eps)**p, K(eps)**(p**2)]
        P = PolynomialRing(K, "X")
        X = P.gen()
        quartic = prod(X - v for v in quartic_roots)
        assert quartic.is_squarefree()
        models = []
        for infinity in quartic_roots:
            finite = [v for v in quartic_roots if v != infinity]
            for targets in itertools.permutations(elliptic_roots):
                reciprocal = lambda v: 1 / (v - infinity)
                D = (targets[0] - targets[1]) / (reciprocal(finite[0]) - reciprocal(finite[1]))
                C = targets[0] - D * reciprocal(finite[0])
                if C + D * reciprocal(finite[2]) != targets[2]:
                    continue
                linear = C * (X - infinity) + D
                cleared = linear**3 * (X - infinity) + K(a) * linear * (X - infinity)**3 + K(b) * (X - infinity)**4
                scale = cleared.leading_coefficient()
                assert cleared == scale * quartic
                models.append(scale)
        assert models
        result["quartic_mobius_models"] = len(models)
        result["quartic_scale_square_beta_field"] = bool(models[0].is_square())
        if K == K3:
            # A quadratic extension of Fp3 makes every Fp3 element a square.
            result["quartic_isomorphic_to_standard_over_Fp6"] = True
        else:
            result["quartic_isomorphic_to_standard_over_Fp6"] = bool(models[0].is_square())
            if not models[0].is_square():
                assert root**q == root
                twists = [root + K(c) for c in range(8) if root + K(c) != 0 and not (root + K(c)).is_square()]
                assert twists
                twist = twists[0]
                assert twist**q == twist and (models[0] / twist).is_square()
                result["required_Fp2_twist_element"] = str(twist)
                result["twisted_quartic_isomorphic_over_Fp6"] = True
    return result


rows = []
for p0 in (7, 11, 13, 17, 19, 23, 29, 31, 43):
    F0 = GF(p0)
    R0 = PolynomialRing(F0, "x")
    x0 = R0.gen()
    selected = 0
    for a0 in range(p0):
        for b0 in range(p0):
            if 4 * F0(a0)**3 + 27 * F0(b0)**2 == 0:
                continue
            if not (x0**3 + F0(a0) * x0 + F0(b0)).is_irreducible():
                continue
            rows.append(screen(f"lex-{p0}-{a0}-{b0}", Integer(p0), Integer(a0), Integer(b0), True))
            selected += 1
            if selected == 4:
                break
        if selected == 4:
            break
    assert selected == 4

for name, (p_hex, b_hex) in STANDARD.items():
    p0 = Integer(p_hex, 16)
    rows.append(screen(name, p0, p0 - 3, Integer(b_hex, 16), True))

out = {
    "status": "verified_bounded_type_i_prime_field_screen",
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "selection": "first four lexicographic nonsingular irreducible (a,b) for each listed small p; exact standardized P-192/P-224/P-256 parameters",
    "small_primes": [7, 11, 13, 17, 19, 23, 29, 31, 43],
    "rows": rows,
    "counts": {
        "total": len(rows),
        "delta_zero": sum(r["discriminant_factor_zero"] for r in rows),
        "delta_nonsquare_Fp": sum(not r["discriminant_factor_square_Fp"] for r in rows),
        "beta_verified": sum(r["beta_verified"] for r in rows),
    },
}
(HERE / "type_i_prime_field_screen.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps({"status": out["status"], "counts": out["counts"], "standard": rows[-3:]}, indent=2, sort_keys=True))
