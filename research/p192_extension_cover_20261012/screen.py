"""Exact standard NIST P-192 cover gates. Run with the cryptanalysis Sage launcher."""

import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, HyperellipticCurve, Integer, PolynomialRing


HERE = Path(__file__).resolve().parent
p = int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF", 16)
a = p - 3
b = int("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1", 16)
gx = int("188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012", 16)
gy = int("07192B95FFC8DA78631011ED6B24CDD573F977A11E794811", 16)
r = int("FFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831", 16)

start = time.monotonic()
F = GF(p)
Rx = PolynomialRing(F, "x")
x = Rx.gen()
f = x**3 + F(a) * x + F(b)
assert f.is_irreducible()
assert f.discriminant() != 0

E = EllipticCurve(F, [F(a), F(b)])
G = E(F(gx), F(gy))
assert not G.is_zero()
assert Integer(r).is_prime()
assert r * G == E(0)
assert E.cardinality() == r

# The 2-division roots live in F_(p^3).  The q=p^2 Frobenius cycles them
# in reverse order.  Every Legendre transform has norm -1 to F_(p^2).
K = F.extension(f, "z")
z = K.gen()
roots = [z, z**p, z**(p * p)]
assert len(set(roots)) == 3
assert all(rr**3 + K(a) * rr + K(b) == 0 for rr in roots)
lam = (roots[2] - roots[0]) / (roots[1] - roots[0])
orbit = [lam, 1 - lam, 1 / lam, 1 / (1 - lam),
         lam / (lam - 1), (lam - 1) / lam]
norms = []
for value in orbit:
    conjugates = [value**(p**(2 * i)) for i in range(3)]
    product = conjugates[0] * conjugates[1] * conjugates[2]
    assert product == K(-1)
    assert value.norm() == F(-1)
    norms.append(-1)

# An exact, same-field genus-2 control H -> E.  It is a genuine subgroup
# embedding, although it does not by itself descend the DLP to a smaller field.
Ru = PolynomialRing(F, "u")
u = Ru.gen()
h = u**6 + F(a) * u**2 + F(b)
assert h == f(u**2)
assert h.is_squarefree()
H = HyperellipticCurve(h)
assert H.genus() == 2
J = H.jacobian()
D = J([u**2 - F(gx), Ru(F(gy))])
assert D != J(0)
assert r * D == J(0)

# The point-count recurrence follows from the independently checked
# #E(F_p)=r above.
t = p + 1 - r
traces = {0: 2, 1: t}
orders = {}
for k in range(2, 7):
    traces[k] = t * traces[k - 1] - p * traces[k - 2]
for k in (1, 2, 3, 5, 6):
    nk = p**k + 1 - traces[k]
    orders[str(k)] = {
        "field_bits": (p**k).bit_length(),
        "order_v2": (nk & -nk).bit_length() - 1,
        "r_divides_order": nk % r == 0,
        "r_divides_multiplicative_order": (p**k - 1) % r == 0,
    }

out = {
    "curve": "NIST P-192 / secp192r1",
    "p_hex": hex(p), "a_hex": hex(a), "b_hex": hex(b),
    "r_hex": hex(r), "generator_x_hex": hex(gx), "generator_y_hex": hex(gy),
    "two_division_cubic_irreducible": True,
    "two_division_roots_in_Fp3": True,
    "r_prime": True,
    "curve_order_equals_r": True,
    "legendre_norms_Fp6_over_Fp2": norms,
    "genus2_control": {
        "equation": "v^2=(u^2)^3+a(u^2)+b",
        "map": "(u,v)->(u^2,v)",
        "degree": 2,
        "genus": 2,
        "squarefree": True,
        "generator_pullback_mumford_u": "u^2-Gx",
        "generator_pullback_mumford_v": "Gy",
        "pullback_nonzero": True,
        "r_times_pullback_zero": True,
        "pushforward_of_pullback_multiplier": 2,
        "pushforward_of_Fp6_Fp3_norm_multiplier": 4,
    },
    "extension_orders": orders,
    "elapsed_s_diagnostic": round(time.monotonic() - start, 3),
}
(HERE / "results.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print(json.dumps(out, indent=2, sort_keys=True))
