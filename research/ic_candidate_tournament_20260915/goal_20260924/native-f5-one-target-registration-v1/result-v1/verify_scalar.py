"""Independent Sage replay of the disclosed n17 target certificate.

Run via /Volumes/SSD990/cryptanalysis/sage -python; no target solver is run.
"""
import json

from sage.all import EllipticCurve, GF, PolynomialRing

F2 = GF(2)
t = PolynomialRing(F2, "t").gen()
F = GF(2**17, "z", modulus=t**17 + t**3 + 1)
z = F.gen()
E = EllipticCurve(F, [1, 1, 0, 0, 1])


def element(bits):
    return sum(((bits >> i) & 1) * z**i for i in range(17))


generator = E(element(43693), element(23339))
target = E(element(61889), element(74818))
order = 65587
scalar = 14668
assert not generator.is_zero() and (order * generator).is_zero()
assert not target.is_zero() and (order * target).is_zero()
assert scalar * generator == target
print(json.dumps({"schema_version": 1, "status": "PASS_INDEPENDENT_SAGE_SCALAR_REPLAY",
    "generator": [43693, 23339], "target": [61889, 74818],
    "subgroup_order": order, "recovered_scalar": scalar}, sort_keys=True))
