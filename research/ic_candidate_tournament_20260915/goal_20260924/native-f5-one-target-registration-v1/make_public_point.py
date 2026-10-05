"""Deterministic public-point fixture; never computes or stores its logarithm.

Run only with /Volumes/SSD990/cryptanalysis/sage -python, the checked launcher.
Fixture construction precedes either target-dependent timing interval.
"""
import hashlib
import json

from sage.all import EllipticCurve, GF, PolynomialRing

bits = 17
order = 65587
seed = b"native-f5-one-target-20261005-public-point-v1"
F2 = GF(2)
T = PolynomialRing(F2, "t").gen()
F = GF(2**bits, "z", modulus=T**17 + T**3 + 1)
z = F.gen()
E = EllipticCurve(F, [1, 1, 0, 0, 1])


def encoded(value):
    return sum(int(coefficient) << i for i, coefficient in enumerate(value.polynomial().list()))


def element(value):
    return sum(((value >> i) & 1) * z**i for i in range(bits))


start_x = int.from_bytes(hashlib.sha256(seed).digest()[:4], "big") % (1 << bits)
for offset in range(1 << bits):
    x_bits = (start_x + offset) % (1 << bits)
    try:
        candidate = E.lift_x(element(x_bits))
    except ValueError:
        continue
    public = 2 * candidate
    if public.is_zero() or not (order * public).is_zero():
        continue
    print(json.dumps({"schema_version": 1, "seed": seed.decode(),
        "start_x": start_x, "offset": offset, "lift_x": x_bits,
        "target": [encoded(public[0]), encoded(public[1])],
        "subgroup_order": order, "scalar_generated": False}, sort_keys=True))
    break
else:
    raise RuntimeError("no nonidentity subgroup point found")
