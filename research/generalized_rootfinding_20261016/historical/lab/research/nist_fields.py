"""Exact NIST prime-field presets using FLINT's proved-prime interface.

This separate backend leaves the frozen Pocklington-only field API intact.
FLINT fmpz_is_prime returns one only for a proof (Pocklington/Morrison/BLS
or APRCL), unlike fmpz_is_probabprime. Its binary identity is recorded.
No field-size or primality check constitutes a hard root-solving result.
"""

import ctypes
import ctypes.util
import hashlib
import json
from math import isqrt
from pathlib import Path

from rootfinder.fields import PrimeField

ROOT = Path(__file__).resolve().parents[2]
PRESETS = Path(__file__).with_name("nist_presets.json")
PROOF_DOCUMENTATION = "https://flintlib.org/doc/fmpz.html#c.fmpz_is_prime"


def _library_path(supplied=None):
    if supplied is not None:
        candidate = Path(supplied).resolve()
        if not candidate.is_file():
            raise FileNotFoundError("FLINT proof library was not found")
        return candidate
    candidates = [ROOT / "vendor/msolve-local/usr/lib/x86_64-linux-gnu/libflint.so.24.0.0"]
    candidates.extend(sorted((ROOT / "vendor/python-flint").glob("**/libflint*.so*")))
    for path in candidates:
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError("Install an exact FLINT backend or supply a FLINT library path")


def prove_prime(integer, flint_library=None):
    """Use the documented deterministic proof entry point; retain provenance."""
    if isinstance(integer, bool) or not isinstance(integer, int) or integer < 2:
        raise ValueError("Prime candidate must be an integer at least two")
    if ctypes.sizeof(ctypes.c_long) != 8:
        raise RuntimeError("This version of the FLINT binding requires a 64-bit long ABI")
    path = _library_path(flint_library)
    library = ctypes.CDLL(str(path))
    pointer = ctypes.POINTER(ctypes.c_long)
    library.fmpz_init.argtypes = [pointer]
    library.fmpz_init.restype = None
    library.fmpz_set_str.argtypes = [pointer, ctypes.c_char_p, ctypes.c_int]
    library.fmpz_set_str.restype = ctypes.c_int
    library.fmpz_is_prime.argtypes = [pointer]
    library.fmpz_is_prime.restype = ctypes.c_int
    library.fmpz_clear.argtypes = [pointer]
    library.fmpz_clear.restype = None
    value = ctypes.c_long()
    library.fmpz_init(ctypes.byref(value))
    try:
        if library.fmpz_set_str(ctypes.byref(value), str(integer).encode("ascii"), 10):
            raise ValueError("FLINT rejected the prime candidate encoding")
        outcome = library.fmpz_is_prime(ctypes.byref(value))
    finally:
        library.fmpz_clear(ctypes.byref(value))
    if outcome != 1:
        raise ValueError("Exact FLINT primality proof rejected the candidate")
    return {"method": "FLINT fmpz_is_prime deterministic proof",
            "integer_hex": hex(integer), "proven_prime": True,
            "library_path": str(path),
            "library_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "documentation": PROOF_DOCUMENTATION,
            "portable_pocklington_certificate_emitted": False}


class CertifiedPrimeField(PrimeField):
    """PrimeField arithmetic after an exact source-identified FLINT proof."""

    __slots__ = ("prime_proof",)

    def __init__(self, characteristic, flint_library=None):
        self.prime_proof = prove_prime(characteristic, flint_library)
        self.characteristic = self.order = characteristic
        self.degree = 1
        self.modulus = None


def preset(name):
    """Return copied standard parameters; accept P192 and P-192 spellings."""
    canonical = str(name).upper().replace("-", "")
    payload = json.loads(PRESETS.read_text())
    record = next((item for item in payload["curves"]
                   if item["name"].replace("-", "") == canonical), None)
    if record is None:
        raise ValueError("Choose P-192, P-224, P-256, P-384 or P-521")
    result = dict(record)
    for key in ("p", "a", "b", "gx", "gy", "n"):
        result[key] = int(record[key + "_hex"], 16)
    result["curve"] = {"a2": 0, "a4": result["a"], "a6": result["b"]}
    result["generator"] = (result["gx"], result["gy"])
    return result


def build_nist_target(name, scalar=19, flint_library=None):
    """Charge exact p,n proofs, G/nG, Hasse uniqueness and ordinary target."""
    from lab.research.odd_solver import on_curve, scalar_multiply
    if isinstance(scalar, bool) or not isinstance(scalar, int) or scalar <= 0:
        raise ValueError("Public target scalar must be a positive integer")
    record = preset(name)
    field = CertifiedPrimeField(record["p"], flint_library)
    order_proof = prove_prime(record["n"], flint_library)
    curve, generator = record["curve"], record["generator"]
    if record["h"] != 1 or record["a"] != record["p"] - 3:
        raise AssertionError("NIST preset cofactor or a coefficient changed")
    if not on_curve(field, curve, generator):
        raise AssertionError("Standard generator is not on the curve")
    if scalar_multiply(field, curve, generator, record["n"]) is not None:
        raise AssertionError("Standard prime order does not kill the generator")
    # A nonidentity point of prime order n makes n divide #E. There is a
    # unique multiple of n in the conservative integer Hasse interval.
    radius = isqrt(4 * record["p"]) + 1
    low, high = record["p"] + 1 - radius, record["p"] + 1 + radius
    if (low + record["n"] - 1) // record["n"] != high // record["n"]:
        raise AssertionError("Hasse interval does not certify a unique order multiple")
    if not low <= record["n"] <= high:
        raise AssertionError("Standard group order is outside the Hasse interval")
    target = scalar_multiply(field, curve, generator, scalar)
    if target is None or not on_curve(field, curve, target):
        raise ValueError("Public scalar gives no finite target")
    return {"field": field, "curve": curve, "generator": generator,
            "target": target, "target_scalar": scalar, "preset": record,
            "proof": {"field": field.prime_proof, "group_order": order_proof,
                      "generator_membership_verified": True, "nG_identity_verified": True,
                      "group_order_equal_n_by_Hasse_verified": True,
                      "rational_two_torsion_absent_verified": bool(record["n"] & 1),
                      "target_generated_from_public_scalar": scalar,
                      "target_generated_from_source_tuple": False}}
