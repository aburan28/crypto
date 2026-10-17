"""Exact prime and polynomial-basis extension fields.

Elements are canonical integers.  In GF(p**n), integer ``a`` represents the
ascending base-p digits of a polynomial in the defining generator.  This is
an encoding, not arithmetic modulo p**n.  Binary fields retain their existing
optimized bitmask implementation.

Characteristics below 2**64 use a deterministic primality test.  Larger
characteristics require an exact Pocklington certificate of the form::

    {"n": p, "factors": [{"q": q, "e": e, "a": witness,
                          "certificate": recursive_certificate}, ...]}

``n`` is optional but must match p when present.  Distinct q are proved prime
(recursively when q >= 2**64), their product F = product(q**e) must divide
p-1 and satisfy F**2 > p.  Each witness proves a**(p-1)=1 mod p and
gcd(a**((p-1)/q)-1,p)=1.  Remaining factors of p-1 need not be known.
No probable-prime result is promoted to a field certificate.
"""

from math import gcd
import random

from .binaryfield import BinaryField, _integer, _prime_divisors
from .finitefield import is_prime


_MAX_SMALL_PRIME = (1 << 64) - 1
_MAX_EXTENSION_DEGREE = 256


def verify_prime_certificate(n, certificate=None, _depth=0):
    """Prove n prime or raise ValueError; return True only after an exact proof."""
    _integer(n, "characteristic")
    if n <= _MAX_SMALL_PRIME:
        if not is_prime(n):
            raise ValueError("characteristic must be prime")
        return True
    if _depth > 128:
        raise ValueError("prime certificate nesting exceeds 128 levels")
    if not isinstance(certificate, dict):
        raise ValueError("characteristics >= 2**64 require a Pocklington prime_certificate")
    if "n" in certificate:
        _integer(certificate["n"], "certificate n")
        if certificate["n"] != n:
            raise ValueError("prime certificate n does not match characteristic")
    factors = certificate.get("factors")
    if not isinstance(factors, (list, tuple)) or not factors:
        raise ValueError("prime certificate must supply nonempty factors")
    seen = set()
    known_product = 1
    for factor in factors:
        if not isinstance(factor, dict):
            raise ValueError("prime certificate factors must be objects")
        q, e, witness = factor.get("q"), factor.get("e"), factor.get("a")
        _integer(q, "certificate factor q")
        _integer(e, "certificate factor exponent e")
        _integer(witness, "certificate witness a")
        if q < 2 or q >= n or q in seen:
            raise ValueError("certificate q must be distinct and between 2 and n-1")
        if e < 1 or e > n.bit_length():
            raise ValueError("certificate exponents must be positive and fit n-1")
        if not 1 <= witness < n:
            raise ValueError("certificate witness must be between 1 and n-1")
        verify_prime_certificate(q, factor.get("certificate"), _depth + 1)
        seen.add(q)
        known_product *= q ** e
        if (n - 1) % known_product:
            raise ValueError("certificate factor product must divide n-1")
        if pow(witness, n - 1, n) != 1:
            raise ValueError("Pocklington Fermat witness check failed")
        if gcd(pow(witness, (n - 1) // q, n) - 1, n) != 1:
            raise ValueError("Pocklington gcd witness check failed")
    if known_product * known_product <= n:
        raise ValueError("Pocklington known factor product must exceed sqrt(n)")
    return True


class PrimeField:
    """GF(p), with canonical element integers 0 <= a < p."""

    __slots__ = ("characteristic", "degree", "order", "modulus")

    def __init__(self, p, prime_certificate=None):
        verify_prime_certificate(p, prime_certificate)
        self.characteristic = self.order = p
        self.degree = 1
        self.modulus = None

    def __repr__(self):
        return f"PrimeField({self.characteristic})"

    def validate(self, value):
        _integer(value, "field element")
        if not 0 <= value < self.order:
            raise ValueError(f"field element must be between 0 and {self.order - 1}")
        return value

    def from_int(self, value):
        return _integer(value, "integer") % self.order

    def coefficients(self, value):
        return (self.validate(value),)

    def element(self, coefficients):
        values = list(coefficients)
        if len(values) > 1:
            raise ValueError("prime-field elements have at most one coefficient")
        return self.validate(values[0]) if values else 0

    def _raw_add(self, a, b):
        return (a + b) % self.order

    def _raw_sub(self, a, b):
        return (a - b) % self.order

    def _raw_mul(self, a, b):
        return a * b % self.order

    def add(self, a, b):
        return self._raw_add(self.validate(a), self.validate(b))

    def sub(self, a, b):
        return self._raw_sub(self.validate(a), self.validate(b))

    def neg(self, a):
        return -self.validate(a) % self.order

    def mul(self, a, b):
        return self._raw_mul(self.validate(a), self.validate(b))

    def square(self, a):
        a = self.validate(a)
        return self._raw_mul(a, a)

    def pow(self, a, exponent):
        a = self.validate(a)
        _integer(exponent, "exponent")
        if exponent < 0 and a == 0:
            raise ZeroDivisionError("zero has no multiplicative inverse")
        return pow(a, exponent, self.order)

    def inv(self, a):
        a = self.validate(a)
        if a == 0:
            raise ZeroDivisionError("zero has no multiplicative inverse")
        return pow(a, -1, self.order)


def _canonical_modulus(modulus, p, degree=None):
    if isinstance(modulus, (str, bytes, dict)):
        raise TypeError("extension modulus must contain ascending prime-field coefficients")
    try:
        values = tuple(modulus)
    except TypeError as exc:
        raise TypeError("extension modulus must contain ascending prime-field coefficients") from exc
    for coefficient in values:
        _integer(coefficient, "modulus coefficient")
        if not 0 <= coefficient < p:
            raise ValueError("modulus coefficients must be canonical prime-field elements")
    if len(values) < 2 or values[-1] != 1:
        raise ValueError("extension modulus must be monic and have positive degree")
    if degree is not None and len(values) != degree + 1:
        raise ValueError("modulus degree must match extension degree")
    return values


def _rabin_irreducible(modulus, base):
    from .polynomial import gcd_poly, powmod_poly, sub, divmod_poly

    degree = len(modulus) - 1
    if degree == 1:
        return True
    if modulus[0] == 0:
        return False
    x = divmod_poly((0, 1), modulus, base)[1]
    power = x
    checkpoints = {degree // prime for prime in _prime_divisors(degree)}
    for exponent in range(1, degree + 1):
        power = powmod_poly(power, base.order, modulus, base)
        if exponent in checkpoints and gcd_poly(sub(power, x, base), modulus, base) != (1,):
            return False
    return power == x


def is_irreducible_modulus(modulus, characteristic, prime_certificate=None):
    """Verify Rabin's exact irreducibility criterion over a certified GF(p)."""
    base = PrimeField(characteristic, prime_certificate)
    modulus = _canonical_modulus(modulus, characteristic)
    if len(modulus) - 1 > _MAX_EXTENSION_DEGREE:
        raise ValueError("irreducibility certification is limited to degree 256")
    return _rabin_irreducible(modulus, base)


class ExtensionField:
    """GF(p**n), using base-p encoded polynomial coefficients, n <= 256.

    The defining modulus is proved irreducible, not merely assumed.  Omitted
    moduli are selected by bounded reproducible search.  Reusing a supplied
    modulus is necessary to match an external polynomial-basis encoding.
    """

    __slots__ = ("characteristic", "degree", "order", "modulus", "_base")

    def __init__(self, p, degree, modulus=None, prime_certificate=None,
                 modulus_search_budget=65_536):
        _integer(degree, "degree")
        if not 1 <= degree <= _MAX_EXTENSION_DEGREE:
            raise ValueError("extension degree must be between 1 and 256")
        _integer(modulus_search_budget, "modulus_search_budget")
        if modulus_search_budget < 0:
            raise ValueError("modulus_search_budget must be nonnegative")
        self._base = PrimeField(p, prime_certificate)
        self.characteristic = p
        self.degree = degree
        self.order = p ** degree
        if modulus is None:
            if degree == 1:
                modulus = (0, 1)
            else:
                rng = random.Random(0x4649454C4453 ^ p ^ (degree << 32))
                for _ in range(modulus_search_budget):
                    candidate = (rng.randrange(1, p),) + tuple(
                        rng.randrange(p) for _ in range(degree - 1)) + (1,)
                    if _rabin_irreducible(candidate, self._base):
                        modulus = candidate
                        break
                if modulus is None:
                    raise RuntimeError("deterministic irreducible-modulus search budget exhausted")
        modulus = _canonical_modulus(modulus, p, degree)
        if not _rabin_irreducible(modulus, self._base):
            raise ValueError("extension modulus must be irreducible over GF(p)")
        self.modulus = modulus

    def __repr__(self):
        return (f"ExtensionField({self.characteristic}, degree={self.degree}, "
                f"modulus={self.modulus!r})")

    def validate(self, value):
        _integer(value, "field element")
        if not 0 <= value < self.order:
            raise ValueError(f"field element must be between 0 and {self.order - 1}")
        return value

    def from_int(self, value):
        """Embed an integer in the prime subfield (distinct from digit encoding)."""
        return _integer(value, "integer") % self.characteristic

    def _digits(self, value):
        values = []
        for _ in range(self.degree):
            value, digit = divmod(value, self.characteristic)
            values.append(digit)
        return values

    def _encode(self, values):
        result = 0
        for value in reversed(values):
            result = result * self.characteristic + value
        return result

    def coefficients(self, value):
        return tuple(self._digits(self.validate(value)))

    def element(self, coefficients):
        if isinstance(coefficients, (str, bytes, dict)):
            raise TypeError("element coefficients must be an iterable of integers")
        values = list(coefficients)
        if len(values) > self.degree:
            raise ValueError("element coefficient count exceeds extension degree")
        for value in values:
            self._base.validate(value)
        return self._encode(values)

    def _raw_add(self, a, b):
        p = self.characteristic
        result = 0
        place = 1
        while a or b:
            a, left = divmod(a, p)
            b, right = divmod(b, p)
            result += (left + right) % p * place
            place *= p
        return result

    def _raw_sub(self, a, b):
        p = self.characteristic
        result = 0
        place = 1
        while a or b:
            a, left = divmod(a, p)
            b, right = divmod(b, p)
            result += (left - right) % p * place
            place *= p
        return result

    def _raw_mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        if a == 1:
            return b
        if b == 1:
            return a
        p, n = self.characteristic, self.degree
        left, right = self._digits(a), self._digits(b)
        product = [0] * (2 * n - 1)
        for i, coefficient in enumerate(left):
            if coefficient:
                for j, other in enumerate(right):
                    if other:
                        product[i + j] = (product[i + j] + coefficient * other) % p
        for i in range(2 * n - 2, n - 1, -1):
            coefficient = product[i]
            if coefficient:
                for j in range(n):
                    product[i - n + j] = (
                        product[i - n + j] - coefficient * self.modulus[j]) % p
        return self._encode(product[:n])

    def add(self, a, b):
        return self._raw_add(self.validate(a), self.validate(b))

    def sub(self, a, b):
        return self._raw_sub(self.validate(a), self.validate(b))

    def neg(self, a):
        return self._raw_sub(0, self.validate(a))

    def mul(self, a, b):
        return self._raw_mul(self.validate(a), self.validate(b))

    def square(self, a):
        a = self.validate(a)
        return self._raw_mul(a, a)

    def pow(self, a, exponent):
        a = self.validate(a)
        _integer(exponent, "exponent")
        if exponent < 0:
            a, exponent = self.inv(a), -exponent
        result = 1
        while exponent:
            if exponent & 1:
                result = self._raw_mul(result, a)
            exponent >>= 1
            if exponent:
                a = self._raw_mul(a, a)
        return result

    def inv(self, a):
        from .polynomial import divmod_poly, mul, sub

        a = self.validate(a)
        if a == 0:
            raise ZeroDivisionError("zero has no multiplicative inverse")
        previous, current = self.modulus, tuple(self._digits(a))
        previous_coefficient, current_coefficient = (), (1,)
        while any(current):
            quotient, remainder = divmod_poly(previous, current, self._base)
            previous, current = current, remainder
            previous_coefficient, current_coefficient = current_coefficient, sub(
                previous_coefficient, mul(quotient, current_coefficient, self._base), self._base)
        inverse_leading = self._base.inv(previous[0])
        inverse = tuple(value * inverse_leading % self.characteristic
                        for value in previous_coefficient)
        inverse = divmod_poly(inverse, self.modulus, self._base)[1]
        return self._encode(inverse)


def make_field(spec=None):
    """Build a certified finite field from an explicit JSON-compatible spec.

    The original {degree, modulus-bitmask} form remains a binary-field spec.
    New specs use kind prime, extension, or binary; extension modulus arrays
    contain ascending prime-field coefficients.  Characteristic-two extension
    specs use the optimized BinaryField implementation.
    """
    if spec is None:
        spec = {}
    if not isinstance(spec, dict):
        raise TypeError("field specification must be an object")
    kind = spec.get("kind")
    if kind is None:
        kind = "extension" if "characteristic" in spec else "binary"
    if kind == "prime":
        return PrimeField(spec.get("characteristic"), spec.get("prime_certificate"))
    if kind == "binary":
        if "characteristic" in spec:
            _integer(spec["characteristic"], "characteristic")
            if spec["characteristic"] != 2:
                raise ValueError("binary fields have characteristic two")
        modulus = spec.get("modulus")
        if isinstance(modulus, str):
            modulus = int(modulus, 16) if modulus.lower().startswith("0x") else int(modulus, 10)
        return BinaryField(spec.get("degree", 83), modulus)
    if kind != "extension":
        raise ValueError("field kind must be prime, extension, or binary")
    p = spec.get("characteristic")
    _integer(p, "characteristic")
    degree = spec.get("degree")
    modulus = spec.get("modulus")
    if p == 2:
        if modulus is not None:
            coefficients = _canonical_modulus(modulus, 2, degree)
            modulus = sum(coefficient << i for i, coefficient in enumerate(coefficients))
        return BinaryField(degree, modulus)
    return ExtensionField(p, degree, modulus, spec.get("prime_certificate"),
                          spec.get("modulus_search_budget", 65_536))
