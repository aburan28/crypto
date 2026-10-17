"""Exact arithmetic and complete distinct-root isolation in GF(2**n).

An integer is a polynomial-basis field element: bit i is the coefficient of
z**i.  ``BinaryField`` verifies its defining polynomial by Rabin's exact
irreducibility criterion.  Polynomial coefficients in the separate root
variable X are ascending tuples; the zero polynomial is ``()``.

Root isolation first computes R = gcd(L, X**q - X) by n modular squarings,
then splits R with absolute traces.  It neither enumerates GF(2**83) nor
computes a Groebner basis.  Its completeness certificate concerns DISTINCT
roots in this field, not multiplicities or roots over extension fields.
"""

from functools import lru_cache
import random


_MAX_EXTENSION_DEGREE = 256
_UNRESTRICTED_ENUMERATION_LIMIT = 100_000
_SQUARE_BYTES = tuple(
    sum(((byte >> bit) & 1) << (2 * bit) for bit in range(8))
    for byte in range(256)
)


def _integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer, not a float or boolean")
    return value


def _bit_remainder(a, modulus):
    """Polynomial remainder in GF(2)[z], represented by integer bitmasks."""
    width = modulus.bit_length()
    while a.bit_length() >= width:
        a ^= modulus << (a.bit_length() - width)
    return a


def _bit_gcd(a, b):
    while b:
        a, b = b, _bit_remainder(a, b)
    return a


def _bit_square(a):
    result = 0
    shift = 0
    while a:
        result |= _SQUARE_BYTES[a & 255] << shift
        a >>= 8
        shift += 16
    return result


def _prime_divisors(n):
    result = set()
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            result.add(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    if n > 1:
        result.add(n)
    return result


@lru_cache(maxsize=512)
def is_irreducible(modulus):
    """Apply Rabin's exact criterion to a monic GF(2) polynomial bitmask.

    Every positive integer bitmask is monic at its highest set bit.  Degree
    zero polynomials are not irreducible.  This is a proof, not a probable
    irreducibility test.
    """
    _integer(modulus, "modulus")
    if modulus < 2:
        return False
    degree = modulus.bit_length() - 1
    if degree > _MAX_EXTENSION_DEGREE:
        raise ValueError("irreducibility certification is limited to degree 256")
    if degree == 1:
        return True
    if not modulus & 1:
        return False
    # A higher-degree irreducible has neither z nor z+1 as a factor.
    if modulus.bit_count() % 2 == 0:
        return False
    checkpoints = {degree // prime for prime in _prime_divisors(degree)}
    x = 2
    power = x
    for exponent in range(1, degree + 1):
        power = _bit_remainder(_bit_square(power), modulus)
        if exponent in checkpoints and _bit_gcd(power ^ x, modulus) != 1:
            return False
    return power == x


@lru_cache(maxsize=256)
def _find_irreducible(degree):
    if degree == 1:
        return 3
    # Deterministic generation; every chosen candidate is proved irreducible.
    rng = random.Random(0x42524F4F5453 ^ degree)
    monic = 1 << degree
    for _ in range(65_536):
        candidate = monic | rng.getrandbits(degree) | 1
        if is_irreducible(candidate):
            return candidate
    raise RuntimeError("deterministic irreducible-polynomial search budget exhausted")


class BinaryField:
    """GF(2**degree), with canonical elements 0 <= a < 2**degree.

    ``modulus`` is the defining monic, irreducible binary polynomial bitmask.
    If omitted, reproducible seeded search selects and verifies one.  This
    representation need not match a caller's external field representation;
    provide that representation's modulus when interoperability is required.
    """

    __slots__ = ("degree", "modulus", "order", "characteristic")

    def __init__(self, degree=83, modulus=None):
        _integer(degree, "degree")
        if not 1 <= degree <= _MAX_EXTENSION_DEGREE:
            raise ValueError("binary extension degree must be between 1 and 256")
        if modulus is None:
            modulus = _find_irreducible(degree)
        _integer(modulus, "modulus")
        if modulus <= 0 or modulus.bit_length() != degree + 1:
            raise ValueError("modulus must be a monic bitmask of the specified degree")
        if not is_irreducible(modulus):
            raise ValueError("modulus must be irreducible over GF(2)")
        self.degree = degree
        self.modulus = modulus
        self.order = 1 << degree
        self.characteristic = 2

    def __repr__(self):
        return f"BinaryField(degree={self.degree}, modulus={hex(self.modulus)})"

    def validate(self, value):
        """Return a canonical element, rejecting silent truncation/coercion."""
        _integer(value, "field element")
        if not 0 <= value < self.order:
            raise ValueError(f"field element must be between 0 and {self.order - 1}")
        return value

    def add(self, a, b):
        return self.validate(a) ^ self.validate(b)

    sub = add

    def neg(self, a):
        return self.validate(a)

    def _raw_mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        if a == 1:
            return b
        if b == 1:
            return a
        if a.bit_length() < b.bit_length():
            a, b = b, a
        product = 0
        while b:
            if b & 1:
                product ^= a
            b >>= 1
            a <<= 1
            if a & self.order:
                a ^= self.modulus
        return product

    def mul(self, a, b):
        return self._raw_mul(self.validate(a), self.validate(b))

    def _raw_square(self, a):
        return _bit_remainder(_bit_square(a), self.modulus)

    def square(self, a):
        return self._raw_square(self.validate(a))

    def pow(self, a, exponent):
        a = self.validate(a)
        _integer(exponent, "exponent")
        if exponent < 0:
            a = self.inv(a)
            exponent = -exponent
        result = 1
        while exponent:
            if exponent & 1:
                result = self._raw_mul(result, a)
            exponent >>= 1
            if exponent:
                a = self._raw_square(a)
        return result

    def inv(self, a):
        a = self.validate(a)
        if a == 0:
            raise ZeroDivisionError("zero has no multiplicative inverse")
        # Extended Euclid in GF(2)[z], without exponentiating by q-2.
        u, v = a, self.modulus
        left, right = 1, 0
        while u != 1:
            shift = u.bit_length() - v.bit_length()
            if shift < 0:
                u, v = v, u
                left, right = right, left
                shift = -shift
            u ^= v << shift
            left ^= right << shift
        return _bit_remainder(left, self.modulus)


def _field(field):
    if not isinstance(field, BinaryField):
        raise TypeError("field must be a BinaryField")
    return field


def _trim(values):
    while values and not values[-1]:
        values.pop()
    return tuple(values)


def trim(coeffs, field=None):
    """Validate ascending integer coefficients and return a trimmed tuple."""
    if isinstance(coeffs, (str, bytes, dict)):
        raise TypeError("coefficients must be an iterable of integers")
    try:
        values = list(coeffs)
    except TypeError as exc:
        raise TypeError("coefficients must be an iterable of integers") from exc
    if field is not None:
        _field(field)
    for value in values:
        _integer(value, "coefficient")
        if value < 0:
            raise ValueError("binary polynomial coefficients must be nonnegative")
        if field is not None:
            field.validate(value)
    return _trim(values)


def _add(a, b):
    return _trim([(a[i] if i < len(a) else 0) ^
                  (b[i] if i < len(b) else 0)
                  for i in range(max(len(a), len(b)))])


def add(a, b, field):
    _field(field)
    return _add(trim(a, field), trim(b, field))


sub = add


def _mul(a, b, field):
    if not a or not b:
        return ()
    product = [0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        if left:
            for j, right in enumerate(b):
                if right:
                    product[i + j] ^= field._raw_mul(left, right)
    return _trim(product)


def mul(a, b, field):
    _field(field)
    return _mul(trim(a, field), trim(b, field), field)


def _divmod(a, b, field):
    if not b:
        raise ZeroDivisionError("division by the zero polynomial")
    if len(a) < len(b):
        return (), a
    remainder = list(a)
    quotient = [0] * (len(a) - len(b) + 1)
    inverse_leading = field.inv(b[-1]) if b[-1] != 1 else 1
    while len(remainder) >= len(b):
        shift = len(remainder) - len(b)
        coefficient = field._raw_mul(remainder[-1], inverse_leading)
        quotient[shift] = coefficient
        for j, value in enumerate(b):
            if value:
                remainder[shift + j] ^= field._raw_mul(coefficient, value)
        while remainder and not remainder[-1]:
            remainder.pop()
    return _trim(quotient), tuple(remainder)


def divmod_poly(a, b, field):
    _field(field)
    return _divmod(trim(a, field), trim(b, field), field)


def _monic(poly, field):
    if not poly or poly[-1] == 1:
        return poly
    inverse_leading = field.inv(poly[-1])
    return tuple(field._raw_mul(value, inverse_leading) for value in poly)


def _gcd(a, b, field):
    while b:
        a, b = b, _divmod(a, b, field)[1]
    return _monic(a, field)


def gcd_poly(a, b, field):
    _field(field)
    return _gcd(trim(a, field), trim(b, field), field)


def _squaremod(poly, modulus, field):
    if not poly:
        return ()
    square = [0] * (2 * len(poly) - 1)
    for i, value in enumerate(poly):
        square[2 * i] = field._raw_square(value)
    return _divmod(_trim(square), modulus, field)[1]


def powmod_poly(base, exponent, modulus, field):
    """Raise a polynomial modulo a nonzero polynomial; exponent is nonnegative."""
    _field(field)
    _integer(exponent, "exponent")
    if exponent < 0:
        raise ValueError("exponent must be nonnegative")
    base, modulus = trim(base, field), trim(modulus, field)
    if not modulus:
        raise ZeroDivisionError("modulus is the zero polynomial")
    base = _divmod(base, modulus, field)[1]
    result = _divmod((1,), modulus, field)[1]
    while exponent:
        if exponent & 1:
            result = _divmod(_mul(result, base, field), modulus, field)[1]
        exponent >>= 1
        if exponent:
            base = _squaremod(base, modulus, field)
    return result


def evaluate_poly(coeffs, x, field):
    _field(field)
    x = field.validate(x)
    result = 0
    for value in reversed(trim(coeffs, field)):
        result = field._raw_mul(result, x) ^ value
    return result


def binary_field_roots(coeffs, field, seed=0, max_trials=128):
    """Find and exactly certify all distinct roots in a binary field.

    A COMPLETE result verifies each root against the original polynomial and
    verifies that their count equals deg gcd(L, X**q-X), which is squarefree.
    A global trace-splitting budget is applied across all factors.  Exhaustion
    reports UNRESOLVED and returns only individually verified partial roots.
    The zero polynomial reports UNRESTRICTED, materializing its domain only
    for q <= 100000.  All ordinary polynomial outputs are ascending tuples;
    this result's sorted ``roots`` is a list suitable for JSON serialization.
    """
    _field(field)
    _integer(seed, "seed")
    _integer(max_trials, "max_trials")
    if max_trials < 0:
        raise ValueError("max_trials must be nonnegative")
    polynomial = trim(coeffs, field)
    degree = len(polynomial) - 1
    if not polynomial:
        enumerated = field.order <= _UNRESTRICTED_ENUMERATION_LIMIT
        return {
            "roots": list(range(field.order)) if enumerated else [],
            "status": "unrestricted",
            "degree": -1,
            "restricted_degree": field.order,
            "domain_size": field.order,
            "roots_enumerated": enumerated,
            "split_trials": 0,
            "reason": "L is zero over GF(2^n), so every field element is allowed.",
        }
    if degree == 0:
        restricted = (1,)
    else:
        frobenius = _divmod((0, 1), polynomial, field)[1]
        for _ in range(field.degree):
            frobenius = _squaremod(frobenius, polynomial, field)
        restricted = _gcd(polynomial, _add(frobenius, (0, 1)), field)
    restricted_degree = len(restricted) - 1
    roots = set()
    # Deterministically extract zero; R is squarefree even when L is not.
    if restricted_degree and restricted[0] == 0:
        roots.add(0)
        restricted = restricted[1:]
    stack = [restricted] if len(restricted) > 1 else []
    rng = random.Random(seed)
    trials = 0
    while stack:
        factor = stack.pop()
        if len(factor) == 2:
            roots.add(field._raw_mul(factor[0], field.inv(factor[1])))
            continue
        while trials < max_trials:
            trials += 1
            scalar = rng.randrange(field.order)
            power = _divmod((0, scalar), factor, field)[1]
            trace = ()
            for iteration in range(field.degree):
                trace = _add(trace, power)
                if iteration + 1 < field.degree:
                    power = _squaremod(power, factor, field)
            divisor = _gcd(trace, factor, field)
            if 1 < len(divisor) < len(factor):
                quotient, remainder = _divmod(factor, divisor, field)
                if remainder:
                    raise ArithmeticError("internal trace split did not divide exactly")
                stack.extend((divisor, _monic(quotient, field)))
                break
        # Pending nonlinear factors remain unresolved when the budget is spent.
    verified = sorted(root for root in roots
                      if evaluate_poly(polynomial, root, field) == 0)
    complete = len(verified) == restricted_degree
    return {
        "roots": verified,
        "status": "complete" if complete else "unresolved",
        "degree": degree,
        "restricted_degree": restricted_degree,
        "domain_size": field.order,
        "roots_enumerated": complete,
        "split_trials": trials,
        "reason": (
            f"All {restricted_degree} distinct roots of gcd(L, X^q-X) verified exactly."
            if complete else
            f"Splitting budget exhausted: {len(verified)} of {restricted_degree} roots verified."
        ),
    }
