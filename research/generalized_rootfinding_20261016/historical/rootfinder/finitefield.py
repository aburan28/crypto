"""Exact polynomial arithmetic and root isolation over prime fields.

Coefficients are dense, ascending: ``(3, 2, 1)`` means 3 + 2*x + x**2.
The zero polynomial is represented by ``()``.  This module uses no Groebner
basis computation.  Root completeness is certified by the degree of
``gcd(L, x**p - x)``, computed with modular Frobenius exponentiation.

Cantor--Zassenhaus splitting is randomized, but every returned root and every
claim of completeness is checked exactly.  Exhausting its global trial budget
returns an explicit unresolved result, never a completeness claim.
"""

from functools import lru_cache
import random


_MAX_PRIME = (1 << 64) - 1
_UNRESTRICTED_ENUMERATION_LIMIT = 100_000
_MILLER_RABIN_BASES = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)


def _integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer, not a float or boolean")
    return value


@lru_cache(maxsize=256)
def _prime_test(n):
    if n < 2:
        return False
    for prime in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % prime == 0:
            return n == prime
    odd_part = n - 1
    shifts = 0
    while odd_part % 2 == 0:
        shifts += 1
        odd_part //= 2
    for base in _MILLER_RABIN_BASES:
        base %= n
        if not base:
            continue
        witness = pow(base, odd_part, n)
        if witness in (1, n - 1):
            continue
        for _ in range(shifts - 1):
            witness = witness * witness % n
            if witness == n - 1:
                break
        else:
            return False
    return True


def is_prime(n):
    """Deterministically test an integer below 2**64 for primality.

    Larger inputs are rejected; extending these seven Miller--Rabin witnesses
    past their proven range would turn exact certificates into probable ones.
    """
    _integer(n, "n")
    if n > _MAX_PRIME:
        raise ValueError("primality certification is limited to integers below 2**64")
    return _prime_test(n)


def _field(p):
    _integer(p, "p")
    if not is_prime(p):
        raise ValueError("p must be prime")
    return p


def trim(coeffs, p=None):
    """Validate coefficients, optionally reduce mod p, and remove trailing zeros."""
    if isinstance(coeffs, (str, bytes, dict)):
        raise TypeError("coefficients must be an iterable of integers")
    try:
        values = list(coeffs)
    except TypeError as exc:
        raise TypeError("coefficients must be an iterable of integers") from exc
    for value in values:
        _integer(value, "coefficient")
    if p is not None:
        _field(p)
        values = [value % p for value in values]
    while values and not values[-1]:
        values.pop()
    return tuple(values)


def add(a, b, p):
    """Add polynomials over F_p."""
    a, b = trim(a, p), trim(b, p)
    return trim(tuple(((a[i] if i < len(a) else 0)
                       + (b[i] if i < len(b) else 0)) % p
                      for i in range(max(len(a), len(b)))))


def sub(a, b, p):
    """Subtract polynomials over F_p."""
    a, b = trim(a, p), trim(b, p)
    return trim(tuple(((a[i] if i < len(a) else 0)
                       - (b[i] if i < len(b) else 0)) % p
                      for i in range(max(len(a), len(b)))))


def mul(a, b, p):
    """Multiply polynomials over F_p using dense schoolbook multiplication."""
    a, b = trim(a, p), trim(b, p)
    if not a or not b:
        return ()
    product = [0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        if not left:
            continue
        for j, right in enumerate(b):
            if right:
                product[i + j] = (product[i + j] + left * right) % p
    return trim(product)


def divmod_poly(a, b, p):
    """Return the quotient and remainder of a/b over F_p."""
    a, b = trim(a, p), trim(b, p)
    if not b:
        raise ZeroDivisionError("division by the zero polynomial")
    if len(a) < len(b):
        return (), a
    remainder = list(a)
    quotient = [0] * (len(a) - len(b) + 1)
    inverse_leading = pow(b[-1], -1, p)
    while len(remainder) >= len(b):
        shift = len(remainder) - len(b)
        coefficient = remainder[-1] * inverse_leading % p
        quotient[shift] = coefficient
        for j, value in enumerate(b):
            remainder[shift + j] = (remainder[shift + j] - coefficient * value) % p
        while remainder and not remainder[-1]:
            remainder.pop()
    return trim(quotient), tuple(remainder)


def _monic(poly, p):
    if not poly:
        return ()
    inverse_leading = pow(poly[-1], -1, p)
    return tuple(value * inverse_leading % p for value in poly)


def gcd_poly(a, b, p):
    """Return the monic polynomial gcd; gcd(0,0) is zero."""
    a, b = trim(a, p), trim(b, p)
    while b:
        a, b = b, divmod_poly(a, b, p)[1]
    return _monic(a, p)


def powmod_poly(base, exponent, modulus, p):
    """Compute base**exponent modulo a nonzero polynomial using squaring."""
    _integer(exponent, "exponent")
    if exponent < 0:
        raise ValueError("exponent must be nonnegative")
    modulus = trim(modulus, p)
    if not modulus:
        raise ZeroDivisionError("modulus is the zero polynomial")
    base = divmod_poly(base, modulus, p)[1]
    result = divmod_poly((1,), modulus, p)[1]
    while exponent:
        if exponent & 1:
            result = divmod_poly(mul(result, base, p), modulus, p)[1]
        exponent >>= 1
        if exponent:
            base = divmod_poly(mul(base, base, p), modulus, p)[1]
    return result


def evaluate_poly(coeffs, x, p):
    """Evaluate an integer polynomial at an integer field representative."""
    _integer(x, "x")
    result = 0
    for coefficient in reversed(trim(coeffs, p)):
        result = (result * x + coefficient) % p
    return result


def field_roots(coeffs, p, seed=0, max_trials=128):
    """Isolate and certify all distinct roots of L in the prime field F_p.

    ``max_trials`` caps the total randomized splitting attempts across all
    factors.  ``status='complete'`` means exactly ``restricted_degree`` roots
    were verified, where that degree is the degree of gcd(L,x**p-x).

    A zero polynomial constrains no field elements: its status is
    ``'unrestricted'`` and ``domain_size`` is p.  The roots list is materialized
    only for p <= 100000; ``roots_enumerated`` distinguishes those cases.
    """
    _field(p)
    _integer(seed, "seed")
    _integer(max_trials, "max_trials")
    if max_trials < 0:
        raise ValueError("max_trials must be nonnegative")
    polynomial = trim(coeffs, p)
    degree = len(polynomial) - 1
    if not polynomial:
        enumerated = p <= _UNRESTRICTED_ENUMERATION_LIMIT
        return {
            "roots": list(range(p)) if enumerated else [],
            "status": "unrestricted",
            "degree": -1,
            "restricted_degree": p,
            "domain_size": p,
            "roots_enumerated": enumerated,
            "split_trials": 0,
            "reason": "L is zero over F_p, so every field element is allowed.",
        }
    if degree == 0:
        restricted = (1,)
    else:
        frobenius = powmod_poly((0, 1), p, polynomial, p)
        restricted = gcd_poly(polynomial, sub(frobenius, (0, 1), p), p)
    restricted_degree = len(restricted) - 1
    roots = set()
    trials = 0
    if p == 2:
        roots.update(x for x in (0, 1) if evaluate_poly(polynomial, x, p) == 0)
    else:
        # Removing x first makes root zero deterministic and reduces factor size.
        if restricted_degree and restricted[0] == 0:
            roots.add(0)
            restricted = restricted[1:]
        stack = [restricted] if len(restricted) > 1 else []
        rng = random.Random(seed)
        while stack:
            factor = stack.pop()
            factor_degree = len(factor) - 1
            if factor_degree == 1:
                roots.add(-factor[0] * pow(factor[1], -1, p) % p)
                continue
            while trials < max_trials:
                trials += 1
                candidate = tuple(rng.randrange(p) for _ in range(factor_degree))
                divisor = gcd_poly(candidate, factor, p)
                if not 1 < len(divisor) < len(factor):
                    character = powmod_poly(candidate, (p - 1) // 2, factor, p)
                    divisor = gcd_poly(sub(character, (1,), p), factor, p)
                if 1 < len(divisor) < len(factor):
                    quotient, remainder = divmod_poly(factor, divisor, p)
                    if remainder:
                        raise ArithmeticError("internal factor split did not divide exactly")
                    stack.extend((divisor, _monic(quotient, p)))
                    break
            # With no trials left, remaining nonlinear factors stay unresolved.
    verified = sorted(root for root in roots if evaluate_poly(polynomial, root, p) == 0)
    complete = len(verified) == restricted_degree
    return {
        "roots": verified,
        "status": "complete" if complete else "unresolved",
        "degree": degree,
        "restricted_degree": restricted_degree,
        "split_trials": trials,
        "reason": (
            f"All {restricted_degree} distinct roots of gcd(L, x^p-x) verified exactly."
            if complete else
            f"Splitting budget exhausted: {len(verified)} of {restricted_degree} roots verified."
        ),
    }
