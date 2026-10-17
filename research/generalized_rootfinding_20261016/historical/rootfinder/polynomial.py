"""Field-protocol polynomial arithmetic and exact finite-field root isolation.

Coefficients are canonical field elements in ascending order, and () is zero.
Binary and small-prime fields dispatch to the existing optimized arithmetic.
Other finite fields use exact arithmetic and Frobenius restriction followed by
randomized, exactly verified splitting.  A budget exhaustion is UNRESOLVED;
root count and root membership are never inferred from a heuristic.
"""

import random

from . import binaryfield, finitefield
from .binaryfield import BinaryField, _integer
from .fields import PrimeField, ExtensionField


_UNRESTRICTED_ENUMERATION_LIMIT = 100_000


def _field(field):
    # Explicit protocol rather than only the built-in classes permits future
    # optimized backends to use the same exact polynomial layer.
    for attribute in ("characteristic", "degree", "order", "validate", "add",
                      "sub", "neg", "mul", "pow", "inv"):
        if not hasattr(field, attribute):
            raise TypeError(f"field does not implement finite-field protocol: {attribute}")
    return field


def _trim(values):
    while values and not values[-1]:
        values.pop()
    return tuple(values)


def trim(coeffs, field=None):
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
            raise ValueError("polynomial coefficients must be canonical nonnegative field elements")
        if field is not None:
            field.validate(value)
    return _trim(values)


def _small_prime(field):
    return isinstance(field, PrimeField) and field.order < (1 << 64)


def _operations(field):
    return (getattr(field, "_raw_add", field.add),
            getattr(field, "_raw_sub", field.sub),
            getattr(field, "_raw_mul", field.mul))


def _add(a, b, field, subtract=False):
    operation = _operations(field)[1 if subtract else 0]
    return _trim([operation(a[i] if i < len(a) else 0,
                            b[i] if i < len(b) else 0)
                  for i in range(max(len(a), len(b)))])


def add(a, b, field):
    _field(field)
    if isinstance(field, BinaryField):
        return binaryfield.add(a, b, field)
    a, b = trim(a, field), trim(b, field)
    if _small_prime(field):
        return finitefield.add(a, b, field.order)
    return _add(a, b, field)


def sub(a, b, field):
    _field(field)
    if isinstance(field, BinaryField):
        return binaryfield.sub(a, b, field)
    a, b = trim(a, field), trim(b, field)
    if _small_prime(field):
        return finitefield.sub(a, b, field.order)
    return _add(a, b, field, subtract=True)


def _mul(a, b, field):
    if not a or not b:
        return ()
    plus, _, multiply = _operations(field)
    result = [0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        if left:
            for j, right in enumerate(b):
                if right:
                    result[i + j] = plus(result[i + j], multiply(left, right))
    return _trim(result)


def mul(a, b, field):
    _field(field)
    if isinstance(field, BinaryField):
        return binaryfield.mul(a, b, field)
    a, b = trim(a, field), trim(b, field)
    if _small_prime(field):
        return finitefield.mul(a, b, field.order)
    return _mul(a, b, field)


def _divmod(a, b, field):
    if not b:
        raise ZeroDivisionError("division by the zero polynomial")
    if len(a) < len(b):
        return (), a
    _, minus, multiply = _operations(field)
    remainder = list(a)
    quotient = [0] * (len(a) - len(b) + 1)
    leading_inverse = field.inv(b[-1]) if b[-1] != 1 else 1
    while len(remainder) >= len(b):
        shift = len(remainder) - len(b)
        coefficient = multiply(remainder[-1], leading_inverse)
        quotient[shift] = coefficient
        for j, value in enumerate(b):
            if value:
                remainder[shift + j] = minus(remainder[shift + j], multiply(coefficient, value))
        while remainder and not remainder[-1]:
            remainder.pop()
    return _trim(quotient), tuple(remainder)


def divmod_poly(a, b, field):
    _field(field)
    if isinstance(field, BinaryField):
        return binaryfield.divmod_poly(a, b, field)
    a, b = trim(a, field), trim(b, field)
    if _small_prime(field):
        return finitefield.divmod_poly(a, b, field.order)
    return _divmod(a, b, field)


def _monic(poly, field):
    if not poly or poly[-1] == 1:
        return poly
    multiply = _operations(field)[2]
    inverse = field.inv(poly[-1])
    return tuple(multiply(coefficient, inverse) for coefficient in poly)


def _gcd(a, b, field):
    while b:
        a, b = b, _divmod(a, b, field)[1]
    return _monic(a, field)


def gcd_poly(a, b, field):
    _field(field)
    if isinstance(field, BinaryField):
        return binaryfield.gcd_poly(a, b, field)
    a, b = trim(a, field), trim(b, field)
    if _small_prime(field):
        return finitefield.gcd_poly(a, b, field.order)
    return _gcd(a, b, field)


def _powmod(base, exponent, modulus, field):
    if not modulus:
        raise ZeroDivisionError("modulus is the zero polynomial")
    base = _divmod(base, modulus, field)[1]
    result = _divmod((1,), modulus, field)[1]
    while exponent:
        if exponent & 1:
            result = _divmod(_mul(result, base, field), modulus, field)[1]
        exponent >>= 1
        if exponent:
            base = _divmod(_mul(base, base, field), modulus, field)[1]
    return result


def powmod_poly(base, exponent, modulus, field):
    _field(field)
    _integer(exponent, "exponent")
    if exponent < 0:
        raise ValueError("exponent must be nonnegative")
    if isinstance(field, BinaryField):
        return binaryfield.powmod_poly(base, exponent, modulus, field)
    base, modulus = trim(base, field), trim(modulus, field)
    if _small_prime(field):
        return finitefield.powmod_poly(base, exponent, modulus, field.order)
    return _powmod(base, exponent, modulus, field)


def evaluate_poly(coeffs, x, field):
    _field(field)
    if isinstance(field, BinaryField):
        return binaryfield.evaluate_poly(coeffs, x, field)
    x = field.validate(x)
    coefficients = trim(coeffs, field)
    if _small_prime(field):
        return finitefield.evaluate_poly(coefficients, x, field.order)
    plus, _, multiply = _operations(field)
    result = 0
    for value in reversed(coefficients):
        result = plus(multiply(result, x), value)
    return result


def field_roots(coeffs, field, seed=0, max_trials=128):
    """Return verified distinct roots and an exact Frobenius-gcd root count.

    ``complete`` means membership and count were checked exactly.  The zero
    polynomial is ``unrestricted`` and large fields are never materialized.
    A global randomized split budget returns ``unresolved`` with individually
    verified roots and the exact distinct-root count when exhausted.
    """
    _field(field)
    _integer(seed, "seed")
    _integer(max_trials, "max_trials")
    if max_trials < 0:
        raise ValueError("max_trials must be nonnegative")
    polynomial = trim(coeffs, field)
    if isinstance(field, BinaryField):
        return binaryfield.binary_field_roots(polynomial, field, seed, max_trials)
    if _small_prime(field):
        result = finitefield.field_roots(polynomial, field.order, seed, max_trials)
        result.setdefault("domain_size", field.order)
        result.setdefault("roots_enumerated", result["status"] == "complete")
        return result
    degree = len(polynomial) - 1
    if not polynomial:
        enumerated = field.order <= _UNRESTRICTED_ENUMERATION_LIMIT
        return {"roots": list(range(field.order)) if enumerated else [],
                "status": "unrestricted", "degree": -1,
                "restricted_degree": field.order, "domain_size": field.order,
                "roots_enumerated": enumerated, "split_trials": 0,
                "reason": "L is zero, so every field element is allowed."}
    if degree == 0:
        restricted = (1,)
    else:
        frobenius = _powmod((0, 1), field.order, polynomial, field)
        restricted = _gcd(polynomial, _add(frobenius, (0, 1), field, True), field)
    count = len(restricted) - 1
    roots = set()
    if count and restricted[0] == 0:
        roots.add(0)
        restricted = restricted[1:]
    stack = [restricted] if len(restricted) > 1 else []
    rng, trials = random.Random(seed), 0
    while stack:
        factor = stack.pop()
        if len(factor) == 2:
            roots.add(field.neg(field.mul(factor[0], field.inv(factor[1]))))
            continue
        while trials < max_trials:
            trials += 1
            if field.characteristic == 2:
                scalar = rng.randrange(field.order)
                power = _divmod((0, scalar), factor, field)[1]
                trace = ()
                for step in range(field.degree):
                    trace = _add(trace, power, field)
                    if step + 1 < field.degree:
                        power = _powmod(power, 2, factor, field)
                divisor = _gcd(trace, factor, field)
            else:
                candidate = _trim([rng.randrange(field.order) for _ in range(len(factor) - 1)])
                divisor = _gcd(candidate, factor, field)
                if not 1 < len(divisor) < len(factor):
                    character = _powmod(candidate, (field.order - 1) // 2, factor, field)
                    divisor = _gcd(_add(character, (1,), field, True), factor, field)
            if 1 < len(divisor) < len(factor):
                quotient, remainder = _divmod(factor, divisor, field)
                if remainder:
                    raise ArithmeticError("internal root split did not divide exactly")
                stack.extend((divisor, _monic(quotient, field)))
                break
    verified = sorted(root for root in roots if evaluate_poly(polynomial, root, field) == 0)
    complete = len(verified) == count
    return {"roots": verified, "status": "complete" if complete else "unresolved",
            "degree": degree, "restricted_degree": count, "domain_size": field.order,
            "roots_enumerated": complete, "split_trials": trials,
            "reason": (f"All {count} distinct roots of gcd(L, X^q-X) verified exactly."
                       if complete else
                       f"Splitting budget exhausted: {len(verified)} of {count} roots verified.")}
