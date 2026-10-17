"""Representation-aware domain construction; no Gröbner computations."""

from .polynomial import field_roots, evaluate_poly, gcd_poly, divmod_poly


def linearized_value(terms, constant, x, field):
    value = constant
    for power, coefficient in terms:
        value = field.add(value, field.mul(coefficient, field.pow(x, field.characteristic ** power)))
    return value


def _odd_linearized_roots(terms, constant, field, max_domain):
    """Exact prime-field Gaussian elimination for affine Frobenius maps."""
    p, n = field.characteristic, field.degree
    merged = {}
    for power, coefficient in terms:
        if type(power) is not int or power < 0:
            raise ValueError("Linearized powers must be nonnegative integers")
        field.validate(coefficient)
        j = power % n
        merged[j] = field.add(merged.get(j, 0), coefficient)
    terms = [(j, a) for j, a in sorted(merged.items()) if a]
    field.validate(constant)

    def digits(value):
        result = []
        for _ in range(n):
            result.append(value % p)
            value //= p
        return result

    rhs = digits(field.neg(constant))
    rows = [[0] * n + [rhs[r]] for r in range(n)]
    for c in range(n):
        column = digits(linearized_value(terms, 0, p ** c, field))
        for r in range(n):
            rows[r][c] = column[r]
    rank, pivots = 0, []
    for c in range(n):
        pivot = next((r for r in range(rank, n) if rows[r][c]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][c], p - 2, p)
        rows[rank] = [x * inverse % p for x in rows[rank]]
        for r in range(n):
            if r != rank and rows[r][c]:
                a = rows[r][c]
                rows[r] = [(x - a * y) % p for x, y in zip(rows[r], rows[rank])]
        pivots.append(c)
        rank += 1
    base = {"method": "prime_field_linear_algebra", "rank": rank,
            "nullity": n - rank, "exact_domain_count": True}
    if any(not any(row[:n]) and row[n] for row in rows):
        return {**base, "status": "complete", "roots": [], "domain_size": 0,
                "reason": "Affine prime-field system is inconsistent"}
    origin = sum(rows[r][n] * p ** c for r, c in enumerate(pivots))
    free = [c for c in range(n) if c not in pivots]
    basis = []
    for c in free:
        vector = p ** c
        for r, pivot_col in enumerate(pivots):
            vector += (-rows[r][c] % p) * p ** pivot_col
        basis.append(vector)
    count = p ** len(free)
    base.update(domain_size=count, affine_origin=origin, affine_basis=basis)
    if count > max_domain:
        return {**base, "status": "budget_exceeded", "roots": [],
                "reason": f"Affine domain has {count} elements; cap is {max_domain}"}
    roots = [origin]
    for vector in basis:
        roots = [field.add(x, field.mul(scalar, vector)) for x in roots for scalar in range(p)]
    roots.sort()
    if any(linearized_value(terms, constant, x, field) for x in roots):
        raise ArithmeticError("Prime-field affine elimination verification failed")
    return {**base, "status": "complete", "roots": roots,
            "reason": "Exact affine solve and evaluation certify the entire domain"}


def linearized_roots(terms, constant, field, max_domain=4096):
    """Solve sum a_j x^(2^j) + b = 0 by exact binary elimination.

    Frobenius exponents reduce modulo the extension degree as maps on K.
    Returns an exact count even when enumerating the affine domain is refused.
    """
    if field.characteristic != 2:
        return _odd_linearized_roots(terms, constant, field, max_domain)
    n = field.degree
    merged = {}
    for power, coefficient in terms:
        if type(power) is not int or power < 0:
            raise ValueError("Linearized powers must be nonnegative integers")
        field.validate(coefficient)
        merged[power % n] = merged.get(power % n, 0) ^ coefficient
    terms = [(j, a) for j, a in sorted(merged.items()) if a]
    field.validate(constant)
    rows = [((constant >> r) & 1) << n for r in range(n)]
    for c in range(n):
        image = linearized_value(terms, 0, 1 << c, field)
        for r in range(n):
            rows[r] |= ((image >> r) & 1) << c
    pivot_cols = []
    rank = 0
    for col in range(n):
        pivot = next((r for r in range(rank, n) if (rows[r] >> col) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for r in range(n):
            if r != rank and ((rows[r] >> col) & 1):
                rows[r] ^= rows[rank]
        pivot_cols.append(col)
        rank += 1
    mask = (1 << n) - 1
    inconsistent = any((row & mask) == 0 and ((row >> n) & 1) for row in rows)
    if inconsistent:
        return {"status": "complete", "roots": [], "domain_size": 0,
                "method": "binary_linear_algebra", "rank": rank,
                "nullity": n - rank, "reason": "Affine binary system is inconsistent"}
    particular = sum(((rows[r] >> n) & 1) << col for r, col in enumerate(pivot_cols))
    free = [c for c in range(n) if c not in pivot_cols]
    basis = []
    for c in free:
        vector = 1 << c
        for r, col in enumerate(pivot_cols):
            vector |= ((rows[r] >> c) & 1) << col
        basis.append(vector)
    count = 1 << len(free)
    result = {"method": "binary_linear_algebra", "rank": rank,
              "nullity": len(free), "domain_size": count,
              "affine_origin": particular, "affine_basis": basis,
              "exact_domain_count": True}
    if count > max_domain:
        result.update(status="budget_exceeded", roots=[],
                      reason=f"Affine domain has {count} elements; cap is {max_domain}")
        return result
    roots = [particular]
    for vector in basis:
        roots += [x ^ vector for x in roots]
    roots.sort()
    if any(linearized_value(terms, constant, x, field) for x in roots):
        raise ArithmeticError("Binary elimination verification failed")
    result.update(status="complete", roots=roots,
                  reason="Exact affine solve and evaluation certify the entire domain")
    return result


def polynomial_domain(numerator, denominator, field, max_domain=4096,
                      max_degree=256, seed=0, max_trials=128):
    if len(numerator) - 1 > max_degree or len(denominator) - 1 > max_degree:
        return {"status": "budget_exceeded", "roots": [],
                "reason": f"Expanded polynomial exceeds degree cap {max_degree}"}
    if not any(denominator):
        raise ValueError("A rational L must have a nonzero denominator polynomial")
    # L is an element of K(X): cancel common factors before identifying poles.
    common = gcd_poly(numerator, denominator, field)
    numerator, rem_a = divmod_poly(numerator, common, field)
    denominator, rem_b = divmod_poly(denominator, common, field)
    if any(rem_a) or any(rem_b):
        raise ArithmeticError("Rational-function cancellation failed")
    result = field_roots(numerator, field, seed=seed, max_trials=max_trials)
    result["method"] = "frobenius_gcd_trace_splitting" if field.characteristic == 2 else "frobenius_gcd_quadratic_character_splitting"
    result["rational_semantics"] = "coprime numerator/denominator; poles excluded after cancellation"
    if result["status"] == "unrestricted":
        if field.order > max_domain:
            return {**result, "status": "budget_exceeded", "roots": [],
                    "reason": "Zero numerator leaves an unmaterialized field-sized domain"}
        result["roots"] = list(range(field.order))
        result["status"] = "complete"
        result["roots_enumerated"] = True
    if result["status"] != "complete":
        return result
    roots = [x for x in result["roots"] if evaluate_poly(denominator, x, field) != 0]
    result["excluded_poles"] = len(result["roots"]) - len(roots)
    result["domain_size"] = len(roots)
    result["exact_domain_count"] = True
    if len(roots) > max_domain:
        result.update(status="budget_exceeded", roots=[],
                      reason=f"Root domain has {len(roots)} elements; cap is {max_domain}")
    else:
        result["roots"] = roots
    return result
