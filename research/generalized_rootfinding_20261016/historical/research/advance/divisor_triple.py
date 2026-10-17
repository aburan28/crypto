"""Exact conditional divisor-norm search for rational binary-curve triples.

This experimental backend returns all ordered triples of *finite rational
points* with abscissae in an affine binary subspace, summing to a finite
target.  It does not claim completeness for unrestricted x-only S4 roots,
novelty, or a measured speedup.  No Groebner basis is used.

For E: y²+xy=x³+A*x²+B (B != 0), put T=-target.  A degree-four
function F=x²+c*x+d+b*y, d=z²+c*z+b*yT, vanishes at T.  Its norm,
divided by X+z, is a monic cubic H.  Its coefficients are conditional
binary-affine maps in c.  Necessary product-span constraints yield an
exact affine family; compact domain divisibility and point/group replay
filter that family.  Repeated roots and b=0 are explicitly supported.

Every rational source triple determines a principal divisor triple+T-4O
and hence such a normalized F in L(4O).  Thus the enumeration is complete
when all enumeration and memory caps permit completion.  Any cap
exhaustion returns unresolved with only a verified solution lower bound.
"""

from itertools import permutations, product

from rootfinder.binaryfield import BinaryField
from rootfinder import polynomial as poly


def _integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _basis(columns):
    original_basis, pivots = [], {}
    for original in columns:
        value = _residual(original, pivots)
        if value:
            pivots[value.bit_length() - 1] = value
            original_basis.append(original)
    return tuple(original_basis), pivots


def _residual(value, pivots):
    # Reducing all pivots, rather than stopping at the first free bit,
    # makes this a fixed linear quotient map, suitable for affine systems.
    for pivot in sorted(pivots, reverse=True):
        if (value >> pivot) & 1:
            value ^= pivots[pivot]
    return value


def _solve(columns, rhs):
    pivots, kernel = {}, []
    for index, original in enumerate(columns):
        value, coordinate = original, 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot][0]
                coordinate ^= pivots[pivot][1]
            else:
                pivots[pivot] = (value, coordinate)
                break
        if not value:
            kernel.append(coordinate)
    particular, value = 0, rhs
    while value:
        pivot = value.bit_length() - 1
        if pivot not in pivots:
            return {"particular": None, "basis": tuple(kernel),
                    "rank": len(pivots), "nullity": len(kernel),
                    "cardinality": 0, "consistent": False}
        value ^= pivots[pivot][0]
        particular ^= pivots[pivot][1]
    return {"particular": particular, "basis": tuple(kernel),
            "rank": len(pivots), "nullity": len(kernel),
            "cardinality": 1 << len(kernel), "consistent": True}


def enumerate_family(family):
    """Enumerate an exact coordinate-bitmask family; callers enforce caps."""
    if not family["consistent"]:
        return
    value = family["particular"]
    yield value
    for counter in range(1, family["cardinality"]):
        bit = (counter & -counter).bit_length() - 1
        value ^= family["basis"][bit]
        yield value


def on_curve(field, curve_a, curve_b, point):
    if point is None:
        return True
    x, y = point
    return (field.square(y) ^ field.mul(x, y) ^ field.mul(field.square(x), x)
            ^ field.mul(curve_a, field.square(x)) ^ curve_b) == 0


def negate_point(point):
    return None if point is None else (point[0], point[0] ^ point[1])


def add_points(field, curve_a, left, right):
    """Ordinary binary Weierstrass group law, with None denoting O."""
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u:
        if y != v or x == 0:
            return None
        slope = x ^ field.mul(y, field.inv(x))
        w = field.square(slope) ^ slope ^ curve_a
        return w, field.square(x) ^ field.mul(slope ^ 1, w)
    slope = field.mul(y ^ v, field.inv(x ^ u))
    w = field.square(slope) ^ slope ^ x ^ u ^ curve_a
    return w, field.mul(slope, x ^ w) ^ w ^ y


def prepare_model(field, curve_a, curve_b, target, offset, basis):
    """Prepare verified support spans and the compact outer parameter family.

    This never enumerates the domain, the field, or the outer family.  The
    returned model permits honest full-width rank profiling via
    conditional_family(model,b), independently of solve_triples.
    """
    if not isinstance(field, BinaryField):
        raise TypeError("field must be a verified polynomial-basis BinaryField")
    verified = BinaryField(field.degree, field.modulus)
    if field.characteristic != 2 or field.order != verified.order:
        raise ValueError("field metadata is inconsistent with its modulus")
    curve_a, curve_b, offset = map(field.validate, (curve_a, curve_b, offset))
    if curve_b == 0:
        raise ValueError("curve_b must be nonzero for a nonsingular ordinary curve")
    if target is None or isinstance(target, (str, bytes, dict)):
        raise ValueError("target must be a finite rational point")
    try:
        target = tuple(target)
    except TypeError as error:
        raise TypeError("target must be a coordinate pair") from error
    if len(target) != 2:
        raise ValueError("target must have exactly two coordinates")
    target = tuple(field.validate(value) for value in target)
    if not on_curve(field, curve_a, curve_b, target):
        raise ValueError("target is not on the curve")
    if isinstance(basis, (str, bytes, dict)):
        raise TypeError("basis must be an iterable of field elements")
    basis = tuple(field.validate(value) for value in basis)
    independent, vp = _basis(basis)
    if len(independent) != len(basis):
        raise ValueError("domain basis must be independent over F2")
    w2, w2p = _basis(field.mul(u, v) for i, u in enumerate(basis)
                     for v in basis[i:])
    w3, w3p = _basis(field.mul(field.mul(u, v), w)
                     for i, u in enumerate(basis)
                     for j in range(i, len(basis)) for v in (basis[j],)
                     for w in basis[j:])
    units = tuple(1 << index for index in range(field.degree))
    outer = _solve(tuple(_residual(field.square(unit) ^ unit, vp)
                         for unit in units), _residual(target[0] ^ offset, vp))
    # Compact linearized L_V: coefficient i belongs to X^(2**i).
    linearized = (1,)
    prefixes, scalars = [linearized], []
    for vector in basis:
        scalar, power = 0, vector
        for coefficient in linearized:
            scalar ^= field.mul(coefficient, power)
            power = field.square(power)
        scalars.append(scalar)
        shifted = (0,) + tuple(field.square(value) for value in linearized)
        linearized = tuple(shifted[i] ^ (field.mul(scalar, linearized[i])
                                        if i < len(linearized) else 0)
                           for i in range(len(shifted)))
        prefixes.append(linearized)
    domain_constant, power = 0, offset
    for coefficient in linearized:
        domain_constant ^= field.mul(coefficient, power)
        power = field.square(power)
    z = target[0]
    unit_squares = tuple(field.square(unit) for unit in units)
    conditional_base_columns = tuple(
        _residual(squared, w2p) | (
            _residual(field.mul(z, squared) ^ field.mul(offset, squared), w3p)
            << field.degree)
        for squared in unit_squares)
    return {"field": field, "curve_a": curve_a, "curve_b": curve_b,
            "target": target, "negative_target": negate_point(target),
            "offset": offset, "basis": basis, "domain_pivots": vp,
            "w2_basis": w2, "w3_basis": w3, "w2_pivots": w2p,
            "w3_pivots": w3p, "units": units,
            "unit_squares": unit_squares,
            "conditional_base_columns": conditional_base_columns,
            "outer_family": outer, "linearized": linearized,
            "linearized_prefixes": tuple(prefixes),
            "linearized_scalars": tuple(scalars),
            "domain_constant": domain_constant}


def conditional_family(model, b):
    """Exact affine c family at b, with measured rank and nullity.

    This family enforces necessary product-span constraints only.  Its
    cardinality is the exact number of relaxed cubics, not the number of
    valid point triples.
    """
    f, a = model["field"], model["offset"]
    b = f.validate(b)
    z, yt = model["negative_target"]
    e1 = f.square(b) ^ b ^ z
    if _residual(e1 ^ a, model["domain_pivots"]):
        return {"particular": None, "basis": (), "rank": 0,
                "nullity": 0, "cardinality": 0, "consistent": False,
                "outer_admitted": False}
    z2, b2, a2 = f.square(z), f.square(b), f.square(a)
    e20 = f.mul(b2, model["curve_a"]) ^ f.mul(z, e1)
    e30 = f.mul(b2, yt ^ z2 ^ f.mul(z, model["curve_a"])) ^ f.mul(z2, z)
    rhs2 = _residual(e20 ^ a2, model["w2_pivots"])
    rhs3 = _residual(e30 ^ f.mul(a2, e1) ^ f.mul(a, e20)
                     ^ f.mul(a2, a), model["w3_pivots"])
    # Multiplication by polynomial-basis alpha is one exact shift/reduce.
    # Share all b-independent Frobenius contributions across outer values,
    # and generate both multiplication columns in O(n) such shifts.
    columns, product_b, product_ab = [], b, f.mul(a, b)
    for base_column in model["conditional_base_columns"]:
        low = _residual(product_b, model["w2_pivots"])
        high = _residual(product_ab, model["w3_pivots"])
        columns.append(base_column ^ low ^ (high << f.degree))
        product_b <<= 1
        product_ab <<= 1
        if product_b & f.order:
            product_b ^= f.modulus
        if product_ab & f.order:
            product_ab ^= f.modulus
    family = _solve(columns, rhs2 | (rhs3 << f.degree))
    family["outer_admitted"] = True
    return family


def cubic_from_parameters(model, b, c):
    """Construct H by literal norm and exact division, without formulas."""
    f = model["field"]
    b, c = f.validate(b), f.validate(c)
    z, yt = model["negative_target"]
    d = f.square(z) ^ f.mul(c, z) ^ f.mul(b, yt)
    ap = (d, c, 1)
    norm = poly.add(poly.mul(ap, ap, f), poly.mul((0, b), ap, f), f)
    b2 = f.square(b)
    curve = (model["curve_b"], 0, model["curve_a"], 1)
    norm = poly.add(norm, tuple(f.mul(b2, value) for value in curve), f)
    quotient, remainder = poly.divmod_poly(norm, (z, 1), f)
    if remainder or len(quotient) != 4 or quotient[-1] != 1:
        raise ArithmeticError("curve-function norm did not yield a monic cubic")
    e1 = b2 ^ b ^ z
    e2 = f.square(c) ^ f.mul(b, c) ^ f.mul(b2, model["curve_a"]) ^ f.mul(z, e1)
    e3 = f.mul(z, f.square(c)) ^ f.mul(b2, yt ^ f.square(z)
             ^ f.mul(z, model["curve_a"])) ^ f.mul(f.square(z), z)
    if quotient != (e3, e2, e1, 1):
        raise ArithmeticError("conditional coefficient formulas disagree with norm")
    return quotient, d


def domain_contains(model, x):
    return _residual(x ^ model["offset"], model["domain_pivots"]) == 0


def _domain_remainder(model, cubic):
    """Reduce the compact linearized vanishing chain in constant degree."""
    f = model["field"]
    remainder = poly.divmod_poly((model["offset"], 1), cubic, f)[1]
    for scalar in model["linearized_scalars"]:
        remainder = poly.add(poly.powmod_poly(remainder, 2, cubic, f),
                             tuple(f.mul(scalar, value) for value in remainder), f)
    return remainder


def cubic_in_domain(model, cubic):
    """Require H | L_D^3, including every repeated root multiplicity.

    L_D is squarefree and has precisely the affine domain as its roots.
    Since deg H=3, divisibility by its cube is equivalent to requiring all
    cubic roots, over an algebraic closure and with multiplicities, in D.
    Requiring H | L_D would incorrectly discard repeated source points.
    """
    return not poly.powmod_poly(_domain_remainder(model, cubic), 3, cubic,
                               model["field"])


def _linearized_value(coefficients, value, field):
    result = 0
    for coefficient in coefficients:
        result ^= field.mul(coefficient, value)
        value = field.square(value)
    return result


def cubic_domain_roots(model, cubic):
    """Certify support and split via nested affine subspaces, deterministically.

    Returns None if some cubic root lies outside D.  Otherwise returns all
    distinct roots and a gcd-step count.  The support certificate gives
    radical(H)=gcd(H,L_D), so neither an 83-step Frobenius restriction nor
    randomized trace splitting is needed.  Each nested L_i cuts the known
    affine domain into its two cosets; a single compact chain is reduced
    modulo the radical, whose degree is at most three.
    """
    f = model["field"]
    remainder = _domain_remainder(model, cubic)
    if poly.powmod_poly(remainder, 3, cubic, f):
        return None
    radical = poly.gcd_poly(cubic, remainder, f)
    remainders = [poly.divmod_poly((model["offset"], 1), radical, f)[1]]
    for scalar in model["linearized_scalars"]:
        last = remainders[-1]
        remainders.append(poly.add(poly.powmod_poly(last, 2, radical, f),
                                   tuple(f.mul(scalar, value) for value in last), f))
    roots, steps = [], 0
    stack = [(radical, len(model["basis"]), 0)]
    while stack:
        factor, level, shift = stack.pop()
        if len(factor) == 2:
            root = f.mul(factor[0], f.inv(factor[1]))
            if not domain_contains(model, root) or poly.evaluate_poly(cubic, root, f):
                raise ArithmeticError("affine subspace splitting produced an invalid root")
            roots.append(root)
            continue
        if len(factor) <= 1:
            continue
        if level == 0:
            raise ArithmeticError("nonlinear factor survived singleton domain splitting")
        level -= 1
        constant = _linearized_value(model["linearized_prefixes"][level], shift, f)
        cut = poly.add(remainders[level], (constant,), f)
        left = poly.gcd_poly(factor, cut, f)
        right, division_remainder = poly.divmod_poly(factor, left, f)
        if division_remainder:
            raise ArithmeticError("affine subspace split did not divide exactly")
        steps += 1
        if len(left) > 1:
            stack.append((left, level, shift))
        if len(right) > 1:
            stack.append((right, level, shift ^ model["basis"][level]))
    if len(roots) != len(radical) - 1 or len(set(roots)) != len(roots):
        raise ArithmeticError("affine subspace factor count disagrees with radical degree")
    return {"roots": sorted(roots), "distinct_root_count_exact": len(roots),
            "gcd_split_steps": steps}


def _ordinate_roots(model, x):
    """Exact rational lifts without randomized polynomial root splitting."""
    f = model["field"]
    rhs = f.mul(f.square(x), x) ^ f.mul(model["curve_a"], f.square(x)) ^ model["curve_b"]
    if x == 0:
        return (f.pow(rhs, 1 << (f.degree - 1)),)
    if f.degree & 1:
        inverse = f.inv(x)
        value = f.mul(rhs, f.square(inverse))
        half_trace, power = 0, value
        for _ in range((f.degree + 1) // 2):
            half_trace ^= power
            power = f.square(f.square(power))
        if f.square(half_trace) ^ half_trace != value:
            return ()
        y = f.mul(x, half_trace)
        return tuple(sorted((y, y ^ x)))
    family = _solve(tuple(squared ^ f.mul(x, unit)
                           for unit, squared in zip(model["units"], model["unit_squares"])), rhs)
    if family["cardinality"] > 2:
        raise ArithmeticError("binary ordinate equation has more than two roots")
    return tuple(sorted(enumerate_family(family)))


def solve_triples(field, curve_a, curve_b, target, offset, basis, *,
                  max_outer=1 << 19, max_candidates=1 << 23,
                  max_output=10_000, max_unique=1_000_000,
                  max_split_trials=2048):
    """Return exact ordered rational point triples, or an honest partial set.

    Output clipping alone retains an exact count.  The separate max_unique
    memory cap bounds the global deduplication set.  Each outer/c family is
    checked against its enumeration cap *before* that family is visited.
    Compact subspace splitting and exact ordinate linear algebra are
    deterministic.  max_split_trials remains a compatibility parameter;
    this backend spends zero randomized split trials.  No uniform-nullity,
    domain-generic, or performance assumption appears in completeness.
    """
    caps = {name: _integer(value, name) for name, value in (
        ("max_outer", max_outer), ("max_candidates", max_candidates),
        ("max_output", max_output), ("max_unique", max_unique),
        ("max_split_trials", max_split_trials))}
    model = prepare_model(field, curve_a, curve_b, target, offset, basis)
    stats = {"domain_dimension": len(model["basis"]),
             "domain_cardinality": 1 << len(model["basis"]),
             "w2_dimension": len(model["w2_basis"]),
             "w3_dimension": len(model["w3_basis"]),
             "outer_rank": model["outer_family"]["rank"],
             "outer_nullity": model["outer_family"]["nullity"],
             "outer_cardinality_exact": model["outer_family"]["cardinality"],
             "outer_visited": 0, "inner_rank_histogram": {},
             "inner_nullity_histogram": {}, "inner_inconsistent": 0,
             "relaxed_candidates_planned": 0, "candidates_visited": 0,
             "refused_family_cardinality": None,
             "candidate_bound_needed_at_stop": None,
             "cubics_supported": 0, "cubics_factored": 0,
             "split_trials": 0, "ordered_point_replays": 0,
             "deterministic_source_gcd_steps": 0, "ordinate_lifts": 0,
             "unique_solutions_verified": 0}
    seen, solutions = set(), []

    def result(complete, reason=None):
        return {"status": "complete" if complete else "unresolved",
                "complete": complete, "contract": "ordered_finite_rational_point_triples",
                "solution_count_exact": len(seen) if complete else None,
                "solution_count_lower_bound": len(seen), "solutions": solutions,
                "output_complete": complete and len(solutions) == len(seen),
                "reason": reason, "stats": stats, "caps": caps}

    if stats["outer_cardinality_exact"] > max_outer:
        return result(False, "outer parameter family exceeds max_outer before enumeration")
    for b in enumerate_family(model["outer_family"]):
        stats["outer_visited"] += 1
        family = conditional_family(model, b)
        for key, value in (("inner_rank_histogram", family["rank"]),
                           ("inner_nullity_histogram", family["nullity"])):
            stats[key][value] = stats[key].get(value, 0) + 1
        if not family["consistent"]:
            stats["inner_inconsistent"] += 1
            continue
        if stats["relaxed_candidates_planned"] + family["cardinality"] > max_candidates:
            stats["refused_family_cardinality"] = family["cardinality"]
            stats["candidate_bound_needed_at_stop"] = (
                stats["relaxed_candidates_planned"] + family["cardinality"])
            return result(False, "conditional candidate family exceeds max_candidates before enumeration")
        stats["relaxed_candidates_planned"] += family["cardinality"]
        for c in enumerate_family(family):
            stats["candidates_visited"] += 1
            cubic, d = cubic_from_parameters(model, b, c)
            isolated = cubic_domain_roots(model, cubic)
            if isolated is None:
                continue
            stats["cubics_supported"] += 1
            stats["cubics_factored"] += 1
            stats["deterministic_source_gcd_steps"] += isolated["gcd_split_steps"]
            remaining, multiset = cubic, []
            for x in isolated["roots"]:
                if not domain_contains(model, x):
                    raise ArithmeticError("compact domain certificate disagrees with membership")
                while poly.evaluate_poly(remaining, x, field) == 0:
                    quotient, remainder = poly.divmod_poly(remaining, (x, 1), field)
                    if remainder:
                        raise ArithmeticError("root multiplicity division failed")
                    multiset.append(x)
                    remaining = quotient
            if remaining != (1,) or len(multiset) != 3:
                raise ArithmeticError("supported cubic did not split completely with multiplicities")
            ordinates = {}
            for x in set(multiset):
                if b:
                    y = field.mul(field.square(x) ^ field.mul(c, x) ^ d, field.inv(b))
                    ordinates[x] = (y,) if on_curve(field, curve_a, curve_b, (x, y)) else ()
                else:
                    ordinates[x] = _ordinate_roots(model, x)
                    stats["ordinate_lifts"] += 1
            for ordered_x in sorted(set(permutations(multiset))):
                for ordered_y in product(*(ordinates[x] for x in ordered_x)):
                    triple = tuple(zip(ordered_x, ordered_y))
                    stats["ordered_point_replays"] += 1
                    total = add_points(field, curve_a, add_points(field, curve_a,
                                      triple[0], triple[1]), triple[2])
                    if total != model["target"] or triple in seen:
                        continue
                    if len(seen) >= max_unique:
                        return result(False, "verified solution deduplication exceeds max_unique")
                    seen.add(triple)
                    stats["unique_solutions_verified"] += 1
                    if len(solutions) < max_output:
                        solutions.append(triple)
    return result(True)
