"""Independent bounded reconstruction of generic-worker factor-base recipes.

This checker uses scalar Python field/group arithmetic, row-reduced binary
linear maps and Berlekamp factorization. It never calls the measured builder.
The production factor scan's degree-24 cutoff and eight-representative sampler
batches are explicit implementation rules, not mathematical assumptions.
"""
import copy
from functools import lru_cache
import json

from generic_query_law import StdRng08, uint
from identity import factor_base_inventory, sha256
from oracle import Curve, irreducible_binary_polynomial, polynomial_remainder, require

MAX_ABSCISSAE = 1 << 20


def poly_product(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


def poly_gcd(a, b):
    while b:
        a, b = b, polynomial_remainder(a, b)
    return a


def poly_quotient(a, b):
    result = 0
    while a.bit_length() >= b.bit_length():
        shift = a.bit_length() - b.bit_length()
        result ^= 1 << shift
        a ^= b << shift
    require(a == 0, 'inexact binary polynomial division')
    return result


def nullspace(columns, dimension):
    """RREF of the transposed coefficient array, with ascending free columns."""
    rows = [sum(((column >> i) & 1) << j for j, column in enumerate(columns))
            for i in range(dimension)]
    pivots, pivot_row = [], 0
    for column in range(len(columns)):
        selected = next((i for i in range(pivot_row, len(rows)) if rows[i] >> column & 1), None)
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        for i in range(len(rows)):
            if i != pivot_row and rows[i] >> column & 1:
                rows[i] ^= rows[pivot_row]
        pivots.append(column)
        pivot_row += 1
    result = []
    for free in range(len(columns)):
        if free not in pivots:
            vector = 1 << free
            for row, pivot in zip(rows, pivots):
                if row >> free & 1:
                    vector |= 1 << pivot
            result.append(vector)
    return result


@lru_cache(maxsize=32)
def complete_factors(n):
    """Factor square-free x^n+1 without scanning 2^degree monic polynomials."""
    require(5 <= n <= 31 and n % 2 == 1, 'base factorization requires bounded odd degree')
    polynomial = (1 << n) | 1
    columns = [polynomial_remainder(1 << (2 * i), polynomial) ^ (1 << i) for i in range(n)]
    basis = nullspace(columns, n)
    factors = [polynomial]
    for vector in basis:
        split = []
        for factor in factors:
            divisor = poly_gcd(factor, vector)
            if 1 < divisor < factor:
                split.extend((divisor, poly_quotient(factor, divisor)))
            else:
                split.append(factor)
        factors = split
    factors.sort(key=lambda value: (value.bit_length(), value))
    product = 1
    for factor in factors:
        require(irreducible_binary_polynomial(factor), 'incomplete independent factorization')
        product = poly_product(product, factor)
    require(product == polynomial and len(set(factors)) == len(factors), 'factorization does not close')
    return tuple(factors)


def field_power(c, value, iterations):
    for _ in range(iterations):
        value = c.fm(value, value)
    return value


def lifts(c, x):
    """Odd-degree half-trace root first, then its negative; x=0 gives torsion."""
    if x == 0:
        return [(0, 1)]
    inv = c.inv(x)
    rhs = x ^ c.a ^ c.fm(inv, inv)
    root, power = rhs, rhs
    for _ in range((c.n - 1) // 2):
        power = field_power(c, power, 2)
        root ^= power
    if c.fm(root, root) ^ root != rhs:
        return []
    y = c.fm(x, root)
    return [(x, y), (x, x ^ y)]


def span(basis):
    require(len(basis) <= 20, 'base reconstruction exceeds declared verifier resource bound')
    points = [0]
    for value in basis:
        points += [p ^ value for p in points]
    require(len(set(points)) == len(points), 'dependent binary basis')
    return points


def x_orbit(c, x):
    result = set()
    for _ in range(c.n):
        result.add(x)
        x = c.fm(x, x)
    return result


def exact_fields(recipe, fields):
    require(type(recipe) is dict and set(recipe) == {'kind', *fields}, 'unknown/missing base recipe fields')


def construct(c, recipe, depth=0):
    require(depth <= 8, 'base recipe nesting exceeds verifier bound')
    require(type(recipe) is dict, 'invalid base recipe')
    kind = recipe.get('kind')
    detail = {'family': kind}
    if kind in ('factor', 'divisor'):
        all_factors = complete_factors(c.n)
        available = [p for p in all_factors if p.bit_length() - 1 <= 24]
        detail.update(factor_scan_max_degree=24,
                      factorization_complete=len(available) == len(all_factors),
                      omitted_factor_degrees=[p.bit_length() - 1 for p in all_factors if p not in available])
        if kind == 'factor':
            exact_fields(recipe, ('index',))
            index = uint(recipe['index'], 64, 'factor index')
            largest = max(p.bit_length() for p in available)
            choices = [p for p in available if p.bit_length() == largest]
            require(index < len(choices), 'factor index outside enumerated family')
            polynomial = choices[index]
        else:
            exact_fields(recipe, ('indices',))
            indices = recipe['indices']
            require(type(indices) is list and len(set(indices)) == len(indices), 'duplicate divisor index')
            polynomial = 1
            for index in indices:
                require(uint(index, 64, 'divisor index') < len(available), 'divisor index outside enumeration')
                polynomial = poly_product(polynomial, available[index])
        dimension = polynomial.bit_length() - 1
        require(0 < dimension < c.n, 'improper invariant divisor')
        columns = []
        for bit in range(c.n):
            value, result = 1 << bit, 0
            for i in range(dimension + 1):
                if polynomial >> i & 1:
                    result ^= value
                value = c.fm(value, value)
            columns.append(result)
        basis = nullspace(columns, c.n)
        require(len(basis) == dimension, 'invariant kernel dimension mismatch')
        xs = span(basis)
        detail.update(divisor_polynomial=polynomial, nominal_dimension=dimension, basis=basis)
    elif kind == 'frobenius_union':
        exact_fields(recipe, ('seed_masks',))
        basis = recipe['seed_masks']
        require(type(basis) is list and 1 <= len(basis) <= 12, 'invalid union seed dimension')
        for value in basis:
            uint(value, c.n, 'union seed mask')
        xs = sorted(set().union(*(x_orbit(c, x) for x in span(basis))))
        detail.update(seed_dimension=len(basis), nominal_dimension=c.n)
    elif kind == 'subgroup_orbits':
        exact_fields(recipe, ('seed', 'points'))
        requested = uint(recipe['points'], 64, 'requested base points')
        require(0 < requested <= MAX_ABSCISSAE, 'sampled base exceeds verifier resource bound')
        rng = StdRng08.seed_from_u64(recipe['seed'])
        seen, abscissae, drawn, batches = set(), set(), 0, 0
        while True:
            added = 0
            while added < 8:
                drawn += 1
                require(drawn <= 1024 * requested, 'sampled base exhausted its declared budget')
                points = lifts(c, rng.nonzero_below(1 << c.n))
                if not points:
                    continue
                point = c.mul(points[0], c.h)
                if point is None or point[0] in seen:
                    continue
                seen.add(point[0])
                abscissae.update(x_orbit(c, point[0]))
                added += 1
            batches += 1
            if 2 * len(abscissae) >= requested:
                xs = sorted(abscissae)
                break
        detail.update(requested_points=requested, sampling='StdRng08-uniform-nonzero-abscissa',
                      projection='cofactor-multiply-first-half-trace-lift', representative_batch=8,
                      draws=drawn, batches=batches, nominal_dimension=c.n)
    elif kind in ('two_torsion_saturated', 'pruned'):
        exact_fields(recipe, ('parent',) if kind == 'two_torsion_saturated'
                     else ('parent', 'retained_abscissa_orbits'))
        parent, parent_detail = construct(c, recipe['parent'], depth + 1)
        detail['parent'] = parent_detail
        if kind == 'two_torsion_saturated':
            require(c.mul((0, 1), c.h) is None, 'torsion saturation needs an even cofactor')
            shifted = [c.add(p, (0, 1)) for p in parent]
            xs = sorted({p[0] for p in parent + shifted if p is not None})
            detail['nominal_dimension'] = c.n
        else:
            retained = recipe['retained_abscissa_orbits']
            require(type(retained) is list and retained, 'empty orbit restriction')
            available = {min(x_orbit(c, p[0])) for p in parent}
            for value in retained:
                uint(value, c.n, 'retained abscissa')
                require(value in available, 'pruned orbit absent from parent')
            xs = sorted(set().union(*(x_orbit(c, value) for value in retained)))
            detail['nominal_dimension'] = parent_detail['nominal_dimension']
    else:
        require(False, 'unknown factor-base construction')
    require(0 < len(xs) <= MAX_ABSCISSAE, 'base domain exceeds verifier resource bound')
    points = [point for x in xs for point in lifts(c, x)]
    require(points and len(set(points)) == len(points), 'empty or duplicated reconstructed base')
    detail['domain_abscissae'] = len(xs)
    return points, detail


def effective_recipe(job, fixture):
    cfg = job['config']
    orbits, cube = cfg.get('factor_base_orbits'), cfg.get('factor_base_cube_root', False)
    require(type(cube) is bool, 'invalid cube-root policy')
    explicit = cfg.get('factor_base')
    require(not (orbits is not None and cube) and
            not (explicit is not None and (orbits is not None or cube)), 'conflicting base policies')
    if orbits is not None or cube:
        parent = job['factor_base']
        require(parent['kind'] == 'subgroup_orbits', 'orbit policy requires sampled subgroup base')
        if orbits is not None:
            require(1 <= uint(orbits, 64, 'orbit policy') <= 8, 'invalid orbit policy')
            points = 2 * fixture['degree'] * orbits
        else:
            bound, points = int(fixture['subgroup_order']) // 2, 1
            while points**3 < bound:
                points += 1
            points = max(points, 2 * fixture['degree'])
        return dict(kind='subgroup_orbits', seed=parent['seed'], points=points)
    return copy.deepcopy(explicit if explicit is not None else job['factor_base'])


@lru_cache(maxsize=128)
def _reconstruct(recipe_json, fixture_json):
    points, detail = construct(Curve(json.loads(fixture_json)), json.loads(recipe_json))
    return tuple(points), json.dumps(detail, sort_keys=True)


def verify_base(report, fixture, job):
    require(job['degree'] == fixture['degree'] and job['curve_a'] == fixture['curve_a'], 'base curve mismatch')
    require(5 <= fixture['degree'] <= 31 and fixture['degree'] % 2 == 1, 'generic base degree unsupported')
    recipe = effective_recipe(job, fixture)
    require(report.get('effective_factor_base') == recipe, 'reported base recipe differs from job')
    points, detail = _reconstruct(json.dumps(recipe, sort_keys=True), json.dumps(fixture, sort_keys=True))
    encoded = [[str(x), str(y)] for x, y in points]
    require(report['factor_base'] == encoded, 'factor-base construction or ordering mismatch')
    inventory = factor_base_inventory(report, fixture)
    return dict(schema_version=1, recipe=recipe, construction=json.loads(detail), inventory=inventory,
                recipe_sha256=sha256(recipe), scope='independent factor-base construction and census',
                promotion_eligible=False)
