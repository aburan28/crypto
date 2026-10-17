"""Complete bounded all-x S4 solver using compact norm fibers and twist closure.

Outputs ordered source abscissa triples in D with S4(x1,x2,x3,z)=0.
Rational norm triples alone are insufficient: Frobenius changes the sign of
nonrational points above K-abscissae. Their contribution must be rational
two-torsion when the chosen target is rational. Exact twist-pair closure
supplies mixed triples; at a two-torsion target twist triples supply the
all-nonrational branch. An initially nonrational target is made rational by
an exact quadratic twist. All arithmetic and checks stay in K.

Finite points only. No claim of novelty or measured Gröbner speedup.
"""

from itertools import permutations, product
import random

from rootfinder import polynomial as poly
from lab.research.odd_family import (prepare_model, source_sums, domain_contains,
    conditional_family, coefficient_maps, cubic_from_parameter)


class _Unresolved(Exception):
    pass


def _count(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f'{name} must be a nonnegative integer')
    return value


def on_curve(field, curve, point):
    if point is None:
        return True
    x, y = point
    return field.square(y) == poly.evaluate_poly((curve['a6'], curve['a4'], curve['a2'], 1), x, field)


def negate(field, point):
    return None if point is None else (point[0], field.neg(point[1]))


def add_points(field, curve, left, right):
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u:
        if field.add(y, v) == 0:
            return None
        if y != v:
            raise ArithmeticError('inconsistent equal-abscissa curve points')
        derivative = field.add(field.add(field.mul(3 % field.characteristic, field.square(x)),
                                         field.mul(field.mul(2 % field.characteristic, curve['a2']), x)),
                               curve['a4'])
        slope = field.mul(derivative, field.inv(field.mul(2 % field.characteristic, y)))
    else:
        slope = field.mul(field.sub(v, y), field.inv(field.sub(u, x)))
    w = field.sub(field.sub(field.sub(field.square(slope), curve['a2']), x), u)
    result = (w, field.sub(field.mul(slope, field.sub(x, w)), y))
    if not on_curve(field, curve, result):
        raise ArithmeticError('group-law replay produced an off-curve point')
    return result


def scalar_multiply(field, curve, point, scalar):
    _count(scalar, 'scalar')
    total = None
    while scalar:
        if scalar & 1:
            total = add_points(field, curve, total, point)
        scalar >>= 1
        if scalar:
            point = add_points(field, curve, point, point)
    return total


def original_s4_value(field, curve, x1, x2, x3, target_x):
    """Direct degree-two resultant verification in original field coordinates."""
    def coefficients(x, y):
        two, four = 2 % field.characteristic, 4 % field.characteristic
        xy, total = field.mul(x, y), field.add(x, y)
        a = field.square(field.sub(x, y))
        inside = field.add(field.mul(total, field.add(xy, curve['a4'])),
            field.add(field.mul(field.mul(two, curve['a2']), xy), field.mul(two, curve['a6'])))
        b = field.neg(field.mul(two, inside))
        c = field.sub(field.square(field.sub(curve['a4'], xy)),
            field.mul(field.mul(four, curve['a6']), field.add(total, curve['a2'])))
        return a, b, c
    a, b, c = coefficients(x1, x2)
    d, e, g = coefficients(x3, target_x)
    first = field.sub(field.mul(a, g), field.mul(c, d))
    second = field.sub(field.mul(a, e), field.mul(b, d))
    third = field.sub(field.mul(b, g), field.mul(c, e))
    return field.sub(field.square(first), field.mul(second, third))


def _domain_parameters(model, scale=1):
    f, domain = model['field'], model['domain']
    result = {'offset': f.mul(scale, model['offset'])}
    if domain['kind'] == 'interval':
        result.update(interval_length=domain['size'], interval_step=f.mul(scale, domain['step']),
                      max_interval=domain['size'])
    else:
        result['basis'] = tuple(f.mul(scale, value) for value in model['basis'])
    return result


def domain_values(model):
    f, domain = model['field'], model['domain']
    if domain['kind'] == 'interval':
        for integer in range(domain['size']):
            yield f.add(model['offset'], f.mul(domain['step'], integer))
    else:
        for coordinates in product(range(f.characteristic), repeat=len(model['basis'])):
            value = model['offset']
            for coefficient, basis in zip(coordinates, model['basis']):
                value = f.add(value, f.mul(coefficient, basis))
            yield value


def twist_curve(field, curve, nonsquare):
    square = field.square(nonsquare)
    return {'a2': field.mul(nonsquare, curve['a2']),
            'a4': field.mul(square, curve['a4']),
            'a6': field.mul(field.mul(square, nonsquare), curve['a6'])}


def solve_all_x(model, *, max_sums=1000, max_domain_scan=1000,
                max_parameter_degree=512, max_scalar_products=5_000_000,
                max_split_trials=20_000, max_nonresidue_trials=128,
                max_unique_solutions=100_000, max_output=10_000, seed=0):
    """Every ordered K-abscissa triple in D satisfying the geometric S4 contract.

    The desired target abscissa is model.fixed_point.x. Either sign of its
    ordinate gives the same x-only root set. A supplied finite rational target
    is required by prepare_model; model adaptation for nonsquare target C(z)
    is supplied separately by prepare_x_model below. All family, splitting,
    scan and memory limits cause explicit unresolved output, never exact UNSAT.
    """
    for name, value in [('max_sums', max_sums), ('max_domain_scan', max_domain_scan),
                        ('max_split_trials', max_split_trials), ('max_nonresidue_trials', max_nonresidue_trials),
                        ('max_unique_solutions', max_unique_solutions), ('max_output', max_output), ('seed', seed)]:
        _count(value, name)
    f, curve = model['field'], model['curve']
    target = negate(f, model['fixed_point'])
    scale_to_original = model.get('scale_to_original', 1)
    inverse_scale = f.inv(scale_to_original)
    stats = {'sum_fibers': 0, 'parameter_candidates': 0, 'repeated_parameter_degree': 12,
             'split_trials': 0, 'nonresidue_trials': 0, 'mixed_domain_elements_scanned': 0,
             'point_triples_replayed': 0, 'field_materialized': False,
             'maximum_actual_parameter_degree': -1, 'quotient_scalar_products': 0}
    result = {'schema': 'odd_curve_all_x_s4_result_v1', 'status': 'unresolved',
              'scope': 'all ordered source abscissa triples of the S4 finite-point summation contract',
              'solutions': [], 'solution_count_exact': None, 'solution_count_lower_bound': 0,
              'search_complete': False, 'output_complete': False, 'verified': True,
              'target_x': f.mul(inverse_scale, target[0]), 'stats': stats, 'reason': None}
    roots = set()
    root_cache = {}

    def all_roots(polynomial):
        polynomial = poly.trim(polynomial, f)
        if polynomial in root_cache:
            return root_cache[polynomial]
        report = poly.field_roots(polynomial, f, seed=seed,
                                  max_trials=max_split_trials-stats['split_trials'])
        stats['split_trials'] += report['split_trials']
        if report['status'] != 'complete':
            raise _Unresolved('field-root splitting budget exhausted')
        values = tuple(report['roots'])
        root_cache[polynomial] = values
        return values

    def emit(xs):
        xs = tuple(f.mul(inverse_scale, x) for x in xs)
        if not all(domain_contains(model, f.mul(scale_to_original, x)) for x in xs):
            raise ArithmeticError('mapped source abscissa left original compact domain')
        if original_s4_value(f, model.get('original_curve', curve), *xs, result['target_x']):
            raise ArithmeticError('mapped all-x root failed original S4 resultant verification')
        if xs not in roots:
            if len(roots) >= max_unique_solutions:
                raise _Unresolved('unique all-x root storage budget exhausted')
            roots.add(xs)

    def rational_branch(branch_model, desired, branch_scale=1):
        if branch_model['domain']['sum_size'] > max_sums:
            raise _Unresolved('exact source-sum cover exceeds max_sums')
        inverse_branch = f.inv(branch_scale)
        for s in source_sums(branch_model):
            stats['sum_fibers'] += 1
            family = conditional_family(branch_model, s,
                max_parameter_degree=max_parameter_degree, max_scalar_products=max_scalar_products,
                max_split_trials=max_split_trials-stats['split_trials'],
                max_parameters=2*branch_model['domain']['size']+12, seed=seed)
            stats['split_trials'] += family['stats']['split_trials']
            if family['status'] != 'complete':
                raise _Unresolved('compact norm family unresolved: ' + str(family['reason']))
            stats['parameter_candidates'] += len(family['parameters'])
            stats['maximum_actual_parameter_degree'] = max(stats['maximum_actual_parameter_degree'],
                *family['stats'].get('actual_parameter_degrees', [-1]))
            stats['quotient_scalar_products'] += family['stats'].get('quotient_scalar_products', 0)
            maps = coefficient_maps(branch_model, s)
            for parameter in family['parameters']:
                cubic = cubic_from_parameter(branch_model, s, parameter)
                residue = cubic
                multiset = []
                for x in all_roots(cubic):
                    while len(residue) > 1 and poly.evaluate_poly(residue, x, f) == 0:
                        residue, remainder = poly.divmod_poly(residue, (f.neg(x), 1), f)
                        if remainder:
                            raise ArithmeticError('root deflation failed')
                        multiset.append(x)
                if residue != (1,) or len(multiset) != 3:
                    raise ArithmeticError('supported norm cubic did not split fully with multiplicity')
                if parameter:
                    c, d = (poly.evaluate_poly(maps[name], parameter, f) for name in ('c', 'd'))
                    inverse_parameter = f.inv(parameter)
                    fibers = {x: (f.neg(f.mul(f.add(f.add(f.square(x), f.mul(c, x)), d), inverse_parameter)),)
                              for x in set(multiset)}
                else:
                    fibers = {x: all_roots((f.neg(poly.evaluate_poly(
                        (branch_model['curve']['a6'], branch_model['curve']['a4'],
                         branch_model['curve']['a2'], 1), x, f)), 0, 1))
                              for x in set(multiset)}
                for xs in set(permutations(multiset)):
                    for ys in product(*(fibers[x] for x in xs)):
                        points = tuple(zip(xs, ys))
                        if not all(on_curve(f, branch_model['curve'], point) for point in points):
                            raise ArithmeticError('norm ordinate reconstruction failed')
                        total = None
                        for point in points:
                            total = add_points(f, branch_model['curve'], total, point)
                        stats['point_triples_replayed'] += 1
                        if total == desired:
                            emit(tuple(f.mul(inverse_branch, x) for x in xs))

    try:
        if model['domain']['size'] > max_domain_scan:
            raise _Unresolved('mixed-twist source-domain scan exceeds max_domain_scan')
        rational_branch(model, target)
        # Find an exactly verified nonsquare; random discovery is capped and
        # does not contribute any heuristic root-completeness assumption.
        rng = random.Random(seed ^ 0x4F44445457495354)
        nonsquare = None
        minus_one = f.neg(1)
        for _ in range(max_nonresidue_trials):
            stats['nonresidue_trials'] += 1
            candidate = rng.randrange(1, f.order)
            if f.pow(candidate, (f.order-1)//2) == minus_one:
                nonsquare = candidate
                break
        if nonsquare is None:
            raise _Unresolved('verified quadratic nonresidue search budget exhausted')
        opposite_curve = twist_curve(f, curve, nonsquare)
        target_opposite = (f.mul(nonsquare, target[0]), 0) if target[1] == 0 else None
        if target_opposite is not None:
            opposite_model = prepare_model(f, opposite_curve, target_opposite,
                                            **_domain_parameters(model, nonsquare))
            rational_branch(opposite_model, target_opposite, nonsquare)
        two_torsion = [None] + [(x, 0) for x in all_roots((curve['a6'], curve['a4'], curve['a2'], 1))]
        inverse_nonsquare = f.inv(nonsquare)
        for q in two_torsion:
            singleton = add_points(f, curve, target, negate(f, q))
            if singleton is None or not domain_contains(model, singleton[0]):
                continue
            opposite_q = None if q is None else (f.mul(nonsquare, q[0]), 0)
            for x in domain_values(model):
                stats['mixed_domain_elements_scanned'] += 1
                rhs = poly.evaluate_poly((curve['a6'], curve['a4'], curve['a2'], 1), x, f)
                if not rhs or f.pow(rhs, (f.order-1)//2) != minus_one:
                    continue
                if opposite_q is None:
                    mate_x = x
                else:
                    opposite_x = f.mul(nonsquare, x)
                    opposite_rhs = poly.evaluate_poly((opposite_curve['a6'], opposite_curve['a4'],
                                                        opposite_curve['a2'], 1), opposite_x, f)
                    ys = all_roots((f.neg(opposite_rhs), 0, 1))
                    if not ys:
                        raise ArithmeticError('verified nonrational source failed twist lifting')
                    point = (opposite_x, ys[0])
                    mate = add_points(f, opposite_curve, opposite_q, negate(f, point))
                    if mate is None or mate[1] == 0:
                        continue
                    mate_x = f.mul(inverse_nonsquare, mate[0])
                    if not domain_contains(model, mate_x):
                        continue
                for xs in set(permutations((x, mate_x, singleton[0]))):
                    emit(xs)
        result.update(search_complete=True, solution_count_exact=len(roots),
                      output_complete=len(roots) <= max_output,
                      status='complete' if len(roots) <= max_output else 'output_truncated',
                      reason='All compact norm fibers, repeated fibers and Frobenius twist closures exhausted')
    except _Unresolved as error:
        result['reason'] = str(error)
    result['solution_count_lower_bound'] = len(roots)
    result['solutions'] = [list(xs) for xs in sorted(roots)[:max_output]]
    return result


def prepare_x_model(field, curve, target_x, *, offset=0, basis=(), interval_length=None,
                    interval_step=1, max_interval=256, seed=0, max_trials=256):
    """Make an arbitrary finite target abscissa rational by an exact twist.

    Initial point/nonresidue discovery is verified and explicitly bounded.
    Caller must include this construction inside the complete solver interval.
    """
    from research.advance.general_chain import GeneralS3Fiber
    checked = GeneralS3Fiber(field, curve)
    if field.characteristic == 2:
        raise ValueError('odd characteristic required')
    curve = checked.curve
    original_curve = dict(curve)
    z = field.validate(target_x)
    rhs = poly.evaluate_poly((curve['a6'], curve['a4'], curve['a2'], 1), z, field)
    roots = poly.field_roots((field.neg(rhs), 0, 1), field, seed=seed, max_trials=max_trials)
    if roots['status'] != 'complete':
        raise _Unresolved('initial target ordinate splitting unresolved')
    scale = 1
    parameters = {'offset': offset, 'basis': basis, 'interval_length': interval_length,
                  'interval_step': interval_step, 'max_interval': max_interval}
    if not roots['roots']:
        rng = random.Random(seed ^ 0x544152474554)
        scale = None
        for _ in range(max_trials):
            candidate = rng.randrange(1, field.order)
            if field.pow(candidate, (field.order-1)//2) == field.neg(1):
                scale = candidate
                break
        if scale is None:
            raise _Unresolved('initial target twist nonresidue search unresolved')
        curve = twist_curve(field, curve, scale)
        z = field.mul(scale, z)
        parameters['offset'] = field.mul(scale, offset)
        parameters['basis'] = tuple(field.mul(scale, value) for value in basis)
        parameters['interval_step'] = field.mul(scale, interval_step)
        rhs = poly.evaluate_poly((curve['a6'], curve['a4'], curve['a2'], 1), z, field)
        roots = poly.field_roots((field.neg(rhs), 0, 1), field, seed=seed, max_trials=max_trials)
        if roots['status'] != 'complete' or not roots['roots']:
            raise _Unresolved('twisted target ordinate splitting unresolved')
    model = prepare_model(field, curve, (z, field.neg(roots['roots'][0])), **parameters)
    model['scale_to_original'] = scale
    model['original_curve'] = original_curve
    return model
