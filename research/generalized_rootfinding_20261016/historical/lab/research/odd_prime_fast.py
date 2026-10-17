"""FLINT prime-field backend for exact compact-domain odd S4 parameter fibers.

A balanced product of interval leaves X-d in the monic cubic quotient
K[B][X]/H_s replaces scalar Horner evaluation. Native fmpz_mod_poly products
and gcd keep parameter degree <=2k. The completed full all-x algorithm is
reused through copies of function globals; frozen modules are never mutated.

Optional certified native roots use gcd(P,B^p-B), verified distinct membership
and equality with the gcd degree. FLINT factorization calls have an explicit
budget; internal random splitting trials are not exposed by FLINT. External
process wall/memory limits must therefore bound each complete invocation.
No novelty or timing claim is made by this implementation.
"""

from itertools import permutations
from pathlib import Path
import sys
from types import FunctionType, SimpleNamespace

from rootfinder import polynomial as reference_poly
from research.advance.general_chain import GeneralS3Fiber
from lab.research import odd_family, odd_solver


_CONTEXTS = {}


def _flint():
    try:
        import flint
    except ImportError:
        local = Path(__file__).resolve().parents[2] / 'vendor/python-flint'
        if not local.is_dir():
            raise ImportError('Install pinned python-flint0.8.0 in project-local vendor/python-flint')
        sys.path.insert(0, str(local))
        import flint
    if flint.__version__ != '0.8.0':
        raise ImportError('The preregistered native backend requires python-flint0.8.0')
    flint.ctx.threads = 1
    return flint


def _context(field):
    if field.degree != 1 or field.characteristic != field.order or field.characteristic <= 2:
        raise ValueError('fast interval backend requires an odd prime field')
    p = field.characteristic
    if p not in _CONTEXTS:
        flint = _flint()
        # This is the proved-prime interface, not is_probable_prime. The
        # context proof is charged once per prime per cold process.
        if not flint.fmpz(p).is_prime():
            raise ValueError('exact FLINT primality proof rejected the parameter modulus')
        _CONTEXTS[p] = flint.fmpz_mod_poly_ctx(p)
    return _CONTEXTS[p]


def dependency_metadata():
    flint = _flint()
    return {'python_flint': flint.__version__, 'flint': flint.__FLINT_VERSION__,
            'prime_proof_api': 'fmpz_is_prime',
            'primary_proof_documentation': 'https://flintlib.org/doc/fmpz.html#c.fmpz_is_prime',
            'polynomial_api': 'fmpz_mod_poly_ctx; multiplication, gcd, pow_mod, roots',
            'native_threads': flint.ctx.threads}


def _tuple(polynomial):
    return tuple(int(coefficient) for coefficient in polynomial.coeffs())


def _clone(function, replacements):
    namespace = dict(function.__globals__)
    namespace.update(replacements)
    copied = FunctionType(function.__code__, namespace, function.__name__,
                          function.__defaults__, function.__closure__)
    copied.__kwdefaults__ = dict(function.__kwdefaults__ or {})
    copied.__doc__ = function.__doc__
    return copied


def _product_tree(values, multiply, identity):
    layer = list(values)
    if not layer:
        return identity, 0
    rounds = 0
    while len(layer) > 1:
        next_layer = []
        for index in range(0, len(layer)-1, 2):
            next_layer.append(multiply(layer[index], layer[index+1]))
        if len(layer) & 1:
            next_layer.append(layer[-1])
        layer = next_layer
        rounds += 1
    return layer[0], rounds


def prepare_prime_model(field, curve, fixed_point, *, offset=0, basis=(),
                        interval_length=None, interval_step=1, max_interval=256):
    """Build the exact interval polynomial with a balanced native product tree."""
    if interval_length is None:
        return odd_family.prepare_model(field, curve, fixed_point, offset=offset, basis=basis,
                                        interval_step=interval_step, max_interval=max_interval)
    if basis:
        raise ValueError('prime interval backend cannot combine an interval with a subspace basis')
    ctx = _context(field)
    checked = GeneralS3Fiber(field, curve)
    curve = dict(checked.curve)
    offset, interval_step = field.validate(offset), field.validate(interval_step)
    odd_family._count(interval_length, 'interval_length')
    odd_family._count(max_interval, 'max_interval')
    if not interval_step or not 1 <= interval_length <= min(field.characteristic, max_interval):
        raise ValueError('interval size/step must fit the distinct prime-field domain and cap')
    fixed_point = tuple(fixed_point)
    if len(fixed_point) != 2:
        raise ValueError('fixed_point must be a coordinate pair')
    z, y = (field.validate(value) for value in fixed_point)
    if field.square(y) != reference_poly.evaluate_poly((curve['a6'], curve['a4'], curve['a2'], 1), z, field):
        raise ValueError('fixed_point must lie on the completed curve')
    p = field.characteristic
    polynomial, rounds = _product_tree(
        (ctx([-(offset + interval_step*integer) % p, 1]) for integer in range(interval_length)),
        lambda left, right: left*right, ctx.one())
    return {'field': field, 'curve': curve, 'fixed_point': (z, y), 'offset': offset, 'basis': (),
            'domain': {'kind': 'interval', 'size': interval_length,
                       'sum_size': min(p, 3*interval_length-2), 'polynomial': _tuple(polynomial),
                       'step': interval_step, 'inverse_step': field.inv(interval_step)},
            '_flint_domain_polynomial': polynomial,
            'fast_preparation': {'domain_leaves_built': interval_length,
                                 'domain_product_tree_rounds': rounds,
                                 'field_materialized': False, **dependency_metadata()}}


def prepare_x_prime_model(field, curve, target_x, **parameters):
    """Original-coordinate target adaptation using the fast interval setup."""
    return _clone(odd_solver.prepare_x_model, {'prepare_model': prepare_prime_model})(
        field, curve, target_x, **parameters)


class _Quotient:
    def __init__(self, ctx, s, maps, max_degree, max_dense_products):
        self.ctx, self.s = ctx, s
        self.t, self.u = ctx(list(maps['e2'])), ctx(list(maps['e3']))
        self.max_degree, self.max_dense_products = max_degree, max_dense_products
        self.dense_products = 0
        self.polynomial_multiplications = 0
        self.quotient_multiplications = 0

    def mul(self, a, b):
        if a.is_zero() or b.is_zero():
            return self.ctx.zero()
        if a.degree()+b.degree() > self.max_degree:
            raise odd_family._Limit('native parameter intermediate degree exceeds guard')
        bound = (a.degree()+1)*(b.degree()+1)
        if self.dense_products + bound > self.max_dense_products:
            raise odd_family._Limit('native construction dense scalar-product upper bound exceeds guard')
        self.dense_products += bound
        self.polynomial_multiplications += 1
        return a*b

    def product(self, left, right):
        self.quotient_multiplications += 1
        values = [self.ctx.zero() for _ in range(5)]
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                values[i+j] += self.mul(a, b)
        for degree in (4, 3):
            value = values[degree]
            if not value.is_zero():
                values[degree-1] += value*self.s
                values[degree-2] -= self.mul(value, self.t)
                values[degree-3] += self.mul(value, self.u)
        return tuple(values[:3])


def parameter_remainder_fast(model, source_sum, *, max_parameter_degree=512,
                             max_scalar_products=5_000_000):
    """Balanced cubic-quotient interval product, retaining the exact2k bound."""
    if model['domain']['kind'] != 'interval':
        return odd_family.parameter_domain_remainder(model, source_sum,
            max_parameter_degree=max_parameter_degree, max_scalar_products=max_scalar_products)
    f, domain = model['field'], model['domain']
    ctx = _context(f)
    s = f.validate(source_sum)
    odd_family._count(max_parameter_degree, 'max_parameter_degree')
    odd_family._count(max_scalar_products, 'max_scalar_products')
    k = domain['size']
    if 2*k > max_parameter_degree:
        raise odd_family._Limit(f'certified2k parameter-degree bound {2*k} exceeds max_parameter_degree')
    maps = odd_family.coefficient_maps(model, s)
    quotient = _Quotient(ctx, s, maps, max_parameter_degree+16, max_scalar_products)
    p, step, offset = f.characteristic, domain['step'], model['offset']
    residue, rounds = _product_tree(
        ((ctx([-(offset+step*integer) % p]), ctx.one(), ctx.zero()) for integer in range(k)),
        quotient.product, (ctx.one(), ctx.zero(), ctx.zero()))
    degrees = [row.degree() for row in residue]
    if any(degree > 2*k for degree in degrees):
        raise ArithmeticError('balanced quotient violated the proved2k parameter degree bound')
    return {'coefficients': tuple(_tuple(row) for row in residue),
            'coefficient_maps': maps, '_native_coefficients': residue,
            'stats': {'domain_degree': k, 'certified_parameter_degree_bound': 2*k,
                      'actual_parameter_degrees': degrees,
                      'quotient_scalar_products': quotient.dense_products,
                      'quotient_scalar_product_scope': 'dense schoolbook upper bound for tracked native polynomial products only; scalar scaling, addition, setup, gcd and roots excluded',
                      'native_polynomial_multiplications': quotient.polynomial_multiplications,
                      'native_quotient_multiplications': quotient.quotient_multiplications,
                      'balanced_product_tree_rounds': rounds,
                      'interval_domain_leaves_built': k,
                      'domain_materialized': False, 'field_materialized': False,
                      'interval_roots_enumerated_in_preprocessing': True,
                      **dependency_metadata()}}


def _native_domain_remainder(model, cubic):
    ctx = _context(model['field'])
    polynomial = model.get('_flint_domain_polynomial')
    if polynomial is None:
        polynomial = ctx(list(model['domain']['polynomial']))
    return polynomial % ctx(list(cubic))


class _NativeRoots:
    def __init__(self, max_calls):
        self.max_calls = odd_family._count(max_calls, 'max_native_factor_calls')
        self.calls = 0

    def roots(self, coefficients, field, seed=0, max_trials=128):
        """Exact native root membership and Frobenius-degree completeness."""
        ctx = _context(field)
        polynomial = ctx(list(coefficients))
        if polynomial.is_zero():
            return {'status': 'unrestricted', 'roots': [], 'restricted_degree': field.order,
                    'split_trials': 0, 'roots_enumerated': False}
        if polynomial.degree() == 0:
            return {'status': 'complete', 'roots': [], 'restricted_degree': 0, 'split_trials': 0}
        variable = ctx.gen()
        restricted = polynomial.gcd(variable.pow_mod(field.order, polynomial)-variable)
        count = restricted.degree()
        if count == 0:
            values = []
        elif count == 1:
            values = [int(-restricted[0]/restricted[1])]
        else:
            if self.calls >= self.max_calls:
                raise odd_solver._Unresolved('native finite-field factorization call budget exhausted')
            self.calls += 1
            values = sorted(int(root) for root in restricted.roots(multiplicities=False))
        if len(set(values)) != len(values) or len(values) != count:
            raise ArithmeticError('native root list and exact Frobenius gcd degree disagree')
        raw = tuple(int(coefficient) for coefficient in polynomial.coeffs())
        if any(reference_poly.evaluate_poly(raw, field.validate(root), field) != 0 for root in values):
            raise ArithmeticError('native root failed independent field-protocol membership verification')
        return {'status': 'complete', 'roots': values, 'restricted_degree': count,
                'degree': polynomial.degree(), 'split_trials': 0,
                'root_backend': 'FLINT Frobenius-gcd-count-certified',
                'native_factor_calls': self.calls,
                'trial_budget_scope': 'native internal splitting trials not exposed; separately capped factorization calls and external process deadline required'}


def conditional_family_fast(model, source_sum, *, max_parameter_degree=512,
                            max_scalar_products=5_000_000, max_split_trials=256,
                            max_parameters=10_000, seed=0, root_finder=None, include_repeated=True):
    """Same exact support-family contract, native tree/gcd, configurable roots."""
    if model['domain']['kind'] != 'interval':
        return odd_family.conditional_family(model, source_sum,
            max_parameter_degree=max_parameter_degree, max_scalar_products=max_scalar_products,
            max_split_trials=max_split_trials, max_parameters=max_parameters, seed=seed)
    f, s = model['field'], model['field'].validate(source_sum)
    for name, value in [('max_parameters', max_parameters), ('max_split_trials', max_split_trials), ('seed', seed)]:
        odd_family._count(value, name)
    ctx = _context(f)
    root_finder = root_finder or reference_poly.field_roots
    result = {'schema': 'odd_curve_exact_domain_family_v1', 'status': 'unresolved',
              'scope': ('all norm parameters whose monic cubic has all roots in the prime interval domain'
                        if include_repeated else 'distinct-root norm parameters only; repeated parameters delegated to certified global closure'),
              'source_sum': s, 'parameters': [], 'solution_count_exact': None,
              'search_complete': False, 'output_complete': False, 'verified': True,
              'stats': {'split_trials': 0, 'parameter_backend': 'FLINT balanced cubic quotient'}, 'reason': None}
    if not odd_family._sum_admitted(model, s):
        result.update(status='complete', solution_count_exact=0, search_complete=True, output_complete=True)
        return result
    try:
        construction = parameter_remainder_fast(model, s, max_parameter_degree=max_parameter_degree,
                                                max_scalar_products=max_scalar_products)
        result['stats'].update(construction['stats'])
        common = ctx.zero()
        for row in construction['_native_coefficients']:
            common = common.gcd(row)
            if common.degree() == 0:
                break
        if common.is_zero():
            raise ArithmeticError('all formal supported-cubic constraints cannot vanish identically')
        delta = ctx(list(odd_family.discriminant(model, s, construction['coefficient_maps'])))
        result['stats']['simple_parameter_polynomial_degree'] = common.degree()
        result['stats']['repeated_parameter_polynomial_degree'] = delta.degree()
        admitted = set()
        branches = [('distinct', common)] + ([('repeated', delta)] if include_repeated else [])
        result['stats']['repeated_parameter_branch_included'] = include_repeated
        for branch, polynomial in branches:
            report = root_finder(_tuple(polynomial), f, seed=seed,
                max_trials=max_split_trials-result['stats']['split_trials'])
            result['stats']['split_trials'] += report['split_trials']
            if report['status'] != 'complete':
                raise odd_family._Limit(branch+' parameter root splitting unresolved')
            for parameter in report['roots']:
                cubic = odd_family.cubic_from_parameter(model, s, parameter)
                remainder = _native_domain_remainder(model, cubic)
                if branch == 'distinct':
                    if not remainder.is_zero() or int(delta(parameter)) == 0:
                        raise ArithmeticError('distinct native parameter failed support/multiplicity replay')
                    admitted.add(parameter)
                elif remainder.pow_mod(3, ctx(list(cubic))).is_zero():
                    admitted.add(parameter)
        values = sorted(admitted)
        result.update(status='complete' if len(values) <= max_parameters else 'output_truncated',
                      parameters=values[:max_parameters], solution_count_exact=len(values),
                      search_complete=True, output_complete=len(values) <= max_parameters,
                      reason='Every native exact compact-domain fiber exhausted' +
                             (' including degree-twelve repeated fibers' if include_repeated else '; repeated fibers delegated to a separately certified global closure'))
    except (odd_family._Limit, odd_solver._Unresolved) as error:
        result['reason'] = str(error)
    return result


def solve_prime_all_x(model, *, max_native_factor_calls=20_000, native_roots=True,
                      batch_repeated=True, **limits):
    """Complete S4 with native arithmetic and guarded global repeated closure.

    If the curve has no rational two-torsion and the target abscissa lies
    outside D, all repeated roots are P+P+(T-2P) for rational P. An exact
    Frobenius-gcd degree-zero certificate enables one O(k) closure scan instead
    of factoring the degree-twelve discriminant at every source sum.
    Other curves retain the original complete discriminant/twist algorithm.
    """
    if model['domain']['kind'] != 'interval':
        raise ValueError('fast prime all-x entry point requires an interval model')
    f = model['field']
    ctx = _context(f)
    factors = _NativeRoots(max_native_factor_calls)
    finder = factors.roots if native_roots else reference_poly.field_roots
    proxy = SimpleNamespace(**{**vars(reference_poly), 'field_roots': finder})
    target = odd_solver.negate(f, model['fixed_point'])
    no_two_torsion = False
    if batch_repeated and native_roots and target[1] and not odd_family.domain_contains(model, target[0]):
        cubic = ctx([model['curve']['a6'], model['curve']['a4'], model['curve']['a2'], 1])
        variable = ctx.gen()
        no_two_torsion = cubic.gcd(variable.pow_mod(f.order, cubic)-variable).degree() == 0
    can_batch = bool(no_two_torsion)
    def family(*args, **kwargs):
        return conditional_family_fast(*args, root_finder=finder,
                                       include_repeated=not can_batch, **kwargs)
    solver = _clone(odd_solver.solve_all_x, {'conditional_family': family,
        'prepare_model': prepare_prime_model, 'poly': proxy})
    requested_output = limits.get('max_output', 10_000)
    unique_cap = limits.get('max_unique_solutions', 100_000)
    if can_batch:
        # Exact x-only deduplication needs the complete distinct-root set
        # while adding the bulk repeated branch. Retention cap stays separate.
        odd_family._count(requested_output, 'max_output')
        limits['max_output'] = unique_cap
    result = solver(model, **limits)
    closure_scans = 0
    if can_batch:
        values = {tuple(row) for row in result['solutions']}
        inverse_scale = f.inv(model.get('scale_to_original', 1))
        try:
            if result['search_complete']:
                for x in odd_solver.domain_values(model):
                    closure_scans += 1
                    rhs = reference_poly.evaluate_poly((model['curve']['a6'], model['curve']['a4'],
                                                         model['curve']['a2'], 1), x, f)
                    report = finder((f.neg(rhs), 0, 1), f)
                    if report['status'] != 'complete':
                        raise odd_solver._Unresolved('global repeated-source ordinate splitting unresolved')
                    for y in report['roots']:
                        point = (x, y)
                        doubled = odd_solver.add_points(f, model['curve'], point, point)
                        mate = odd_solver.add_points(f, model['curve'], target, odd_solver.negate(f, doubled))
                        if mate is None or not odd_family.domain_contains(model, mate[0]):
                            continue
                        for row in set(permutations((x, x, mate[0]))):
                            row = tuple(f.mul(inverse_scale, value) for value in row)
                            if odd_solver.original_s4_value(f, model.get('original_curve', model['curve']),
                                                           *row, result['target_x']):
                                raise ArithmeticError('bulk repeated root failed original S4 verification')
                            if row not in values:
                                if len(values) >= unique_cap:
                                    raise odd_solver._Unresolved('global repeated-root unique storage budget exhausted')
                                values.add(row)
                result.update(solution_count_exact=len(values),
                    output_complete=len(values) <= requested_output,
                    status='complete' if len(values) <= requested_output else 'output_truncated',
                    reason='All distinct norm fibers and the certified global repeated point closure exhausted')
        except odd_solver._Unresolved as error:
            result.update(status='unresolved', search_complete=False, output_complete=False,
                          solution_count_exact=None, reason=str(error))
        result['solution_count_lower_bound'] = len(values)
        result['solutions'] = [list(row) for row in sorted(values)[:requested_output]]
    result['stats'].update(parameter_backend='FLINT balanced interval quotient/gcd',
        root_backend='FLINT Frobenius-degree-certified' if native_roots else 'reference certified splitting',
        native_factor_calls=factors.calls,
        native_factor_call_budget=factors.max_calls,
        native_internal_split_trials_exposed=False,
        original_modules_mutated=False,
        batch_repeated_closure=can_batch,
        rational_two_torsion_absence_proved_by_frobenius_gcd=no_two_torsion,
        global_repeated_domain_scans=closure_scans,
        **dependency_metadata())
    return result
