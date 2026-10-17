"""Exact compact-domain odd-curve norm-parameter families without Gröbner bases.

For Y²=C(X), F=X²+cX+d+BY vanishes at fixed rational T=(z,y),
d=-z²-cz-By. Write norm(F)/(X-z)=X³-sX²+E2(B)X-E3(B), fixing s.
E2 has degree four and E3 has degree at most four in B. Evaluate the source-domain vanishing
polynomial modulo this cubic over K[B], and gcd the three remainder
coefficients. This is an exact distinct-root support test with parameter
degree at most 2*|D|. Repeated roots use a separate degree-twelve cubic
 discriminant and an exact P_D³ divisibility check. No field enumeration or
source-domain enumeration is used for affine-subspace models.

The result certifies supported norm cubics, not ordered rational point
triples or unrestricted summation-polynomial roots. Point lifting, group
replay and a matched performance experiment are separate subsequent steps.
"""

from itertools import product

from rootfinder import polynomial as poly
from research.advance.general_chain import GeneralS3Fiber


class _Limit(Exception):
    pass


def _count(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f'{name} must be a nonnegative integer')
    return value


def _scalar(field, integer):
    return integer % field.characteristic


def _linearized_value(coefficients, value, field):
    total, power = 0, value
    for coefficient in coefficients:
        total = field.add(total, field.mul(coefficient, power))
        power = field.pow(power, field.characteristic)
    return total


def _check_basis(field, basis):
    rows = {}
    for value in basis:
        vector = list(field.coefficients(value))
        for pivot in sorted(rows):
            factor = vector[pivot]
            if factor:
                vector = [(a - factor*b) % field.characteristic for a, b in zip(vector, rows[pivot])]
        pivot = next((i for i, value in enumerate(vector) if value), None)
        if pivot is None:
            raise ValueError('domain basis must be independent over the prime subfield')
        inverse = pow(vector[pivot], -1, field.characteristic)
        rows[pivot] = tuple(value*inverse % field.characteristic for value in vector)


def prepare_model(field, curve, fixed_point, *, offset=0, basis=(), interval_length=None, interval_step=1,
                  max_interval=256):
    """Prepare compact affine-subspace or interval domain and nonsingular curve.

    ``fixed_point`` is a finite rational point where F must vanish. A later
    source-point reconstruction sums to its negative. ``basis`` specifies an
    Fp-subspace D=offset+V; alternatively ``interval_length`` specifies
    D=offset+{0,...,length-1} embedded in the prime subfield without duplicates.
    Interval polynomial construction is charged preprocessing, not a free step.
    """
    if field.characteristic == 2:
        raise ValueError('this component requires odd characteristic')
    validated = GeneralS3Fiber(field, curve)
    curve = dict(validated.curve)
    offset = field.validate(offset)
    interval_step = field.validate(interval_step)
    fixed_point = tuple(fixed_point)
    if len(fixed_point) != 2:
        raise ValueError('fixed_point must be a coordinate pair')
    z, y = (field.validate(value) for value in fixed_point)
    cubic = (curve['a6'], curve['a4'], curve['a2'], 1)
    if field.square(y) != poly.evaluate_poly(cubic, z, field):
        raise ValueError('fixed_point must lie on the completed curve')
    basis = tuple(field.validate(value) for value in basis)
    _count(max_interval, 'max_interval')
    if interval_length is not None:
        _count(interval_length, 'interval_length')
        if basis:
            raise ValueError('choose an affine subspace or an interval, not both')
        if not interval_step:
            raise ValueError('interval_step must be nonzero')
        if not 1 <= interval_length <= min(field.characteristic, max_interval):
            raise ValueError('interval length must be positive, fit Fp and respect max_interval')
        domain_polynomial = (1,)
        for integer in range(interval_length):
            root = field.add(offset, field.mul(interval_step, integer))
            domain_polynomial = poly.mul(domain_polynomial, (field.neg(root), 1), field)
        descriptor = {'kind': 'interval', 'size': interval_length,
                      'sum_size': min(field.characteristic, 3*interval_length-2),
                      'polynomial': domain_polynomial, 'step': interval_step,
                      'inverse_step': field.inv(interval_step)}
    else:
        _check_basis(field, basis)
        coefficients, scalars = (1,), []
        for vector in basis:
            delta = _linearized_value(coefficients, vector, field)
            if not delta:
                raise ArithmeticError('independent-domain annihilator construction failed')
            multiplier = field.pow(delta, field.characteristic-1)
            shifted = (0,) + tuple(field.pow(coefficient, field.characteristic)
                                    for coefficient in coefficients)
            coefficients = tuple(field.sub(shifted[i], field.mul(multiplier, coefficients[i])
                                               if i < len(coefficients) else 0)
                                 for i in range(len(shifted)))
            scalars.append(multiplier)
        descriptor = {'kind': 'affine_subspace', 'size': field.characteristic ** len(basis),
                      'sum_size': field.characteristic ** len(basis),
                      'linearized': coefficients, 'prefix_scalars': tuple(scalars),
                      'domain_constant': _linearized_value(coefficients, offset, field)}
    return {'field': field, 'curve': curve, 'fixed_point': (z, y), 'offset': offset,
            'basis': basis, 'domain': descriptor}


def domain_contains(model, value):
    f, domain = model['field'], model['domain']
    shifted = f.sub(f.validate(value), model['offset'])
    if domain['kind'] == 'interval':
        return f.mul(shifted, domain['inverse_step']) < domain['size']
    return _linearized_value(domain['linearized'], shifted, f) == 0


def source_sums(model):
    """Enumerate the exact sum cover, whose size must be admitted externally."""
    f, domain = model['field'], model['domain']
    offset = f.mul(_scalar(f, 3), model['offset'])
    if domain['kind'] == 'interval':
        for integer in range(domain['sum_size']):
            yield f.add(offset, f.mul(domain['step'], integer))
    else:
        for coordinates in product(range(f.characteristic), repeat=len(model['basis'])):
            value = offset
            for coefficient, basis in zip(coordinates, model['basis']):
                value = f.add(value, f.mul(coefficient, basis))
            yield value


def _sum_admitted(model, s):
    f, domain = model['field'], model['domain']
    shifted = f.sub(s, f.mul(_scalar(f, 3), model['offset']))
    if domain['kind'] == 'interval':
        return f.mul(shifted, domain['inverse_step']) < domain['sum_size']
    return _linearized_value(domain['linearized'], shifted, f) == 0


class _Ring:
    """Parameter-polynomial operations with bounded intermediate degrees."""
    def __init__(self, field, max_degree, max_products):
        self.field, self.max_degree, self.max_products = field, max_degree, max_products
        self.products = 0

    def check(self, value):
        if len(value)-1 > self.max_degree:
            raise _Limit('parameter intermediate degree exceeds configured guard')
        return value

    def add(self, a, b): return self.check(poly.add(a, b, self.field))
    def sub(self, a, b): return self.check(poly.sub(a, b, self.field))

    def mul(self, a, b):
        products = sum(bool(value) for value in a)*sum(bool(value) for value in b)
        if self.products + products > self.max_products:
            raise _Limit('quotient-construction scalar-product budget exhausted')
        self.products += products
        if a and b and len(a)+len(b)-2 > self.max_degree:
            raise _Limit('parameter intermediate degree exceeds configured guard')
        return self.check(poly.mul(a, b, self.field))

    def scale(self, a, coefficient): return self.mul(a, (coefficient,))

    def frobenius(self, a):
        if not a:
            return ()
        degree = (len(a)-1)*self.field.characteristic
        if degree > self.max_degree:
            raise _Limit('parameter Frobenius intermediate degree exceeds configured guard')
        result = [0]*(degree+1)
        for index, value in enumerate(a):
            if value:
                result[index*self.field.characteristic] = self.field.pow(value, self.field.characteristic)
        return poly.trim(result)


def coefficient_maps(model, source_sum):
    """Return c(B), d(B), E2(B), E3(B), each as ordinary parameter polynomials."""
    f, s = model['field'], model['field'].validate(source_sum)
    z, y = model['fixed_point']
    a2, a4 = model['curve']['a2'], model['curve']['a4']
    inverse_two = f.inv(_scalar(f, 2))
    c = (f.neg(f.mul(inverse_two, f.add(z, s))), 0, inverse_two)
    d = poly.sub((f.neg(f.square(z)), f.neg(y)), tuple(f.mul(z, value) for value in c), f)
    c_square = poly.mul(c, c, f)
    e2 = poly.sub(c_square, (f.square(z), f.mul(_scalar(f, 2), y), f.add(a2, z)), f)
    r = poly.add(c, (z,), f)
    e3 = tuple(f.mul(z, value) for value in poly.mul(r, r, f))
    coefficient = f.add(a4, f.add(f.mul(a2, z), f.square(z)))
    e3 = poly.add(e3, (0, 0, coefficient), f)
    e3 = poly.add(e3, poly.mul((0, f.mul(_scalar(f, 2), y)), r, f), f)
    return {'c': c, 'd': d, 'e2': e2, 'e3': e3}


def discriminant(model, s, maps=None):
    """Cubic discriminant in B; exact degree twelve, leading coefficient -1/16."""
    f = model['field']
    maps = maps or coefficient_maps(model, s)
    t, u = maps['e2'], maps['e3']
    result = tuple(f.mul(f.square(s), value) for value in poly.mul(t, t, f))
    for coefficient, term in [(_scalar(f, 4), poly.mul(poly.mul(t, t, f), t, f)),
                              (f.mul(_scalar(f, 4), f.pow(s, 3)), u),
                              (_scalar(f, 27), poly.mul(u, u, f))]:
        result = poly.sub(result, tuple(f.mul(coefficient, value) for value in term), f)
    result = poly.add(result, tuple(f.mul(f.mul(_scalar(f, 18), s), value)
                                    for value in poly.mul(t, u, f)), f)
    if len(result)-1 != 12 or result[-1] != f.neg(f.inv(_scalar(f, 16))):
        raise ArithmeticError('cubic discriminant degree/leading-coefficient identity failed')
    return result


def _residue_mul(left, right, s, maps, ring):
    values = [()] * 5
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if a and b:
                values[i+j] = ring.add(values[i+j], ring.mul(a, b))
    for degree in (4, 3):
        value = values[degree]
        if value:
            values[degree-1] = ring.add(values[degree-1], ring.scale(value, s))
            values[degree-2] = ring.sub(values[degree-2], ring.mul(value, maps['e2']))
            values[degree-3] = ring.add(values[degree-3], ring.mul(value, maps['e3']))
    return tuple(values[:3])


def _residue_power_x(exponent, s, maps, ring):
    value, base = ((1,), (), ()), ((), (1,), ())
    while exponent:
        if exponent & 1:
            value = _residue_mul(value, base, s, maps, ring)
        exponent >>= 1
        if exponent:
            base = _residue_mul(base, base, s, maps, ring)
    return value


def parameter_domain_remainder(model, source_sum, *, max_parameter_degree=512,
                               max_scalar_products=5_000_000):
    """Return three polynomial-in-B coefficients of P_D mod H_s exactly."""
    _count(max_parameter_degree, 'max_parameter_degree')
    _count(max_scalar_products, 'max_scalar_products')
    f, domain, s = model['field'], model['domain'], f_validate(model, source_sum)
    k = domain['size']
    if 2*k > max_parameter_degree:
        raise _Limit(f'certified 2k parameter-degree bound {2*k} exceeds max_parameter_degree')
    # Coefficient products may temporarily exceed the final 2k degree bound
    # before exact cancellations; keep this extra degree guard explicit.
    extra = 4*f.characteristic if domain['kind'] == 'affine_subspace' and model['basis'] else 0
    ring = _Ring(f, max_parameter_degree + extra + 16, max_scalar_products)
    maps = coefficient_maps(model, s)
    if domain['kind'] == 'interval':
        residue = ((), (), ())
        x = ((), (1,), ())
        for coefficient in reversed(domain['polynomial']):
            residue = _residue_mul(residue, x, s, maps, ring)
            residue = (ring.add(residue[0], (coefficient,)), residue[1], residue[2])
    else:
        residue = ((), (1,), ())
        if model['basis']:
            matrix = (((1,), (), ()),
                      _residue_power_x(f.characteristic, s, maps, ring),
                      _residue_power_x(2*f.characteristic, s, maps, ring))
            for multiplier in domain['prefix_scalars']:
                next_residue = [(), (), ()]
                for column, coefficient in zip(matrix, residue):
                    coefficient = ring.frobenius(coefficient)
                    for index, value in enumerate(column):
                        if coefficient and value:
                            next_residue[index] = ring.add(next_residue[index], ring.mul(coefficient, value))
                residue = tuple(ring.sub(value, ring.scale(previous, multiplier))
                                for value, previous in zip(next_residue, residue))
        residue = (ring.sub(residue[0], (domain['domain_constant'],)), residue[1], residue[2])
    if any(len(value)-1 > 2*k for value in residue):
        raise ArithmeticError('proved final parameter-degree bound was violated')
    return {'coefficients': residue, 'coefficient_maps': maps,
            'stats': {'domain_degree': k, 'certified_parameter_degree_bound': 2*k,
                      'actual_parameter_degrees': [len(value)-1 for value in residue],
                      'quotient_scalar_products': ring.products,
                      'intermediate_degree_guard': ring.max_degree,
                      'domain_materialized': False, 'field_materialized': False,
                      'interval_roots_enumerated_in_preprocessing': domain['kind'] == 'interval'}}


def f_validate(model, value):
    return model['field'].validate(value)


def cubic_from_parameter(model, s, parameter):
    f = model['field']
    parameter = f.validate(parameter)
    maps = coefficient_maps(model, s)
    t, u = (poly.evaluate_poly(maps[name], parameter, f) for name in ('e2', 'e3'))
    return (f.neg(u), t, f.neg(s), 1)


def numeric_domain_remainder(model, cubic):
    f, domain = model['field'], model['domain']
    if domain['kind'] == 'interval':
        return poly.divmod_poly(domain['polynomial'], cubic, f)[1]
    residue = (0, 1)
    for multiplier in domain['prefix_scalars']:
        residue = poly.sub(poly.powmod_poly(residue, f.characteristic, cubic, f),
                           tuple(f.mul(multiplier, value) for value in residue), f)
    residue = poly.sub(residue, (domain['domain_constant'],), f)
    return poly.divmod_poly(residue, cubic, f)[1]


def supported_cubic(model, cubic):
    """Exact repeated-root support test H | P_D³; no domain enumeration."""
    remainder = numeric_domain_remainder(model, cubic)
    return not poly.powmod_poly(remainder, 3, cubic, model['field'])


def conditional_family(model, source_sum, *, max_parameter_degree=512,
                       max_scalar_products=5_000_000, max_split_trials=256,
                       max_parameters=10_000, seed=0):
    """Every parameter B whose norm cubic splits completely in D, with repeats.

    Exact compact remainder and discriminant roots cover disjoint simple and
    repeated support branches. A splitting or construction budget exhaustion
    returns unresolved, never a complete empty family. This family itself is
    an exact intermediate component; it has no rational-point output claim.
    """
    f, s = model['field'], f_validate(model, source_sum)
    for name, value in [('max_split_trials', max_split_trials), ('max_parameters', max_parameters),
                        ('seed', seed)]:
        _count(value, name)
    result = {'schema': 'odd_curve_exact_domain_family_v1', 'status': 'unresolved',
              'scope': 'all norm parameters whose monic cubic has all roots in the compact source domain',
              'source_sum': s, 'parameters': [], 'solution_count_exact': None,
              'search_complete': False, 'output_complete': False, 'verified': True,
              'stats': {'split_trials': 0}, 'reason': None}
    if not _sum_admitted(model, s):
        result.update(status='complete', solution_count_exact=0, search_complete=True,
                      output_complete=True, reason='source sum is outside the exact sum cover')
        return result
    try:
        construction = parameter_domain_remainder(model, s, max_parameter_degree=max_parameter_degree,
                                                 max_scalar_products=max_scalar_products)
        result['stats'].update(construction['stats'])
        common = ()
        for coefficient in construction['coefficients']:
            common = poly.gcd_poly(common, coefficient, f)
            if common == (1,):
                break
        if not common:
            raise ArithmeticError('supported-cubic parameter constraints cannot all vanish identically')
        delta = discriminant(model, s, construction['coefficient_maps'])
        result['stats']['simple_parameter_polynomial_degree'] = len(common)-1
        result['stats']['repeated_parameter_polynomial_degree'] = len(delta)-1
        admitted = set()
        for branch, polynomial in [('distinct', common), ('repeated', delta)]:
            report = poly.field_roots(polynomial, f, seed=seed,
                                      max_trials=max_split_trials-result['stats']['split_trials'])
            result['stats']['split_trials'] += report['split_trials']
            if report['status'] != 'complete':
                raise _Limit(f'{branch} parameter root splitting unresolved; exact degree '
                             f"{report['restricted_degree']} but {len(report['roots'])} roots isolated")
            for parameter in report['roots']:
                cubic = cubic_from_parameter(model, s, parameter)
                if branch == 'distinct':
                    if numeric_domain_remainder(model, cubic):
                        raise ArithmeticError('parameter remainder root failed exact support replay')
                    if not poly.evaluate_poly(delta, parameter, f):
                        raise ArithmeticError('distinct support branch contained a repeated cubic')
                    admitted.add(parameter)
                elif supported_cubic(model, cubic):
                    admitted.add(parameter)
        values = sorted(admitted)
        result.update(status='complete' if len(values) <= max_parameters else 'output_truncated',
                      parameters=values[:max_parameters], solution_count_exact=len(values),
                      search_complete=True, output_complete=len(values) <= max_parameters,
                      reason='Exact compact-domain remainder and every degree-twelve repeated fiber exhausted')
    except _Limit as error:
        result['reason'] = str(error)
    return result
