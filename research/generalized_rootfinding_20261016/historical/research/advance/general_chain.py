"""Exact five-source S3 separator join over arbitrary finite fields.

For odd characteristic the curve is the completed form
Y^2 = C(X) = X^3 + a2*X^2 + a4*X + a6, with gcd(C,C')=1.
For characteristic two this uses the ordinary binary curve S3 with b != 0.
For odd characteristic a zero fiber is impossible: if its leading coefficient
(x-y)^2 vanishes, x=y=t; its linear coefficient is -4*C(t). If that
also vanishes, its constant coefficient is C'(t)^2, nonzero because C is
squarefree. The binary b != 0 case likewise has no unrestricted fiber.
The output is every finite-coordinate polynomial-chain root, not a claim of
rational point lifting, a curve-group decomposition, novelty or timing gain.

A full-field intermediate separator joins three source coordinates on the
left with two on the right. Exact quadratic fibers replace enumeration of
three intermediate field variables. Runtime is O(k^3 + k^2) root calls plus
root extraction and explicit output, and storage is at most 4*k^2 records.
"""

import math
import time

from rootfinder import polynomial as poly
from rootfinder.binaryfield import BinaryField
from research.advance.finite_fiber import BinaryQuadraticFiber


class _Unresolved(Exception):
    pass


def _budget(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


class GeneralS3Fiber:
    """Verified nonsingular curve S3 and exact degree-at-most-two fibers.

    Generic finite fields use the existing certified finite-field root
    isolator. Its exhausted splitting budget is an unresolved result, never
    an empty root set. BinaryField uses deterministic cached Artin–Schreier
    inversion, including even extension degrees.
    """

    def __init__(self, field, curve=None, *, max_split_trials=128, seed=0):
        for name in ('characteristic', 'degree', 'order', 'validate', 'add',
                     'sub', 'neg', 'mul', 'pow', 'inv', 'square'):
            if not hasattr(field, name):
                raise TypeError(f"finite field backend is missing {name}")
        self.field = field
        self.max_split_trials = _budget(max_split_trials, 'max_split_trials')
        self.seed = _budget(seed, 'seed')
        self.split_trials = 0
        self._binary = None
        if curve is None:
            curve = {'b': 1} if field.characteristic == 2 else {'a2': 0, 'a4': 0, 'a6': 1}
        if not isinstance(curve, dict):
            raise TypeError('curve must be a dictionary of canonical field coefficients')
        if field.characteristic == 2:
            if set(curve) - {'b', 'a2'}:
                raise ValueError('binary curve parameters are b and optional a2')
            self.curve = {'b': field.validate(curve.get('b', 1)),
                          'a2': field.validate(curve.get('a2', 0))}
            if self.curve['b'] == 0:
                raise ValueError('binary ordinary curve requires b != 0 for bounded fibers')
        else:
            if set(curve) - {'a2', 'a4', 'a6'}:
                raise ValueError('completed odd-characteristic curve parameters are a2, a4, a6')
            self.curve = {name: field.validate(curve.get(name, default))
                          for name, default in (('a2', 0), ('a4', 0), ('a6', 1))}
            a2, a4, a6 = (self.curve[name] for name in ('a2', 'a4', 'a6'))
            cubic = (a6, a4, a2, 1)
            derivative = (a4, field.mul(2 % field.characteristic, a2),
                          3 % field.characteristic)
            if poly.gcd_poly(cubic, derivative, field) != (1,):
                raise ValueError("curve C must be squarefree: gcd(C,C') != 1")

    def coefficients(self, x, y):
        """Return ascending coefficients (c,b,a) of S3(x,y,Z)."""
        f = self.field
        x, y = f.validate(x), f.validate(y)
        xy = f.mul(x, y)
        if f.characteristic == 2:
            return (f.add(f.square(xy), self.curve['b']), xy,
                    f.square(f.add(x, y)))
        a2, a4, a6 = (self.curve[name] for name in ('a2', 'a4', 'a6'))
        two, four = 2 % f.characteristic, 4 % f.characteristic
        s = f.add(x, y)
        a = f.square(f.sub(x, y))
        inside = f.add(f.mul(s, f.add(xy, a4)),
                       f.add(f.mul(f.mul(two, a2), xy), f.mul(two, a6)))
        b = f.neg(f.mul(two, inside))
        c = f.sub(f.square(f.sub(a4, xy)),
                  f.mul(f.mul(four, a6), f.add(s, a2)))
        return c, b, a

    def evaluate(self, x, y, z):
        """Evaluate symmetric S3 directly, independently of the fiber formula."""
        f = self.field
        x, y, z = (f.validate(value) for value in (x, y, z))
        xy = f.mul(x, y)
        if f.characteristic == 2:
            return f.add(f.add(f.mul(f.square(f.add(x, y)), f.square(z)),
                               f.mul(xy, z)),
                         f.add(f.square(xy), self.curve['b']))
        a2, a4, a6 = (self.curve[name] for name in ('a2', 'a4', 'a6'))
        pair_sum = f.add(xy, f.add(f.mul(x, z), f.mul(y, z)))
        total = f.add(f.add(f.add(x, y), z), a2)
        product_plus = f.add(f.mul(xy, z), a6)
        return f.sub(f.square(f.sub(a4, pair_sum)),
                     f.mul(4 % f.characteristic, f.mul(total, product_plus)))

    def roots(self, x, y):
        """All distinct exact field roots, or _Unresolved for exhausted splitting."""
        f = self.field
        x, y = f.validate(x), f.validate(y)
        c, b, a = self.coefficients(x, y)
        if a == b == c == 0:
            raise ArithmeticError('nonsingular S3 produced an unrestricted quadratic fiber')
        if isinstance(f, BinaryField):
            if self._binary is None:
                self._binary = BinaryQuadraticFiber(f)
            roots = self._binary.roots(a, b, c)
        else:
            report = poly.field_roots((c, b, a), f,
                                      seed=self.seed ^ (x * 0x9E3779B1 + y),
                                      max_trials=self.max_split_trials)
            self.split_trials += report['split_trials']
            if report['status'] != 'complete':
                raise _Unresolved('exact quadratic splitting budget exhausted; '
                                  f"{len(report['roots'])} of {report['restricted_degree']} roots isolated")
            roots = tuple(report['roots'])
            if len(roots) != report['restricted_degree']:
                raise ArithmeticError('complete quadratic report has inconsistent root count')
        if len(roots) > 2 or len(set(roots)) != len(roots):
            raise ArithmeticError('quadratic distinct-root bound or uniqueness failed')
        if any(self.evaluate(x, y, root) != 0 for root in roots):
            raise ArithmeticError('quadratic root failed direct symmetric S3 verification')
        return tuple(sorted(roots))


def solve_general_chain(field, domain, target, *, curve=None,
                        max_domain=4096, max_source_assignments=1_000_000,
                        max_root_calls=5_000_000, max_records=1_000_000,
                        max_output=10_000, max_split_trials=128, seed=0,
                        max_seconds=None):
    """Find/count all ordered (x1,x2,x3,x4,x5,u1,u2,u3) chain roots.

    Sources lie in the supplied explicit common domain; intermediates range
    over the entire field. All four original symmetric S3 links are enforced.
    Output caps preserve full bucket multiplicities and the exact count.
    Preflight caps include worst-case quadratic-call and right-record bounds.
    """
    started = time.perf_counter()
    budgets = {name: _budget(value, name) for name, value in (
        ('domain_size', max_domain), ('source_assignments', max_source_assignments),
        ('root_calls', max_root_calls), ('right_records', max_records))}
    max_output = _budget(max_output, 'max_output')
    if max_seconds is not None and (isinstance(max_seconds, bool)
            or not isinstance(max_seconds, (int, float))
            or not math.isfinite(max_seconds) or max_seconds < 0):
        raise ValueError('max_seconds must be finite and nonnegative or None')
    fiber = GeneralS3Fiber(field, curve, max_split_trials=max_split_trials, seed=seed)
    target = field.validate(target)
    domain = tuple(domain)
    for value in domain:
        field.validate(value)
    if len(set(domain)) != len(domain):
        raise ValueError('domain must contain distinct canonical field elements')
    domain = tuple(sorted(domain))
    k = len(domain)
    planned = {'domain_size': k, 'source_assignments': k + 2 * k * k + k ** 3,
               'root_calls': k + 3 * k * k + 2 * k ** 3, 'right_records': 4 * k * k}
    stats = {'planned': planned, 'budgets': budgets, 'source_assignments': 0,
             'root_calls': 0, 'right_records': 0, 'right_keys': 0,
             'left_partial_records': 0, 'joined_assignments': 0,
             'split_trials': 0, 'output_cap': max_output}
    result = {'schema': 'general-finite-field-s3-chain-result-v1',
              'algorithm': 'exact_general_s3_separator_fiber_join',
              'scope': 'all finite-coordinate polynomial-chain roots; no curve-group decomposition claim',
              'field': {'characteristic': field.characteristic, 'degree': field.degree,
                        'order': field.order, 'modulus': getattr(field, 'modulus', None)},
              'curve': fiber.curve, 'target': target, 'domain': domain,
              'status': 'unresolved', 'search_complete': False, 'output_complete': False,
              'solution_count_exact': None, 'solution_count_lower_bound': 0,
              'solutions': [], 'stats': stats, 'unresolved_reason': None,
              'verified': True,
              'certificate': 'global exhaustion not established'}

    def finish():
        stats['split_trials'] = fiber.split_trials
        stats['pipeline_seconds'] = time.perf_counter() - started
        result['solutions'].sort()
        return result

    exceeded = [name for name, cap in budgets.items() if planned[name] > cap]
    if exceeded:
        result['unresolved_reason'] = 'preflight bound exceeds ' + ', '.join(exceeded) + ' budget'
        return finish()

    def check_time():
        if max_seconds is not None and time.perf_counter() - started >= max_seconds:
            raise _Unresolved('charged pipeline time budget exhausted')

    def roots(left, right):
        check_time()
        if stats['root_calls'] >= budgets['root_calls']:
            raise ArithmeticError('preflight quadratic-call bound was violated')
        stats['root_calls'] += 1
        return fiber.roots(left, right)

    table = {}
    count = 0
    try:
        check_time()
        for x5 in domain:
            stats['source_assignments'] += 1
            u3_values = roots(x5, target)
            for x4 in domain:
                stats['source_assignments'] += 1
                for u3 in u3_values:
                    for u2 in roots(u3, x4):
                        table.setdefault(u2, []).append((x4, x5, u3))
                        stats['right_records'] += 1
                        if stats['right_records'] > budgets['right_records']:
                            raise ArithmeticError('preflight right-record bound was violated')
        stats['right_keys'] = len(table)
        for x1 in domain:
            for x2 in domain:
                stats['source_assignments'] += 1
                u1_values = roots(x1, x2)
                for x3 in domain:
                    stats['source_assignments'] += 1
                    for u1 in u1_values:
                        for u2 in roots(u1, x3):
                            stats['left_partial_records'] += 1
                            matches = table.get(u2, ())
                            # Distinct source/intermediate witnesses under one
                            # separator key are retained with full multiplicity.
                            count += len(matches)
                            capacity = max_output - len(result['solutions'])
                            for x4, x5, u3 in matches[:capacity]:
                                check_time()
                                links = ((x1, x2, u1), (u1, x3, u2),
                                         (u2, x4, u3), (u3, x5, target))
                                if any(fiber.evaluate(*link) != 0 for link in links):
                                    raise ArithmeticError('retained join tuple failed original S3 verification')
                                result['solutions'].append([x1, x2, x3, x4, x5, u1, u2, u3])
        check_time()
        result.update(search_complete=True, solution_count_exact=count,
                      output_complete=len(result['solutions']) == count,
                      certificate='all exact quadratic fibers on both sides of full-field u2 exhausted')
        result['status'] = 'complete' if result['output_complete'] else 'output_truncated'
    except _Unresolved as error:
        result['unresolved_reason'] = str(error)
    result['solution_count_lower_bound'] = count
    stats['joined_assignments'] = count
    stats['right_keys'] = len(table)
    return finish()
