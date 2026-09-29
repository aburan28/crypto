"""Small-field, phase-matched Semaev encoding comparison. No challenge inputs."""
import collections
import heapq
import itertools
import json
import math
from pathlib import Path
import random
import resource
import sys
import time

sys.path.insert(0, str(Path(__file__).parent / 'vendor'))
from frobenius_cubic_quotient import CyclotomicFactorSpace
from frobenius_curve_arithmetic import Field, Curve
from test_cubic_quotient import NormalCoordinates

MODULI = {5: 0x25, 13: 0x201b, 19: 0x80027}


class Budget(Exception):
    pass


def check(deadline):
    if time.perf_counter() >= deadline:
        raise Budget('time budget')


def addpoly(*items):
    result = {}
    for item in items:
        for m, c in item.items():
            value = result.get(m, 0) ^ c
            if value:
                result[m] = value
            elif m in result:
                del result[m]
    return result


class Algebra:
    """Field-valued Boolean ANFs; x_i^2=x_i, no truth-table construction."""
    def __init__(self, field, counts, deadline=math.inf):
        self.f, self.counts, self.deadline = field, counts, deadline

    def mul(self, a, b):
        result = {}
        for i, (u, c) in enumerate(a.items()):
            if i % 8 == 0:
                check(self.deadline)
            for v, d in b.items():
                m = u | v
                result[m] = result.get(m, 0) ^ self.f.mul(c, d)
        self.counts['anf_term_products'] += len(a) * len(b)
        return {m: c for m, c in result.items() if c}

    def square(self, a):
        # Frobenius is linear over the Boolean coefficients.
        self.counts['anf_coefficient_squares'] += len(a)
        return {m: self.f.sq(c) for m, c in a.items()}

    def power(self, a, exponent):
        result = {0: 1}
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            exponent >>= 1
            if exponent:
                a = self.square(a)
        return result

    def s4(self, x, y, z, t):
        # Res_z (a*z^2+b*z+c, d*z^2+e*z+f), characteristic two.
        # Each f3(a,b,c)=(ab+ac+bc)^2+abc+1 for E_0.
        b = self.mul(x, y)
        a = self.square(addpoly(x, y))
        c = addpoly(self.square(b), {0: 1})
        d = self.square(addpoly(z, {0: t}))
        e = self.mul(z, {0: t})
        f = addpoly(self.square(e), {0: 1})
        afcd = addpoly(self.mul(a, f), self.mul(c, d))
        aebd = addpoly(self.mul(a, e), self.mul(b, d))
        bfce = addpoly(self.mul(b, f), self.mul(c, e))
        return addpoly(self.square(afcd), self.mul(aebd, bfce))

    def s3(self, x, y, z):
        xy = self.mul(x, y)
        return addpoly(self.square(addpoly(xy, self.mul(x,z), self.mul(y,z))),
                       self.mul(xy,z), {0:1})


def evaluate(poly, assignment):
    value = 0
    for m, c in poly.items():
        if assignment & m == m:
            value ^= c
    return value


def normal_polynomial(poly, space, normal, phase):
    return {m: normal.to_field(space.q.rotate(space.q.unproject(0, c), phase))
            for m, c in poly.items()}


def map_polynomial(space, normal, phase, offset, arm, counts):
    s = space.s
    if arm == 'linear_payload':
        return {1 << (offset + i): normal.to_field(space.q.rotate(word, phase))
                for i, word in enumerate(space.basis)}
    # u ranges LINEARLY over W in the same W basis used for z payloads.
    # z=u^(p^-1 mod (2^s-1)); this is the actual quotient inverse on W.
    u = {1 << (offset + i): space.q.project(word)[1] for i, word in enumerate(space.basis)}
    inverse = pow(space.p, -1, (1 << s) - 1)
    poly = Algebra(space.q.field, counts).power(u, inverse)
    return normal_polynomial(poly, space, normal, phase)


def group_sum(curve, points):
    total = None
    for point in points:
        total = curve.add(total, point)
    return total


def proper(curve, points):
    return all(p != curve.neg(q) for i, p in enumerate(points) for q in points[:i])


def lift(curve, xs, target):
    choices = []
    for x in xs:
        if not x:
            return None
        points = [p for p in curve.lift(x) if curve.scale(p, curve.r) is None]
        if not points:
            return None
        choices.append(points)
    valid = [tuple(points) for points in itertools.product(*choices)
             if proper(curve, points) and group_sum(curve, points) == target]
    return min(valid) if valid else None


def rank_add(matrix, original, prime):
    row = list(original)
    for column in range(len(row)):
        if not row[column]:
            continue
        if column in matrix:
            coeff = row[column]
            row = [(a - coeff * b) % prime for a, b in zip(row, matrix[column])]
        else:
            coeff = pow(row[column], -1, prime)
            matrix[column] = [(a * coeff) % prime for a in row]
            return True
    return False


class Setup:
    def __init__(self, n, s, seed, arm):
        self.n, self.s, self.seed, self.arm = n, s, seed, arm
        self.counts = collections.Counter()
        self.curve = Curve(Field(n, MODULI[n]))
        self.normal = NormalCoordinates(self.curve.f, seed)
        self.space = CyclotomicFactorSpace(n, s)
        self.phases = (0, 1, 2)
        self.maps = [map_polynomial(self.space, self.normal, phase, i*s, arm, self.counts)
                     for i, phase in enumerate(self.phases)]
        # Common factor-base bookkeeping, charged in BOTH arms. This is not
        # used to constrain payloads or give an existence oracle to the solver.
        self.reps, self.lookup, self.domains = [], {}, []
        seed_points = []
        for payload in range(1, 1 << s):
            x = self.normal.to_field(self.space.encode(payload))
            points = [p for p in self.curve.lift(x) if self.curve.scale(p, self.curve.r) is None]
            if not points:
                continue
            rep = min(points)
            column = len(self.reps)
            self.reps.append(rep)
            seed_points.extend(points)
            point, coeff = rep, 1
            for phase in range(n):
                for signed, value in ((point, coeff), (self.curve.neg(point), -coeff % self.curve.r)):
                    assert signed not in self.lookup or self.lookup[signed] == (column, value)
                    self.lookup[signed] = column, value
                point = self.curve.frob(point)
                coeff = coeff * self.curve.lam % self.curve.r
        self.domains = [[self.curve.frob(p, phase) for p in seed_points] for phase in self.phases]
        self.matrix = {}

    def certificate(self):
        return {'degree': self.n, 'payload_bits': self.s, 'seed': self.seed,
                'modulus': MODULI[self.n], 'normal_generator': self.normal.beta,
                'normal_basis': self.normal.basis, 'generator': self.curve.g,
                'subgroup_order': self.curve.r, 'frobenius_eigenvalue': self.curve.lam,
                'phases': self.phases, 'representatives': self.reps,
                'signed_base_points': len(self.lookup),
                'phase_domain_sizes': [len(x) for x in self.domains],
                'inverse_exponent': pow(self.n, -1, (1 << self.s)-1),
                'coordinate_map_degree': max((m.bit_count() for p in self.maps for m in p), default=0),
                'coordinate_map_terms': [len(p) for p in self.maps]}

    def row(self, witness):
        result = [0] * len(self.reps)
        for point in witness:
            column, coeff = self.lookup[tuple(point)]
            result[column] = (result[column] + coeff) % self.curve.r
        assert group_sum(self.curve, [self.curve.scale(p, a) for p, a in zip(self.reps, result)]) == group_sum(self.curve, witness)
        return result

    def inputs(self, count=6, planted=2):
        rng = random.Random(self.seed * 1000 + self.n * 10 + self.s)
        records = []
        # Uniform nonzero scalars; sampler has no access to coverage labels.
        for i in range(count):
            scalar = rng.randrange(1, self.curve.r)
            records.append({'kind': 'natural', 'index': i, 'scalar': scalar,
                            'target': self.curve.scale(self.curve.g, scalar)})
        if all(self.domains):
            for i in range(planted):
                for _ in range(1000):
                    points = [rng.choice(d) for d in self.domains]
                    target = group_sum(self.curve, points)
                    if target is not None and proper(self.curve, points):
                        records.append({'kind': 'planted_control', 'index': i, 'target': target})
                        break
        return records


def boolean_equations(poly, n):
    return [{m for m, c in poly.items() if c >> bit & 1} for bit in range(n)]


def build_system(setup, target, deadline=math.inf):
    algebra=Algebra(setup.curve.f,setup.counts,deadline)
    intermediate={1 << (3*setup.s+i):1<<i for i in range(setup.n)}
    x,y,z=setup.maps
    polys=[algebra.s3(x,y,intermediate),algebra.s3(z,{0:target[0]},intermediate)]
    equations=[p for poly in polys for p in boolean_equations(poly,setup.n)]
    return polys,equations


def leading(poly):
    return max(poly, key=lambda m: (m.bit_count(), m))


def times_monomial(poly, multiplier):
    result = set()
    for m in poly:
        t = m | multiplier
        if t in result:
            result.remove(t)
        else:
            result.add(t)
    return result


def boolean_groebner_reference(equations, seconds=1):
    """Deterministic Boolean Buchberger, including x_i^2+x_i S-pairs.

    This is a bounded algebraic diagnostic, NOT a fast F4 implementation.
    Field pairs are (x_i+1)*g when x_i divides LM(g). Pair products are
    immediately reduced by x_i^2=x_i. Ordering is total-degree then bitmask.
    """
    deadline = time.perf_counter() + seconds
    basis, lms, queue = [], [], []
    peak_degree = max([2] + [m.bit_count() for p in equations for m in p])
    counts = collections.Counter()
    def reduce(poly):
        p, rem = set(poly), set()
        while p:
            check(deadline)
            lm = leading(p)
            for g, lg in zip(basis, lms):
                if lm & lg == lg:
                    p.symmetric_difference_update(times_monomial(g, lm ^ lg))
                    counts['reductions'] += 1
                    break
            else:
                rem.add(lm)
                p.remove(lm)
        return rem
    def insert(p):
        lm = leading(p)
        i = len(basis)
        for j, other in enumerate(lms):
            common = lm | other
            heapq.heappush(queue, (common.bit_count(), 0, i, j))
        for bit in range(lm.bit_length()):
            if lm >> bit & 1:
                heapq.heappush(queue, (lm.bit_count()+1, 1, i, bit))
        basis.append(p)
        lms.append(lm)
    status = 'complete'
    start = time.perf_counter()
    try:
        for p in equations:
            reduced = reduce(p)
            if reduced:
                insert(reduced)
            if reduced == {0}:
                queue.clear()
                break
        while queue:
            degree, kind, i, j = heapq.heappop(queue)
            peak_degree = max(peak_degree, degree)
            counts['critical_pairs'] += 1
            if kind:
                candidate = times_monomial(basis[i], 1 << j) ^ basis[i]
            else:
                common = lms[i] | lms[j]
                candidate = times_monomial(basis[i], common ^ lms[i]) ^ times_monomial(basis[j], common ^ lms[j])
            reduced = reduce(candidate)
            if reduced:
                insert(reduced)
            if reduced == {0}:
                queue.clear()
                break
    except Budget:
        status = 'timeout'
    return {'status': status, 'seconds': time.perf_counter()-start,
            'boolean_buchberger_solving_degree': peak_degree if status == 'complete' else None,
            'observed_peak_degree': peak_degree, 'basis_count': len(basis),
            'basis_max_degree': max([0]+[m.bit_count() for p in basis for m in p]),
            'remaining_pairs': len(queue), 'counts': dict(counts)}, basis


def boolean_groebner(equations, seconds=1):
    """Same Boolean Buchberger order as reference, with packed GF(2) rows."""
    from functools import lru_cache
    start = time.perf_counter()
    deadline = start + seconds
    nv = max([0] + [m.bit_length() for p in equations for m in p])
    if nv > 18:
        # Do not allocate an exponential packed monomial universe for large
        # chained systems. Both implementations have the same pair order.
        return boolean_groebner_reference(equations,seconds)
    ordered = sorted(range(1 << nv), key=lambda m: (m.bit_count(), m))
    positions = [0] * len(ordered)
    for i, m in enumerate(ordered):
        positions[m] = i
    def pack(p):
        value = 0
        for m in p:
            value ^= 1 << positions[m]
        return value
    def unpack(p):
        result = set()
        while p:
            bit = p & -p
            result.add(ordered[bit.bit_length()-1])
            p ^= bit
        return result
    basis, supports, lms, queue = [], [], [], []
    counts = collections.Counter()
    peak_degree = max([2] + [m.bit_count() for p in equations for m in p])
    @lru_cache(maxsize=1024)
    def multiple(i, multiplier):
        return pack(times_monomial(supports[i], multiplier))
    def reduce(p):
        rem = 0
        while p:
            check(deadline)
            bitpos = p.bit_length()-1
            lm = ordered[bitpos]
            for i, lg in enumerate(lms):
                if lm & lg == lg:
                    p ^= multiple(i, lm ^ lg)
                    counts['reductions'] += 1
                    break
            else:
                rem ^= 1 << bitpos
                p ^= 1 << bitpos
        return rem
    def insert(p):
        lm = ordered[p.bit_length()-1]
        i = len(basis)
        for j, other in enumerate(lms):
            heapq.heappush(queue, ((lm | other).bit_count(), 0, i, j))
        for bit in range(lm.bit_length()):
            if lm >> bit & 1:
                heapq.heappush(queue, (lm.bit_count()+1, 1, i, bit))
        basis.append(p); supports.append(unpack(p)); lms.append(lm)
    status = 'complete'
    try:
        for p in equations:
            reduced = reduce(pack(p))
            if reduced:
                insert(reduced)
            if reduced == 1:
                queue.clear(); break
        while queue:
            degree, kind, i, j = heapq.heappop(queue)
            peak_degree = max(peak_degree, degree)
            counts['critical_pairs'] += 1
            if kind:
                candidate = multiple(i, 1 << j) ^ basis[i]
            else:
                common = lms[i] | lms[j]
                candidate = multiple(i, common ^ lms[i]) ^ multiple(j, common ^ lms[j])
            reduced = reduce(candidate)
            if reduced:
                insert(reduced)
            if reduced == 1:
                queue.clear(); break
    except Budget:
        status = 'timeout'
    return {'status': status, 'seconds': time.perf_counter()-start,
            'boolean_buchberger_solving_degree': peak_degree if status == 'complete' else None,
            'observed_peak_degree': peak_degree, 'basis_count': len(basis),
            'basis_max_degree': max([0]+[m.bit_count() for p in supports for m in p]),
            'remaining_pairs': len(queue), 'counts': dict(counts)}, supports


def xor_tree(z3, values, ctx=None):
    if not values:
        return z3.BoolVal(False, ctx=ctx)
    values = list(values)
    while len(values) > 1:
        values = [z3.Xor(values[i], values[i+1]) if i+1<len(values) else values[i]
                  for i in range(0, len(values), 2)]
    return values[0]


def solve_query(setup, record, encoding_seconds=20, query_seconds=2, diagnostic=False, context_mode="shared", solver_seed=260926):
    import z3
    output = dict(record)
    counts_before = setup.counts.copy()
    start = time.perf_counter()
    output['phase_seconds'] = {}
    try:
        deadline = start + encoding_seconds
        polys,equations = build_system(setup,record['target'],deadline)
        monomial_support=set().union(*(p.keys() for p in polys))
        output['anf_degree'] = max([0] + [m.bit_count() for m in monomial_support])
        output['anf_field_terms'] = sum(map(len,polys))
        output['anf_boolean_terms'] = sum(map(len, equations))
        output['phase_seconds']['anf_expansion'] = time.perf_counter()-start
        ctx = z3.Context() if context_mode == 'fresh' else z3.main_ctx()
        variables = [z3.Bool('b%d' % i, ctx=ctx) for i in range(getattr(setup,'nvars',3*setup.s+setup.n))]
        solver = z3.Solver(ctx=ctx)
        solver.set(random_seed=solver_seed)
        monomials = {}
        satstart = time.perf_counter()
        for i, m in enumerate(sorted(monomial_support)):
            if i % 128 == 0:
                check(deadline)
            monomials[m] = z3.And(*[variables[j] for j in range(len(variables)) if m >> j & 1]) if m else z3.BoolVal(True, ctx=ctx)
        for eq in equations:
            check(deadline)
            solver.add(z3.Not(xor_tree(z3, [monomials[m] for m in sorted(eq)], ctx=ctx)))
        output['phase_seconds']['sat_setup'] = time.perf_counter()-satstart
    except Budget:
        output.update(status='encoding_timeout', verified=False, rank_gain=0)
        output['total_query_seconds'] = time.perf_counter()-start
        output['counts'] = dict(setup.counts-counts_before)
        return output
    deadline = start+encoding_seconds if query_seconds == 0 else time.perf_counter()+query_seconds
    attempts, sat_ns, lift_ns = 0, 0, 0
    while True:
        left = deadline-time.perf_counter()
        if left <= 0:
            output.update(status='timeout', verified=False, rank_gain=0)
            break
        solver.set(timeout=max(1,int(left*1000)))
        t = time.perf_counter()
        status = solver.check()
        sat_ns += time.perf_counter()-t
        if status == z3.unknown:
            output.update(status='timeout', reason=solver.reason_unknown(), verified=False, rank_gain=0)
            break
        if status == z3.unsat:
            output.update(status='no_decomposition', verified=False, rank_gain=0)
            break
        model = solver.model()
        assignment = sum(int(z3.is_true(model.eval(v, model_completion=True))) << i for i,v in enumerate(variables))
        assert all(evaluate(poly, assignment) == 0 for poly in polys)
        t = time.perf_counter()
        xs = [evaluate(p, assignment) for p in setup.maps]
        witness = lift(setup.curve, xs, tuple(record['target']))
        attempts += 1
        if witness is not None:
            row = setup.row(witness)
            independent = rank_add(setup.matrix, row, setup.curve.r) if record['kind']=='natural' else False
            output.update(status='verified', verified=True, witness=witness, row=row,
                          assignment=assignment, rank_gain=int(independent))
            lift_ns += time.perf_counter()-t
            break
        solver.add(z3.Or(*[v != z3.BoolVal(bool(assignment >> i & 1), ctx=ctx) for i,v in enumerate(variables)]))
        lift_ns += time.perf_counter()-t
    output['phase_seconds'].update(sat_search=sat_ns, lifting_rank_and_blocking=lift_ns)
    output['candidate_models'] = attempts
    output['solver_statistics'] = {k: v for k,v in solver.statistics()}
    output['total_query_seconds'] = time.perf_counter()-start
    output['counts'] = dict(setup.counts-counts_before)
    output['rank_after'] = len(setup.matrix)
    output['process_peak_rss_bytes_before_diagnostic'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    if diagnostic:
        output['groebner'], _ = boolean_groebner(equations, 1)
    return output


def main():
    import argparse
    import z3
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--s', type=int, required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--arm', choices=['linear_payload','quotient_inverse'], required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--diagnostic', action='store_true')
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (1073741824,1073741824))
    started = time.perf_counter()
    setup = Setup(args.n,args.s,args.seed,args.arm)
    setup_seconds = time.perf_counter()-started
    t = time.perf_counter()
    inputs = setup.inputs(1 if args.pilot else 6, 0 if args.pilot else 2)
    target_seconds = time.perf_counter()-t
    header = {'type':'setup','arm':args.arm,'certificate':setup.certificate(),
              'setup_seconds':setup_seconds,'target_preparation_seconds':target_seconds,
              'setup_counts':dict(setup.counts),'z3_version':z3.get_version_string(),
              'python_version':sys.version,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}
    with args.out.open('x') as stream:
        stream.write(json.dumps(header,sort_keys=True)+'\n');stream.flush()
        for item in inputs:
            if args.diagnostic:
                if item['kind']!='natural' or item['index']>=2:
                    continue
                t=time.perf_counter()
                try:
                    polys,equations=build_system(setup,item['target'],t+20)
                    diagnostic,_=boolean_groebner(equations,1)
                    output=dict(item,groebner=diagnostic,anf_degree=max([0]+[m.bit_count() for p in polys for m in p]))
                except Budget:
                    output=dict(item,status='encoding_timeout',groebner=None)
                output['diagnostic_process_peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
            else:
                output = solve_query(setup,item)
            stream.write(json.dumps(output,sort_keys=True)+'\n');stream.flush()
    print(json.dumps({'output':str(args.out),'seconds':time.perf_counter()-started,'rank':len(setup.matrix)}))


if __name__ == '__main__':
    main()
