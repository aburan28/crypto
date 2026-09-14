"""Exact signed-mass, trace-fiber and implementation-work bound audit."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import inspect
import json
from math import comb
from pathlib import Path
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / 'structured_01'))
import run as prior


def fraction(value):
    return {'numerator': value.numerator, 'denominator': value.denominator,
            'decimal': float(value)}


def order(n):
    if n % 2:
        raise ValueError('F4 curve requires even field degree')
    a, b = 2, -1
    for _ in range(1, n // 2):
        a, b = b, -b - 4 * a
    return 2 ** n + 1 - b


def factors(value):
    out = {}
    p = 2
    while p * p <= value:
        while value % p == 0:
            out[p] = out.get(p, 0) + 1
            value //= p
        p += 1
    if value > 1:
        out[value] = out.get(value, 0) + 1
    return out


class Quotient:
    def __init__(self, f, c, k):
        self.f, self.c, self.k = f, c, k
        columns = [f.toCoords(f.add(f.frob(f.fromCoords(1 << i), k), f.fromCoords(1 << i)))
                   for i in range(f.m)]
        subspace = prior.space.Space(f, prior.space.kernel(columns))
        if subspace.d != k:
            raise ArithmeticError('incorrect subfield dimension')
        self.points = [None]
        for x in subspace.values:
            p = c.pointFromX(x)
            if p is not None:
                for q in (p,) if not x else (p, c.neg(p)):
                    if any(f.frob(a, k) != a for a in q):
                        raise ArithmeticError('lift escaped odd-degree coefficient subfield')
                    self.points.append(q)
        if len(self.points) != order(k):
            raise ArithmeticError('Frobenius recurrence disagrees with point count')
        self.indices = {p: i for i, p in enumerate(self.points)}
        self.add = [[self.indices[c.add(p, q)] for q in self.points] for p in self.points]
        self.neg = [self.indices[c.neg(p)] for p in self.points]
        self.cache = {None: 0}

    def trace(self, p):
        if p not in self.cache:
            q = None
            for j in range(self.f.m // self.k):
                q = self.c.add(q, tuple(self.f.frob(x, self.k * j) for x in p))
            if q not in self.indices:
                raise ArithmeticError('group trace escaped quotient')
            self.cache[p] = self.indices[q]
        return self.cache[p]

    def tripleMass(self, base):
        size = len(self.points)
        dp = [[0] * size for _ in range(4)]
        dp[0][0] = 1
        for points in base.values():
            classes = [self.trace(p) for p in points]
            for degree in (3, 2, 1):
                for a, count in enumerate(dp[degree - 1]):
                    if count:
                        for b in classes:
                            dp[degree][self.add[a][b]] += count
        return dp[3]


def masses(oracle, quotient):
    f, c = oracle.f, oracle.c
    allMass = quotient.tripleMass(oracle.base)
    zeroIncidences = 0
    collisions = [0] * len(quotient.points)
    for x, points in oracle.base.items():
        for p in points:
            zeroIncidences += sum(x not in pair for pair in oracle.pairs.get(c.neg(p), []))
            count = sum(x not in pair for pair in oracle.pairs.get(c.neg(c.dbl(p)), []))
            collisions[quotient.trace(c.neg(p))] += count
    if zeroIncidences % 3:
        raise ArithmeticError('infinity multiplicity is not three')
    infinity = zeroIncidences // 3
    valid = [a - b for a, b in zip(allMass, collisions)]
    valid[0] -= infinity
    k = len(oracle.base)
    if sum(allMass) != 8 * comb(k, 3) or min(valid) < 0:
        raise ArithmeticError('incorrect signed triple mass')
    n = order(f.m)
    if n % len(quotient.points):
        raise ArithmeticError('nonintegral trace fiber size')
    h = n // len(quotient.points)
    capacities = [h - (i == 0) for i in range(len(quotient.points))]
    ceiling = Fraction(sum(min(a, b) for a, b in zip(valid, capacities)), n - 1)
    old = min(Fraction(1), Fraction(sum(allMass), n - 1))
    corrected = min(Fraction(1), Fraction(sum(valid), n - 1))
    if not 0 <= ceiling <= corrected <= old <= 1:
        raise ArithmeticError('bound monotonicity failure')
    return {'K_abscissas': k, 'M_signed_points': 2 * k, 'curve_order': n,
            'order_factorization': factors(n), 'all_signed_triples': sum(allMass),
            'infinity_triples': infinity, 'target_x_collision_triples': sum(collisions),
            'valid_signed_triples': sum(valid), 'signed_mean': fraction(Fraction(sum(valid), n - 1)),
            'old_probability_ceiling': fraction(old), 'corrected_probability_ceiling': fraction(corrected),
            'trace_probability_ceiling': fraction(ceiling),
            'trace_over_old_ceiling': fraction(ceiling / old) if old else None,
            'mean_uniform_attempts_lower_bound': fraction(1 / ceiling) if ceiling else None,
            'quotient_order': len(quotient.points), 'fiber_capacity': capacities,
            'all_mass_by_trace': allMass, 'target_x_collisions_by_trace': collisions,
            'valid_mass_by_trace': valid,
            'impossible_trace_classes': [i for i, count in enumerate(valid) if not count]}


def absoluteTrace(f, x):
    total = 0
    for j in range(f.m):
        total = f.add(total, f.frob(x, j))
    if total not in (0, f.one()):
        raise ArithmeticError('invalid absolute trace')
    return int(bool(total))


def branches(f, v, r, traces=None):
    zs = set(v.values) - {0, r}
    tr = absoluteTrace(f, r)
    traces = traces if traces is not None else {h: absoluteTrace(f, h) for h in v.values}
    return sum((2 - (h == r)) * (len(zs) - (h in zs))
               for h in v.values if traces[h] == tr)


def observedBranches(f, c, v, target):
    code = prior.solvers.Search.candidates.__code__
    lines, start = inspect.getsourcelines(code)
    visitLine = start + next(i for i, line in enumerate(lines) if 'k = scalar.polyEval' in line)
    visits = 0
    def trace(frame, event, arg):
        nonlocal visits
        if frame.f_code is code:
            if event == 'line' and frame.f_lineno == visitLine:
                visits += 1
            return trace
        return None
    search = prior.solvers.Search(f, c, v, target, time.perf_counter() + 60)
    previous = sys.gettrace()
    try:
        sys.settrace(trace)
        for _ in search.candidates():
            pass
    finally:
        sys.settrace(previous)
    return visits


def validate(seed):
    f = prior.field.Onb(6)
    alpha = prior.algebra.subfieldGenerator(f, 2)
    c = prior.algebra.Curve(f, alpha, alpha)
    quotient = Quotient(f, c, 6)
    small = Quotient(f, c, 2)
    rng = random.Random(seed)
    bases = [(kind, prior.space.makeBasis(f, 4, kind, alpha)[0]) for kind in ('prefix', 'f4-stable')]
    for i in range(3):
        basis = []
        while len(basis) < 4:
            basis = prior.space.independent(basis + [rng.randrange(1, 64)])
        bases.append(('fresh-' + str(i), basis))
    checks = 0
    rows = []
    for name, basis in bases:
        v = prior.space.Space(f, basis)
        oracle = prior.Oracle(f, c, v)
        result = masses(oracle, quotient)
        allMass = Counter()
        valid = Counter()
        projected = {}
        infinity = collision = 0
        xs = sorted(oracle.base)
        for i, x in enumerate(xs):
            for j in range(i + 1, len(xs)):
                for z in xs[j + 1:]:
                    for mask in range(8):
                        total = None
                        for bit, xx in enumerate((x, xs[j], z)):
                            total = c.add(total, oracle.base[xx][mask >> bit & 1])
                        allMass[quotient.indices[total]] += 1
                        if total is None:
                            infinity += 1
                        elif f.toCoords(total[0]) in (x, xs[j], z):
                            collision += 1
                        else:
                            valid[quotient.indices[total]] += 1
                            projected.setdefault(total, set()).add((x, xs[j], z))
        if result['all_mass_by_trace'] != [allMass[i] for i in range(len(quotient.points))]:
            raise ArithmeticError('DP differs from signed triples')
        if result['valid_mass_by_trace'] != [valid[i] for i in range(len(quotient.points))]:
            raise ArithmeticError('corrected mass differs from signed triples')
        if (result['infinity_triples'], result['target_x_collision_triples']) != (infinity, collision):
            raise ArithmeticError('exclusion formula failure')
        exactProbability = Fraction(len(projected), order(6) - 1)
        if result['trace_probability_ceiling'] != fraction(exactProbability):
            raise ArithmeticError('identity quotient is not exact')
        for target in quotient.points[1:]:
            if oracle.expected(target) != projected.get(target, set()):
                raise ArithmeticError('projected oracle differs from signed triples')
            if observedBranches(f, c, v, target) != branches(f, v, target[0]):
                raise ArithmeticError('hybrid branch formula failed')
            checks += 1
        rows.append({'base_kind': name, 'basis': basis, 'exact_probability': fraction(exactProbability),
                     'signed_vs_projected_duplicate_count': sum(valid.values()) - sum(map(len, projected.values())),
                     'bounds': result})
    return {'tiny_target_space_checks': checks, 'base_curve_order': len(small.points),
            'GF64_curve_order': len(quotient.points), 'fresh_spaces': 3,
            'failures': 0, 'spaces': rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    contract = json.loads((HERE / 'contract.json').read_text())
    if args.validate_only:
        result = validate(contract['seed'])
        print(json.dumps({k: v for k, v in result.items() if k != 'spaces'}))
        return
    oldPath = ROOT / contract['input']
    if hashlib.sha256(oldPath.read_bytes()).hexdigest() != contract['input_sha256']:
        raise ArithmeticError('frozen measurement changed')
    old = [json.loads(s) for s in oldPath.read_text().splitlines()]
    sources = [HERE / 'run.py', HERE / 'contract.json', HERE / 'PROOFS.md']
    sources += [ROOT / name for name in old[0]['sha256']]
    provenance = {'kind': 'provenance', 'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                  'prior_raw_sha256': contract['input_sha256'], 'scope': contract['scope']}
    start = time.perf_counter()
    with (HERE / 'raw.jsonl').open('x') as out:
        def record(row):
            out.write(json.dumps(row) + '\n')
            out.flush()
        record(provenance)
        validation = validate(contract['seed'])
        record({'kind': 'validation', **validation})
        print('validation passed', validation['tiny_target_space_checks'], flush=True)
        for descriptor in [r for r in old if r['kind'] == 'space']:
            begin = time.perf_counter()
            f, c, v = prior.context(descriptor, False)
            oracle = prior.Oracle(f, c, v)
            quotient = Quotient(f, c, contract['trace_quotient_degree'])
            result = masses(oracle, quotient)
            key = tuple(descriptor[k] for k in ('n', 'base_kind', 'd'))
            batch = next(r for r in old if r['kind'] == 'batch' and tuple(r[k] for k in ('n', 'base_kind', 'd')) == key)
            k, m = result['K_abscissas'], result['M_signed_points']
            skips = sum(f.fromCoords(q['target'][0]) in {p[0][0] for p in oracle.base.values()} for q in batch['queries'])
            floor = 7 * comb(k, 2) + 2 * (8 * m - 2 * skips)
            measured = batch['counters']['field_api']['totals']['multiplications']
            if measured < floor or batch['table']['table_entries'] != k * (k - 1):
                raise ArithmeticError('saved S3 batch violates work/space floor')
            targets = [(r['target'], r['cohort']) for r in old if r['kind'] == 'instance' and tuple(r[a] for a in ('n', 'base_kind', 'd')) == key]
            targets += [(q['target'], 'frozen-batch-' + q['stratum']) for q in batch['queries']]
            rng = random.Random(contract['seed'] + descriptor['n'] * 1000 + descriptor['d'] * 10 + (descriptor['base_kind'] == 'f4-stable'))
            for _ in range(contract['fresh_targets_per_space']):
                p = prior.previous.uniform(f, c, rng)
                targets.append(([f.toCoords(x) for x in p], 'fresh-bound-holdout'))
            checks = []
            traces = {h: absoluteTrace(f, h) for h in v.values}
            for target, cohort in targets:
                p = tuple(f.fromCoords(x) for x in target)
                expected = oracle.expected(p)
                traceClass = quotient.trace(p)
                if len(expected) > result['valid_mass_by_trace'][traceClass]:
                    raise ArithmeticError('oracle exceeds signed fiber mass')
                predicted = branches(f, v, p[0], traces)
                checks.append({'target': target, 'cohort': cohort, 'trace_class': traceClass,
                               'exact_projected_relations': len(expected), 'hybrid_complete_branches': predicted,
                               'hybrid_complete_multiplication_floor': 7 * predicted})
            n = result['curve_order']
            mean = Fraction(result['valid_signed_triples'], n - 1)
            lower = 2 * (Fraction(m) - Fraction(2 * m, n - 1)) / mean if mean else None
            hybridFloor = sum(r['hybrid_complete_multiplication_floor'] for r in checks if r['cohort'].startswith('frozen-batch-'))
            record({'kind': 'space_bounds', **{a: descriptor[a] for a in ('n', 'base_kind', 'd', 'basis', 'coefficients')},
                    **result, 'target_checks': checks, 'batch_multiplication_floor': floor,
                    'batch_measured_multiplications': measured, 'batch_measured_over_floor': fraction(Fraction(measured, floor)),
                    'same_batch_hybrid_multiplication_floor': hybridFloor,
                    'hybrid_floor_over_measured_s3_multiplications': fraction(Fraction(hybridFloor, measured)),
                    'uniform_query_multiplications_per_output_lower_bound': fraction(lower) if lower else None,
                    'quotient_points': [None if p is None else [f.toCoords(x) for x in p] for p in quotient.points],
                    'harness_seconds': time.perf_counter() - begin,
                    'common_operation_speedup': None, 'full_dlp_S': None, 'rho_ratio': None})
            print(*key, 'K', k, 'ceilings', result['old_probability_ceiling']['decimal'], result['trace_probability_ceiling']['decimal'], 'zero classes', len(result['impossible_trace_classes']), flush=True)
        record({'kind': 'completion', 'harness_seconds': time.perf_counter() - start, 'failures': 0})


if __name__ == '__main__':
    main()
