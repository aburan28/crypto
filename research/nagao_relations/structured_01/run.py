"""Matched frozen/holdout trials and explicit eight-target batch accounting."""
import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CODE = ROOT / 'ecc2k130/codegen'
OLD = HERE.parent / 'subfield_01'
sys.path[:0] = [str(HERE), str(OLD), str(CODE)]
import algebra
import field
import hybrid
import semaev
import cnf
import nagaocompare
import pycryptosat
import space
import solvers

spec = importlib.util.spec_from_file_location('subfield_previous', OLD / 'run.py')
previous = importlib.util.module_from_spec(spec)
spec.loader.exec_module(previous)


def encoded(f, c, v, target, instance, variant):
    pvars = semaev.encode(f, c, f.m, target, instance, variant)
    for xs in pvars:
        coefficients = [instance.newVar() for _ in v.basis]
        for i, x in enumerate(xs):
            instance.addXor([x] + [b for b, mask in zip(coefficients, v.basisCoords) if mask >> i & 1], False)
    return pvars


class Oracle(previous.PairOracle):
    def __init__(self, f, c, v):
        self.f, self.c, self.d = f, c, v.d
        self.base = {}
        for x in v.values[1:]:
            p = c.pointFromX(x)
            if p:
                self.base[f.toCoords(x)] = [p, c.neg(p)]
        invs = dict(zip(v.values[1:], hybrid.batchInverse(f, v.values[1:])))
        self.pairs = {}
        xs = sorted(self.base)
        for i, x in enumerate(xs):
            p = self.base[x][0]
            for y in xs[i + 1:]:
                for q in self.base[y]:
                    dx = f.add(p[0], q[0])
                    lam = f.mul(f.add(p[1], q[1]), invs[dx])
                    xx = f.add(f.add(f.add(f.sqr(lam), lam), dx), c.a2)
                    yy = f.add(f.add(f.mul(lam, f.add(p[0], xx)), xx), p[1])
                    total = xx, yy
                    # Direct affine group law with batched denominator inverses;
                    # all four signs are represented by these two sums and negatives.
                    for point in (total, c.neg(total)):
                        self.pairs.setdefault(point, []).append((x, y))


def context(instance, counted=True):
    f = solvers.Field(instance['n']) if counted else field.Onb(instance['n'])
    a, b = [f.fromCoords(x) for x in instance['coefficients']]
    c = solvers.Curve(f, a, b) if counted else algebra.Curve(f, a, b)
    v = space.Space(f, instance['basis'])
    return f, c, v


def cell(instance, variant, mode, budget):
    start = time.perf_counter()
    deadline = start + budget
    f, c, v = context(instance)
    target = tuple(f.fromCoords(x) for x in instance['target'])
    phases = {'setup': time.perf_counter() - start, 'search': 0., 'verification': 0.}
    solutions = {}
    complete = False
    first = None
    rejected = candidates = 0
    tableMeta = None
    search = None
    try:
        if variant in ('hybrid-reference', 'hybrid-filtered'):
            t = time.perf_counter()
            search = solvers.Search(f, c, v, target, deadline, variant == 'hybrid-filtered')
            phases['setup'] += time.perf_counter() - t
            seen = set()
            iterator = iter(search.candidates())
            while True:
                f.phase = 'search'
                t = time.perf_counter()
                try:
                    a, b, invB, z = next(iterator)
                except StopIteration:
                    phases['search'] += time.perf_counter() - t
                    complete = True
                    break
                except TimeoutError:
                    phases['search'] += time.perf_counter() - t
                    raise
                phases['search'] += time.perf_counter() - t
                candidates += 1
                if (a, b) in seen:
                    continue
                seen.add((a, b))
                f.phase = 'verification'
                t = time.perf_counter()
                result = search.recover(a, b, invB, z)
                if result:
                    solutions[result[0]] = result[1]
                    if first is None:
                        first = time.perf_counter() - start
                phases['verification'] += time.perf_counter() - t
                if solutions and mode == 'first':
                    break
        elif variant == 'pair-invariants-s3':
            t = time.perf_counter()
            table = solvers.PairTable(f, c, v, deadline)
            phases['setup'] += time.perf_counter() - t
            tableMeta = table.metadata()
            f.phase = 'search_and_verification'
            t = time.perf_counter()
            for xs, points in table.solve(target, mode, deadline):
                solutions[xs] = points
                if first is None:
                    first = time.perf_counter() - start
            complete = mode == 'enumerate' or not solutions
            phases['search'] += time.perf_counter() - t
        else:
            t = time.perf_counter()
            cn = cnf.Cnf()
            pvars = encoded(f, c, v, target, cn, variant)
            solver = pycryptosat.Solver(threads=1)
            for clause in cn.clauses:
                solver.add_clause(clause)
            for row, rhs in cn.xors:
                solver.add_xor_clause(row, rhs)
            tableMeta = cn.stats()
            phases['setup'] += time.perf_counter() - t
            while time.perf_counter() < deadline:
                t = time.perf_counter()
                ok, model = solver.solve(time_limit=max(.001, deadline - t))
                phases['search'] += time.perf_counter() - t
                if ok is None:
                    break
                if not ok:
                    complete = True
                    break
                t = time.perf_counter()
                vals = {i: value for i, value in enumerate(model) if i}
                if not nagaocompare.cnfSatisfied(cn, vals):
                    raise ArithmeticError('SAT model fails its CNF')
                xs = tuple(sum(int(vals[x]) << i for i, x in enumerate(vs)) for vs in pvars)
                if len(set(xs)) != 3 or any(f.fromCoords(x) not in v.indices or x == 0 or x == instance['target'][0] for x in xs):
                    raise ArithmeticError('SAT support mismatch')
                f.phase = 'verification'
                points = previous.certificate(f, c, xs, target)
                if points:
                    solutions[xs] = points
                    if first is None:
                        first = time.perf_counter() - start
                else:
                    rejected += 1
                solver.add_clause([-x if vals[x] else x for vs in pvars for x in vs])
                phases['verification'] += time.perf_counter() - t
                if solutions and mode == 'first':
                    break
    except TimeoutError:
        pass
    elapsed = time.perf_counter() - start
    # Timeout work interrupted inside a setup/query phase is retained here.
    phases['other_or_interrupted_phase'] = max(0., elapsed - sum(phases.values()))
    status = 'first' if mode == 'first' and solutions else 'complete' if complete else 'timeout'
    return {'variant': variant, 'mode': mode, 'status': status, 'within_budget': elapsed <= budget,
            'all_phase_seconds': elapsed, 'first_verified_seconds': first, 'phase_seconds': phases,
            'solutions': [list(x) for x in sorted(solutions)],
            'certificates': [{'xs': list(x), 'points': solutions[x]} for x in sorted(solutions)],
            'verified_unique_relations': len(solutions), 'candidates': candidates, 'rejected_sat_models': rejected,
            'early_rejected_functions': search.earlyRejected if search else None, 'table': tableMeta,
            'counters': f.fullReport(), 'common_operations': None, 'full_dlp_S': None,
            'rho_ratio': None, 'speedup': None, 'budget_seconds': budget}


def check(row, instance, oracle):
    target = tuple(oracle.f.fromCoords(x) for x in instance['target'])
    expected = oracle.expected(target)
    got = {tuple(x) for x in row['solutions']}
    if not got <= expected or row['status'] == 'complete' and got != expected:
        raise ArithmeticError('complete-set/oracle mismatch')
    for cert in row['certificates']:
        if sorted(p[0] for p in cert['points']) != cert['xs']:
            raise ArithmeticError('certificate projection mismatch')
        total = None
        for pp in cert['points']:
            p = tuple(oracle.f.fromCoords(x) for x in pp)
            if not oracle.c.onCurve(p):
                raise ArithmeticError('invalid point certificate')
            total = oracle.c.add(total, p)
        if total != target:
            raise ArithmeticError('invalid sum certificate')
    return len(expected)


def batch(instance, targets, oracle, budget):
    start = time.perf_counter()
    f, c, v = context(instance)
    table = solvers.PairTable(f, c, v, start + budget)
    setupSeconds = time.perf_counter() - start
    queries = []
    for item in targets:
        f.phase = 'queries'
        t = time.perf_counter()
        target = tuple(f.fromCoords(x) for x in item['target'])
        solutions = dict(table.solve(target, 'enumerate', start + budget))
        row = {'status': 'complete', 'solutions': [list(x) for x in sorted(solutions)],
               'certificates': [{'xs': list(x), 'points': solutions[x]} for x in sorted(solutions)],
               'all_phase_seconds': time.perf_counter() - t, 'target': item['target'], 'stratum': item['stratum']}
        row['expected_count'] = check(row, {**instance, **item}, oracle)
        queries.append(row)
    return {'variant': 'pair-invariants-s3', 'mode': 'batch-enumerate', 'batch_targets': len(queries),
            'setup_seconds': setupSeconds, 'all_phase_seconds': time.perf_counter() - start,
            'queries': queries, 'counters': f.fullReport(), 'table': table.metadata(),
            'common_operations': None, 'full_dlp_S': None, 'rho_ratio': None, 'speedup': None}


def validate():
    result = semaev.validateIdentity()
    checks = targets = parityChecks = 0
    for kind in ('prefix', 'f4-stable'):
        plain = field.Onb(6)
        alpha = algebra.subfieldGenerator(plain, 2)
        basis, _ = space.makeBasis(plain, 4, kind, alpha)
        instance = {'n': 6, 'd': 4, 'basis': basis, 'coefficients': [plain.toCoords(alpha)] * 2}
        f, c, v = context(instance)
        oracle = Oracle(f, c, v)
        cache = solvers.ImageCache(f, v, True, time.perf_counter() + 30)
        for u, image in cache.images.items():
            for i in range(64):
                value = f.fromCoords(i)
                if image.contains(value) != (image.preimage(value) is not None):
                    raise ArithmeticError('compiled parity filter mismatch')
                parityChecks += 1
        pair = solvers.PairTable(f, c, v, time.perf_counter() + 30)
        # Every S3 pair-table entry agrees with direct group arithmetic.
        for q, entries in pair.table.items():
            for x, y in entries:
                possible = {c.add(px, py)[0] for px in (pair.base[x], c.neg(pair.base[x]))
                            for py in (pair.base[y], c.neg(pair.base[y]))}
                if q not in possible:
                    raise ArithmeticError('pair invariant/group law mismatch')
        # Independent exhaustive signed triples verify the batched group oracle.
        truth = {}
        xs = sorted(oracle.base)
        for i, x in enumerate(xs):
            for j in range(i + 1, len(xs)):
                for z in xs[j + 1:]:
                    for mask in range(8):
                        total = None
                        for bit, xx in enumerate((x, xs[j], z)):
                            total = c.add(total, oracle.base[xx][mask >> bit & 1])
                        if total is not None and f.toCoords(total[0]) not in (x, xs[j], z):
                            truth.setdefault(total, set()).add((x, xs[j], z))
        for x in range(64):
            p = c.pointFromX(f.fromCoords(x))
            if p is None:
                continue
            for target in sorted({p, c.neg(p)}):
                if oracle.expected(target) != truth.get(target, set()):
                    raise ArithmeticError('pair/triple oracle disagreement')
                item = {**instance, 'target': [f.toCoords(x) for x in target]}
                for variant in ('hybrid-reference', 'hybrid-filtered', 'pair-invariants-s3', 'chained-s3', 's4-symmetric'):
                    row = cell(item, variant, 'enumerate', 30)
                    check(row, item, oracle)
                    if row['status'] != 'complete':
                        raise ArithmeticError('tiny validation timeout')
                if kind == 'f4-stable' and row['certificates']:
                    conjugate = tuple(f.frob(x, 2) for x in target)
                    if conjugate != target and 'fixed_target_frobenius_counterexample' not in result:
                        total = None
                        for point in row['certificates'][0]['points']:
                            pp = tuple(f.frob(f.fromCoords(x), 2) for x in point)
                            if pp[0] not in v.indices or not c.onCurve(pp):
                                raise ArithmeticError('Frobenius did not preserve structured support')
                            total = c.add(total, pp)
                        if total != conjugate or total == target:
                            raise ArithmeticError('Frobenius target transformation mismatch')
                        result['fixed_target_frobenius_counterexample'] = {
                            'target': item['target'], 'conjugate_target': [f.toCoords(x) for x in conjugate],
                            'conjugated_factors_remain_in_base': True,
                            'conjugated_relation_is_for_original_target': False}
                search = solvers.Search(f, c, v, target, time.perf_counter() + 30, True, cache)
                reference = solvers.Search(f, c, v, target, time.perf_counter() + 30)
                for args in search.candidates():
                    if search.recover(*args) != reference.recover(*args):
                        raise ArithmeticError('early rejection changed an admissible function')
                    checks += 1
                if kind == 'prefix':
                    old = previous.runCell(6, 2, item['coefficients'], 4, item['target'], 'hybrid-image', 'enumerate', 30)
                    if {tuple(x) for x in old['solutions']} != truth.get(target, set()):
                        raise ArithmeticError('frozen predecessor differs')
                targets += 1
    result.update(exhaustive_target_space_pairs=targets, candidate_rejection_checks=checks,
                  compiled_parity_checks=parityChecks, all_five_solvers_match_triples=True,
                  frozen_prefix_predecessor_matches=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    if args.validate_only:
        print(json.dumps(validate()), flush=True)
        return
    contract = json.loads((HERE / 'contract.json').read_text())
    paths = list(HERE.glob('*.py')) + list(OLD.glob('*.py')) + list(CODE.glob('*.py')) + [HERE / 'contract.json', OLD / 'raw.jsonl']
    provenance = {'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  'python': platform.python_version(), 'pycryptosat': importlib.metadata.version('pycryptosat'),
                  'hardware': platform.uname()._asdict(), 'worker_threads': 1}
    raw = HERE / 'raw.jsonl'
    with raw.open('x') as out:
        def record(row):
            out.write(json.dumps(row) + '\n')
            out.flush()
        record({'kind': 'provenance', **provenance})
        record({'kind': 'validation', **validate()})
        print('validation passed', flush=True)
        old = [json.loads(x) for x in (OLD / 'raw.jsonl').read_text().splitlines()]
        for n in contract['field_degrees']:
            for baseKind, dimensions in contract['bases'].items():
                for d in dimensions:
                    f = field.Onb(n)
                    alpha = algebra.subfieldGenerator(f, 2)
                    t = time.perf_counter()
                    basis, construction = space.makeBasis(f, d, baseKind, alpha)
                    descriptor = {'n': n, 'd': d, 'base_kind': baseKind, 'basis': basis,
                                  'coefficients': [f.toCoords(alpha)] * 2}
                    _, c, v = context(descriptor, False)
                    oracle = Oracle(f, c, v)
                    record({'kind': 'space', **descriptor, 'construction': construction,
                            'invariants': v.metadata(alpha), 'signed_base_size': 2 * len(oracle.base),
                            'construction_and_oracle_seconds': time.perf_counter() - t})
                    if baseKind == 'prefix' and d == 8:
                        targets = [{k: r[k] for k in ('target', 'stratum', 'cohort')} for r in old
                                   if r['kind'] == 'instance' and r['n'] == n and r['d'] == d]
                        for r in targets:
                            r['cohort'] = 'frozen-' + r['cohort']
                    else:
                        rng = random.Random(contract['seed'] + n * 1000 + d * 10 + int(baseKind == 'f4-stable'))
                        targets = []
                        for stratum in ('uniform', 'known_decomposable'):
                            p = previous.uniform(f, c, rng) if stratum == 'uniform' else oracle.supported(rng)
                            targets.append({'target': [f.toCoords(x) for x in p], 'stratum': stratum, 'cohort': 'fresh-holdout'})
                    for item in targets:
                        instance = {**descriptor, **item}
                        digest = hashlib.sha256(json.dumps(instance, sort_keys=True).encode()).hexdigest()
                        expected = oracle.expected(tuple(f.fromCoords(x) for x in item['target']))
                        record({'kind': 'instance', 'instance_sha256': digest, 'expected': [list(x) for x in sorted(expected)], **instance})
                        for mode in ('first', 'enumerate'):
                            variants = contract['variants'][:]
                            random.Random(digest + mode).shuffle(variants)
                            for variant in variants:
                                row = cell(instance, variant, mode, contract['cold_budget_seconds'])
                                count = check(row, instance, oracle)
                                record({'kind': 'trial', 'instance_sha256': digest, 'expected_count': count, **descriptor, **item, **row})
                                print(n, baseKind, d, item['cohort'], item['stratum'], mode, variant,
                                      row['status'], len(row['solutions']), round(row['all_phase_seconds'], 3), flush=True)
                    rng = random.Random(contract['batch_seed'] + n * 1000 + d * 10 + int(baseKind == 'f4-stable'))
                    requests = []
                    for stratum in ('uniform', 'known_decomposable'):
                        for _ in range(contract['batch_targets_per_stratum']):
                            p = previous.uniform(f, c, rng) if stratum == 'uniform' else oracle.supported(rng)
                            requests.append({'target': [f.toCoords(x) for x in p], 'stratum': stratum})
                    row = batch(descriptor, requests, oracle, contract['batch_budget_seconds'])
                    record({'kind': 'batch', **descriptor, **row})
                    print('batch', n, baseKind, d, len(row['queries']), round(row['all_phase_seconds'], 3), flush=True)


if __name__ == '__main__':
    main()
