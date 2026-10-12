"""Frozen exploratory subfield and dimension panel; no full-DLP cost claim."""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import random
import resource
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CODE = ROOT / 'ecc2k130/codegen'
sys.path.insert(0, str(CODE))
import field
import cnf
import nagaocompare
import nagaoannihilator as scalar
import pycryptosat
import algebra
import hybrid
import semaev


def coords(f, p):
    return None if p is None else [f.toCoords(v) for v in p]


def makeCurve(n, k, coefficients=None, counted=False):
    f = scalar.CountedField(n) if counted else field.Onb(n)
    if coefficients is None:
        t = algebra.subfieldGenerator(f, k)
        coefficients = [f.toCoords(t)] * 2
    c = algebra.Curve(f, *(f.fromCoords(v) for v in coefficients))
    return f, c


def certificate(f, c, xs, target):
    points = [c.pointFromX(f.fromCoords(x)) for x in xs]
    if any(p is None for p in points):
        return None
    for mask in range(8):
        signed = [c.neg(p) if mask >> i & 1 else p for i, p in enumerate(points)]
        total = None
        for p in signed:
            if not c.onCurve(p):
                raise ArithmeticError('off-curve lift')
            total = c.add(total, p)
        if total == target:
            return [coords(f, p) for p in signed]
    return None


class PairOracle:
    """Group-law pair lookup only; no RR coefficients or Semaev equations."""
    def __init__(self, f, c, d):
        self.f, self.c, self.d = f, c, d
        self.base = {}
        for x in range(1, 1 << d):
            p = c.pointFromX(f.fromCoords(x))
            if p is not None:
                self.base[x] = [p, c.neg(p)]
        self.pairs = {}
        xs = sorted(self.base)
        for i, x in enumerate(xs):
            for y in xs[i + 1:]:
                for p in self.base[x]:
                    for q in self.base[y]:
                        self.pairs.setdefault(c.add(p, q), []).append((x, y))

    def expected(self, target):
        r = self.f.toCoords(target[0])
        out = set()
        for z, points in self.base.items():
            if z == r:
                continue
            for p in points:
                residual = self.c.add(target, self.c.neg(p))
                for x, y in self.pairs.get(residual, []):
                    if z not in (x, y) and r not in (x, y):
                        out.add(tuple(sorted((x, y, z))))
        return out

    def supported(self, rng):
        if len(self.base) < 3:
            raise ValueError('base too small for supported targets')
        while True:
            xs = rng.sample(sorted(self.base), 3)
            total = None
            for x in xs:
                total = self.c.add(total, self.base[x][rng.randrange(2)])
            if total is not None and self.f.toCoords(total[0]) not in xs:
                return total


def uniform(f, c, rng):
    # Rejection from all (x,sign); the unique x=0 point gets one slot.
    while True:
        x, sign = rng.randrange(1 << f.m), rng.randrange(2)
        if not x and sign:
            continue
        p = c.pointFromX(f.fromCoords(x))
        if p is not None:
            return c.neg(p) if sign else p


def runCell(n, k, coefficients, d, targetCoords, variant, mode, budget):
    start = time.perf_counter()
    deadline = start + budget
    f, curve = makeCurve(n, k, coefficients, counted=True)
    target = tuple(f.fromCoords(v) for v in targetCoords)
    if not curve.onCurve(target):
        raise ValueError('target not on curve')
    solutions = {}
    complete = False
    first = None
    candidates = duplicates = rejected = calls = 0
    phases = {'setup': time.perf_counter() - start, 'search': 0.0,
              'support_extract_verify': 0.0}
    cnfStats = None
    try:
        if variant == 'hybrid-image':
            before = time.perf_counter()
            search = hybrid.Search(f, curve, d, target, deadline)
            phases['setup'] += time.perf_counter() - before
            seen = set()
            it = iter(search.candidates())
            while True:
                f.phase = 'search'
                before = time.perf_counter()
                try:
                    a, b, invB, z = next(it)
                except StopIteration:
                    phases['search'] += time.perf_counter() - before
                    complete = True
                    break
                except TimeoutError:
                    phases['search'] += time.perf_counter() - before
                    raise
                phases['search'] += time.perf_counter() - before
                candidates += 1
                if (a, b) in seen:
                    duplicates += 1
                    continue
                seen.add((a, b))
                f.phase = 'support_extract_verify'
                before = time.perf_counter()
                found = search.recover(a, b, invB, z)
                if found is not None:
                    xs, points = found
                    solutions[xs] = points
                    if first is None:
                        first = time.perf_counter() - start
                phases['support_extract_verify'] += time.perf_counter() - before
                if solutions and mode == 'first':
                    break
        else:
            before = time.perf_counter()
            c = cnf.Cnf()
            pvars = semaev.encode(f, curve, d, target, c, variant)
            solver = pycryptosat.Solver(threads=1)
            for clause in c.clauses:
                solver.add_clause(clause)
            for row, rhs in c.xors:
                solver.add_xor_clause(row, rhs)
            cnfStats = c.stats()
            phases['setup'] += time.perf_counter() - before
            while time.perf_counter() < deadline:
                before = time.perf_counter()
                ok, model = solver.solve(time_limit=max(.001, deadline - before))
                phases['search'] += time.perf_counter() - before
                calls += 1
                if ok is None:
                    break
                if not ok:
                    complete = True
                    break
                before = time.perf_counter()
                f.phase = 'support_extract_verify'
                vals = {i: value for i, value in enumerate(model) if i}
                if not nagaocompare.cnfSatisfied(c, vals):
                    raise ArithmeticError('CNF model failed verification')
                xs = tuple(sum(int(vals[v]) << i for i, v in enumerate(vs)) for vs in pvars)
                if xs != tuple(sorted(set(xs))) or len(xs) != 3 or any(
                        x == 0 or x >= 1 << d or x == targetCoords[0] for x in xs):
                    raise ArithmeticError('Semaev domain mismatch')
                points = certificate(f, curve, xs, target)
                candidates += 1
                if points is not None:
                    if xs in solutions:
                        duplicates += 1
                    solutions[xs] = points
                    if first is None:
                        first = time.perf_counter() - start
                else:
                    rejected += 1
                # Block the projection, including algebraic roots with no signed lift.
                solver.add_clause([-v if vals[v] else v for vs in pvars for v in vs])
                phases['support_extract_verify'] += time.perf_counter() - before
                if solutions and mode == 'first':
                    break
    except TimeoutError:
        pass
    elapsed = time.perf_counter() - start
    phases['orchestration'] = max(0.0, elapsed - sum(phases.values()))
    status = 'first' if mode == 'first' and solutions else 'complete' if complete else 'timeout'
    return {'n': n, 'k': k, 'd': d, 'coefficients': coefficients, 'target': targetCoords,
            'variant': variant, 'mode': mode, 'budget_seconds': budget, 'status': status,
            'within_budget': elapsed <= budget, 'solutions': [list(x) for x in sorted(solutions)],
            'certificates': [{'xs': list(x), 'points': solutions[x]} for x in sorted(solutions)],
            'verified_unique_relations': len(solutions), 'first_verified_seconds': first,
            'all_phase_seconds': elapsed, 'phase_seconds': phases, 'candidate_count': candidates,
            'duplicates': duplicates, 'rejected_candidates': rejected, 'solver_calls': calls,
            'cnf': cnfStats, 'field_api_counts': f.report(),
            'field_counter_scope': 'complete_field_API' if variant == 'hybrid-image' else 'setup_and_lifting_only_excludes_SAT_and_encoding',
            'peak_process_rss_bytes_cumulative': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
            'common_operations': None, 'solver_operations': None, 'full_dlp_S': None,
            'rho_ratio': None, 'floor_ratio': None, 'speedup': None}


def checkRow(row, f, c, expected):
    got = {tuple(x) for x in row['solutions']}
    if not got <= expected or (row['status'] == 'complete' and got != expected):
        raise ArithmeticError('independent oracle mismatch')
    for cert in row['certificates']:
        total = None
        for coordsP in cert['points']:
            p = tuple(f.fromCoords(v) for v in coordsP)
            if not c.onCurve(p):
                raise ArithmeticError('certificate curve equation mismatch')
            total = c.add(total, p)
        if coords(f, total) != row['target']:
            raise ArithmeticError('certificate signed sum mismatch')


def validate():
    result = semaev.validateIdentity()
    f, c = makeCurve(6, 2)
    # All x,y pairs independently check the lift set, including the x=0 point.
    brute = {(f.fromCoords(x), f.fromCoords(y)) for x in range(64) for y in range(64)
             if c.onCurve((f.fromCoords(x), f.fromCoords(y)))}
    lifted = set()
    for x in range(64):
        p = c.pointFromX(f.fromCoords(x))
        if p is not None:
            lifted.update((p, c.neg(p)))
    if brute != lifted:
        raise ArithmeticError('brute curve/lift mismatch')
    for n, k in ((6, 2), (9, 3)):
        ff, cc = makeCurve(n, k)
        for v in range(1 << n):
            value = ff.fromCoords(v)
            roots = cc.asRoots(value)
            if bool(roots) != (ff.trace(value) == 0):
                raise ArithmeticError('AS solvability mismatch')
            for w in roots:
                if ff.add(ff.sqr(w), w) != value:
                    raise ArithmeticError('AS root mismatch')
    oracle = PairOracle(f, c, 4)
    # Exhaustive signed triples are independent of the pair lookup structure.
    truth = {}
    xs = sorted(oracle.base)
    for i, x in enumerate(xs):
        for j in range(i + 1, len(xs)):
            for z in xs[j + 1:]:
                triple = (x, xs[j], z)
                for mask in range(8):
                    total = None
                    for bit, xx in enumerate(triple):
                        total = c.add(total, oracle.base[xx][mask >> bit & 1])
                    if total is not None and f.toCoords(total[0]) not in triple:
                        truth.setdefault(total, set()).add(triple)
    supportChecks = 0
    coeffs = [f.toCoords(c.a2), f.toCoords(c.a6)]
    lv = scalar.subspacePolynomial(f, 4)
    for target in sorted(brute):
        expected = oracle.expected(target)
        if expected != truth.get(target, set()):
            raise ArithmeticError('pair/triple oracle mismatch')
        for variant in ('hybrid-image', 'chained-s3', 's4-symmetric'):
            row = runCell(6, 2, coeffs, 4, coords(f, target), variant, 'enumerate', 30)
            checkRow(row, f, c, expected)
            if row['status'] != 'complete':
                raise ArithmeticError('tiny exhaustive validation timed out')
        search = hybrid.Search(f, c, 4, target, time.perf_counter() + 30)
        for a, b, invB, z in search.candidates():
            found = search.recover(a, b, invB, z)
            h, _ = algebra.residualNorm(f, c, target, a, b)
            remainder = scalar.annihilatorRemainder(f, lv, h)
            valid = not any(remainder) and bool(h[0]) and bool(scalar.polyEval(f, h, target[0]))
            if (found is not None) != valid:
                raise ArithmeticError('image/modular support mismatch')
            supportChecks += 1
    result.update({'all_affine_six_bit_targets': len(brute), 'signed_base_size': 2 * len(xs),
                   'all_three_solvers_equal_exhaustive_triples': True,
                   'support_candidate_checks': supportChecks, 'AS_inputs_checked': 64 + 512,
                   'brute_curve_lift_pairs_checked': 64 * 64})
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    if args.validate_only:
        print(json.dumps(validate()), flush=True)
        return
    contract = json.loads((HERE / 'contract.json').read_text())
    paths = list(CODE.glob('*.py')) + list(HERE.glob('*.py')) + [HERE / 'contract.json']
    provenance = {'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  'python': platform.python_version(), 'pycryptosat': importlib.metadata.version('pycryptosat'),
                  'hardware': platform.uname()._asdict(), 'worker_threads': 1,
                  'cpu_count_available': os.cpu_count(), 'memory_limit_bytes': None,
                  'command': [sys.executable] + sys.argv}
    rows = []
    with (HERE / 'raw.jsonl').open('x') as out:
        def record(row):
            rows.append(row)
            out.write(json.dumps(row) + '\n')
            out.flush()
        record({'kind': 'provenance', **provenance})
        result = validate()
        record({'kind': 'validation', **result})
        print('validation', result, flush=True)
        for n in contract['field_degrees']:
            k = contract['coefficient_subfield_degree']
            f, c = makeCurve(n, k)
            coeffs = [f.toCoords(c.a2), f.toCoords(c.a6)]
            t = time.perf_counter()
            small = PairOracle(f, c, min(contract['dimensions']))
            targets = []
            for cohort in ('development', 'holdout'):
                rng = random.Random(contract['seeds'][cohort] + n)
                for stratum in ('uniform', 'known_decomposable'):
                    target = uniform(f, c, rng) if stratum == 'uniform' else small.supported(rng)
                    targets.append((cohort, stratum, target))
            record({'kind': 'target_generation', 'n': n, 'k': k, 'coefficients': coeffs,
                    'setup_validation_seconds': time.perf_counter() - t,
                    'targets': [{'cohort': co, 'stratum': st, 'target': coords(f, p)} for co, st, p in targets]})
            for d in contract['dimensions']:
                t = time.perf_counter()
                oracle = small if d == min(contract['dimensions']) else PairOracle(f, c, d)
                record({'kind': 'oracle', 'n': n, 'd': d, 'signed_base_size': 2 * len(oracle.base),
                        'base_abscissas': sorted(oracle.base),
                        'pair_entries': sum(len(v) for v in oracle.pairs.values()),
                        'setup_validation_seconds': time.perf_counter() - t})
                for cohort, stratum, target in targets:
                    expected = oracle.expected(target)
                    targetCoords = coords(f, target)
                    instance = {'n': n, 'k': k, 'd': d, 'coefficients': coeffs, 'target': targetCoords,
                                'cohort': cohort, 'stratum': stratum,
                                'expected': [list(x) for x in sorted(expected)]}
                    instanceHash = hashlib.sha256(json.dumps(instance, sort_keys=True).encode()).hexdigest()
                    record({'kind': 'instance', 'sha256': instanceHash, **instance})
                    for mode in ('first', 'enumerate'):
                        variants = contract['variants'][:]
                        random.Random(instanceHash + mode).shuffle(variants)
                        for variant in variants:
                            row = runCell(n, k, coeffs, d, targetCoords, variant, mode, contract['budget_seconds'])
                            checkRow(row, f, c, expected)
                            record({'kind': 'trial', 'cohort': cohort, 'stratum': stratum,
                                    'instance_sha256': instanceHash, 'expected_count': len(expected), **row})
                            print(n, d, cohort, stratum, mode, variant, row['status'],
                                  len(row['solutions']), round(row['all_phase_seconds'], 3), flush=True)
                    if d == max(contract['dimensions']):
                        row = runCell(n, k, coeffs, d, targetCoords, 'hybrid-image', 'enumerate',
                                      contract['supplemental_hybrid_budget_seconds'])
                        checkRow(row, f, c, expected)
                        record({'kind': 'supplemental', 'cohort': cohort, 'stratum': stratum,
                                'instance_sha256': instanceHash, 'expected_count': len(expected), **row})
                        print('supplemental', n, d, cohort, stratum, row['status'], len(row['solutions']),
                              round(row['all_phase_seconds'], 3), flush=True)
    groups = []
    for kind in ('trial', 'supplemental'):
        keys = sorted({(r['n'], r['d'], r['variant'], r['mode'], r['stratum']) for r in rows if r['kind'] == kind})
        for n, d, variant, mode, stratum in keys:
            rs = [r for r in rows if r['kind'] == kind and
                  (r['n'], r['d'], r['variant'], r['mode'], r['stratum']) == (n, d, variant, mode, stratum)]
            groups.append({'kind': kind, 'n': n, 'd': d, 'variant': variant, 'mode': mode, 'stratum': stratum,
                           'attempts': len(rs), 'resolved': sum(r['status'] != 'timeout' for r in rs),
                           'resolved_within_budget': sum(r['status'] != 'timeout' and r['within_budget'] for r in rs),
                           'verified_relations': sum(r['verified_unique_relations'] for r in rs),
                           'all_phase_seconds': sum(r['all_phase_seconds'] for r in rs),
                           'field_api_operations': sum(r['field_api_counts']['totals']['fieldOperations'] for r in rs)
                           if variant == 'hybrid-image' else None,
                           'common_operations': None, 'speedup': None, 'full_dlp_S': None, 'rho_ratio': None})
    (HERE / 'summary.json').write_text(json.dumps({'provenance': provenance, 'validation': result,
                                                'groups': groups, 'limits': contract['limits']}, indent=2) + '\n')


if __name__ == '__main__':
    main()
