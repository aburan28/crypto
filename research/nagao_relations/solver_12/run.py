"""The null-object control: brute-force pair enumeration on solver_10's instances.

solver_10 established that the Riemann-Roch solvers beat the Semaev controls on
the real ECC2K-130 curve.  That is a comparison between two algebraic encodings
under one SAT solver.  It is not evidence that either beats the dumbest oracle
that could possibly work: enumerate pairs of factor-base points and look at the
residual.

This runs that oracle on solver_10's exact targets, in solver_10's units, and
checks both oracles against each other on every instance solver_10 completed.
"""
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CODE = ROOT / 'ecc2k130/codegen'
sys.path.insert(0, str(CODE))
import curves
import nagaoannihilator as scalar

CHALLENGE_R = 680564733841876926932320129493409985129
PANEL = HERE.parent / 'solver_10/raw.jsonl'


def factorBase(f, curve, d):
    """Both lifts of every admissible abscissa in V \\ {0}."""
    base = []
    for i in range(1, 1 << d):
        p = curve.pointFromX(f.fromCoords(i))
        if p is None:
            continue
        neg = curve.neg(p)
        if not curve.onCurve(p) or not curve.onCurve(neg):
            raise ArithmeticError('factor-base point off curve')
        base.append(p)
        if neg != p:
            base.append(neg)
    return base


def pairOracle(f, curve, d, target, base, deadline):
    """Every unordered pair, residual tested for factor-base membership."""
    r = target[0]
    xr = f.toCoords(r)
    solutions = set()
    pairs = 0
    first = None
    start = time.perf_counter()
    for i in range(len(base)):
        if time.perf_counter() >= deadline:
            return solutions, pairs, first, False
        p1 = base[i]
        x1 = f.toCoords(p1[0])
        for j in range(i + 1, len(base)):
            p2 = base[j]
            x2 = f.toCoords(p2[0])
            if x1 == x2:
                continue
            pairs += 1
            q = curve.add(target, curve.neg(curve.add(p1, p2)))
            if q is None:
                continue
            x3 = f.toCoords(q[0])
            if x3 == 0 or x3 >= 1 << d or x3 in (x1, x2, xr):
                continue
            if x1 == xr or x2 == xr:
                continue
            solutions.add(tuple(sorted((x1, x2, x3))))
            if first is None:
                first = time.perf_counter() - start
    return solutions, pairs, first, True


def verify(f, curve, d, target, xs):
    """Re-derive the relation from the abscissas alone, as solver_10 does."""
    if len(set(xs)) != 3:
        raise ArithmeticError('repeated abscissa')
    if any(x == 0 or x >= 1 << d or x == f.toCoords(target[0]) for x in xs):
        raise ArithmeticError('abscissa outside the factor base or excluded')
    points = []
    for x in xs:
        p = curve.pointFromX(f.fromCoords(x))
        if p is None:
            raise ArithmeticError('abscissa is not on the curve')
        points.append(p)
    for signs in range(8):
        total = None
        for i, p in enumerate(points):
            total = curve.add(total, curve.neg(p) if signs >> i & 1 else p)
        if total == target:
            return True
    raise ArithmeticError('signed sum does not reach the target')


def main():
    raw = HERE / 'raw.jsonl'
    if raw.exists():
        raise SystemExit('Evidence exists; choose a new sibling folder')
    contract = json.loads((HERE / 'contract.json').read_text())
    panelRows = [json.loads(x) for x in PANEL.read_text().splitlines()]
    sources = list(CODE.glob('*.py')) + [HERE / 'run.py', HERE / 'contract.json', PANEL]
    provenance = {
        'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'command': [sys.executable] + sys.argv,
        'python': platform.python_version(),
        'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    }
    rows = []
    with raw.open('x') as out:
        def record(row):
            rows.append(row)
            out.write(json.dumps(row) + '\n')
            out.flush()

        record({'kind': 'provenance', **provenance})
        for panel in contract['panels']:
            n, d = panel['n'], panel['d']
            if curves.curveOrder(n) != 4 * CHALLENGE_R:
                raise SystemExit('group order mismatch: not the ECC2K-130 curve')
            source = next(r for r in panelRows if r['kind'] == 'panel' and r['d'] == d)
            for stratum, entries in source['targets'].items():
                for entry in entries:
                    f = scalar.CountedField(n)
                    curve = curves.Curve(f)
                    f.phase = 'setup'
                    base = factorBase(f, curve, d)
                    target = tuple(f.fromCoords(v) for v in entry['target'])
                    f.phase = 'search'
                    start = time.perf_counter()
                    solutions, pairs, first, complete = pairOracle(
                        f, curve, d, target, base, start + contract.get('seconds_per_instance', 600))
                    elapsed = time.perf_counter() - start
                    f.phase = 'verify'
                    for xs in solutions:
                        verify(f, curve, d, target, xs)
                    planted = entry['planted']
                    if planted is not None and complete and tuple(planted) not in solutions:
                        raise ArithmeticError('complete pair enumeration missed the planted triple')
                    # Cross-oracle equality against solver_10's completed runs.
                    cross = None
                    for t in panelRows:
                        if (t['kind'] == 'trial' and t['d'] == d and t['mode'] == 'enumerate'
                                and t['variant'] == 'quadratic-image' and t['status'] == 'complete'
                                and t['target'] == entry['target']):
                            theirs = {tuple(x) for x in t['solutions']}
                            if complete and theirs != solutions:
                                raise ArithmeticError('pair oracle and quadratic-image disagree on a complete set')
                            cross = {'compared_against': 'quadratic-image', 'sets_equal': theirs == solutions,
                                     'their_count': len(theirs), 'our_count': len(solutions)}
                    record({'kind': 'trial', 'n': n, 'd': d, 'stratum': stratum,
                            'target': entry['target'], 'planted': planted,
                            'factor_base_points': len(base), 'pairs_enumerated': pairs,
                            'status': 'complete' if complete else 'timeout',
                            'solutions': [list(x) for x in sorted(solutions)],
                            'verified_unique_relations': len(solutions),
                            'first_verified_seconds': first, 'all_phase_seconds': elapsed,
                            'field_api_counts': f.report(), 'cross_oracle': cross})
                    print(d, stratum, 'pairs=%d' % pairs, 'rel=%d' % len(solutions),
                          '%.3fs' % elapsed, 'cross=%s' % (cross and cross['sets_equal']), flush=True)

    trials = [r for r in rows if r['kind'] == 'trial']
    groups = []
    for panel in contract['panels']:
        d = panel['d']
        for st in ('uniform', 'known_decomposable'):
            rr = [r for r in trials if r['d'] == d and r['stratum'] == st]
            if not rr:
                continue
            ours = sum(r['field_api_counts']['totals']['fieldOperations'] for r in rr)
            theirs = 0
            for r in rr:
                for t in panelRows:
                    if (t['kind'] == 'trial' and t['d'] == d and t['mode'] == 'enumerate'
                            and t['variant'] == 'quadratic-image' and t['target'] == r['target']):
                        theirs += t['field_api_counts']['totals']['fieldOperations']
            groups.append({'n': panel['n'], 'd': d, 'stratum': st, 'attempted': len(rr),
                           'complete': sum(r['status'] == 'complete' for r in rr),
                           'verified_relations': sum(r['verified_unique_relations'] for r in rr),
                           'pairs_enumerated': sum(r['pairs_enumerated'] for r in rr),
                           'pair_oracle_ops': ours,
                           'quadratic_image_ops': theirs,
                           'pair_over_image': ours / theirs if theirs else None,
                           'all_phase_seconds': sum(r['all_phase_seconds'] for r in rr)})
    (HERE / 'summary.json').write_text(json.dumps(
        {'provenance': provenance, 'groups': groups,
         'cross_oracle_disagreements': 0,
         'limits': contract['limits']}, indent=2) + '\n')
    for g in groups:
        print('d=%d %s pair/image ops ratio %s' % (g['d'], g['stratum'], g['pair_over_image']), flush=True)


if __name__ == '__main__':
    main()
