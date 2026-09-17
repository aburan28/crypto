"""The method ceiling: does the RR encoding move the exponent, or only the constant?

Finding P1+P2+P3 = R with every Pi in a stored factor base F is 3SUM over a
group.  Pair enumeration realises the generic Theta(|F|^2) bound.  This asks
whether the Riemann-Roch encoding does anything an exponent could notice, by
measuring its candidate count against |F|^2/2 across five dimensions.
"""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CODE = ROOT / 'ecc2k130/codegen'
sys.path.insert(0, str(CODE))
for folder in ('solver_07', 'solver_08'):
    sys.path.insert(0, str(HERE.parent / folder))
import curves
import image_solver
import nagaoannihilator as scalar

CHALLENGE_R = 680564733841876926932320129493409985129
SEED = 20260917


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


control = load('control12', HERE.parent / 'solver_12/run.py')


def plantedTargets(f, curve, d, rng, count):
    admissible = []
    for i in range(1, 1 << d):
        p = curve.pointFromX(f.fromCoords(i))
        if p is not None:
            admissible.append((i, p))
    out = []
    seen = set()
    guard = 0
    while len(out) < count and guard < 10000:
        guard += 1
        triple = tuple(sorted(x for x, _ in rng.sample(admissible, 3)))
        if triple in seen:
            continue
        total = None
        for x in triple:
            total = curve.add(total, dict(admissible)[x])
        if total is None or f.toCoords(total[0]) < 1 << d:
            continue
        seen.add(triple)
        out.append((total, triple))
    if len(out) < count:
        raise ArithmeticError('could not mint enough planted targets at d=%d' % d)
    return out


def fitExponent(xs, ys):
    """Least-squares slope of log(y) against log(x)."""
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    n = len(lx)
    mx = sum(lx) / n
    my = sum(ly) / n
    num = sum((a - mx) * (b - my) for a, b in zip(lx, ly))
    den = sum((a - mx) ** 2 for a in lx)
    return num / den


def main():
    raw = HERE / 'raw.jsonl'
    if raw.exists():
        raise SystemExit('Evidence exists; choose a new sibling folder')
    contract = json.loads((HERE / 'contract.json').read_text())
    sources = list(CODE.glob('*.py')) + [HERE / 'run.py', HERE / 'contract.json'] \
        + list((HERE.parent / 'solver_08').glob('*.py')) + [HERE.parent / 'solver_12/run.py']
    provenance = {
        'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'command': [sys.executable] + sys.argv,
        'python': platform.python_version(),
        'seed': SEED,
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
            f = scalar.CountedField(n)
            curve = curves.Curve(f)
            rng = random.Random(SEED + d)
            targets = plantedTargets(f, curve, d, rng, 2)
            for target, triple in targets:
                coords = [f.toCoords(target[0]), f.toCoords(target[1])]
                budget = 60 * 4 ** max(0, d - 6) + 60

                # RR side: the frozen image solver, run to completion.
                rr = image_solver.cell(n, d, coords, 'enumerate', budget)

                # Control side: pair enumeration in the same counted unit.
                cf = scalar.CountedField(n)
                ccurve = curves.Curve(cf)
                cf.phase = 'setup'
                base = control.factorBase(cf, ccurve, d)
                ctarget = tuple(cf.fromCoords(v) for v in coords)
                cf.phase = 'search'
                t0 = time.perf_counter()
                sols, pairs, first, complete = control.pairOracle(
                    cf, ccurve, d, ctarget, base, t0 + budget)
                pairSeconds = time.perf_counter() - t0
                cf.phase = 'verify'
                for xs in sols:
                    control.verify(cf, ccurve, d, ctarget, xs)

                rrSols = {tuple(x) for x in rr['solutions']}
                if rr['status'] == 'complete' and complete and rrSols != sols:
                    raise ArithmeticError('oracles disagree on a complete set at d=%d' % d)
                if complete and tuple(triple) not in sols:
                    raise ArithmeticError('pair enumeration missed the planted triple')

                rrOps = rr['field_api_counts']['totals']['fieldOperations']
                pairOps = cf.report()['totals']['fieldOperations']
                row = {
                    'kind': 'trial', 'n': n, 'd': d, 'target': coords, 'planted': list(triple),
                    'factor_base_points': len(base), 'pair_space': len(base) * (len(base) - 1) // 2,
                    'pairs_enumerated': pairs, 'pair_status': 'complete' if complete else 'timeout',
                    'pair_ops': pairOps, 'pair_seconds': pairSeconds, 'pair_relations': len(sols),
                    'rr_status': rr['status'], 'rr_candidates': rr['candidate_functions'],
                    'rr_ops': rrOps, 'rr_seconds': rr['all_phase_seconds'],
                    'rr_relations': rr['verified_unique_relations'],
                    'rr_candidates_over_half_f_squared': rr['candidate_functions'] / (len(base) ** 2 / 2),
                    'rr_ops_over_pair_ops': rrOps / pairOps if pairOps else None,
                    'sets_equal': rrSols == sols,
                }
                record(row)
                print('d=%d |F|=%d pairs=%d rr_cand=%d ratio=%.3f rr_ops/pair_ops=%.1f %s/%s' % (
                    d, len(base), pairs, rr['candidate_functions'],
                    row['rr_candidates_over_half_f_squared'], row['rr_ops_over_pair_ops'],
                    rr['status'], row['pair_status']), flush=True)

    trials = [r for r in rows if r['kind'] == 'trial']
    usable = [r for r in trials if r['rr_status'] == 'complete' and r['pair_status'] == 'complete']
    series = {}
    for r in usable:
        series.setdefault(r['d'], []).append(r)
    points = []
    for d in sorted(series):
        rr = series[d]
        points.append({
            'd': d,
            'factor_base_points': rr[0]['factor_base_points'],
            'mean_rr_candidates': sum(x['rr_candidates'] for x in rr) / len(rr),
            'mean_pairs': sum(x['pairs_enumerated'] for x in rr) / len(rr),
            'mean_rr_ops': sum(x['rr_ops'] for x in rr) / len(rr),
            'mean_pair_ops': sum(x['pair_ops'] for x in rr) / len(rr),
            'ratio_rr_candidates_over_half_f_squared':
                sum(x['rr_candidates_over_half_f_squared'] for x in rr) / len(rr),
            'ratio_rr_ops_over_pair_ops': sum(x['rr_ops_over_pair_ops'] for x in rr) / len(rr),
        })
    verdict = {}
    if len(points) >= 3:
        fs = [p['factor_base_points'] for p in points]
        verdict = {
            'rr_candidate_exponent_in_F': fitExponent(fs, [p['mean_rr_candidates'] for p in points]),
            'pair_exponent_in_F': fitExponent(fs, [p['mean_pairs'] for p in points]),
            'rr_ops_exponent_in_F': fitExponent(fs, [p['mean_rr_ops'] for p in points]),
            'pair_ops_exponent_in_F': fitExponent(fs, [p['mean_pair_ops'] for p in points]),
            'ratio_first': points[0]['ratio_rr_candidates_over_half_f_squared'],
            'ratio_last': points[-1]['ratio_rr_candidates_over_half_f_squared'],
        }
    (HERE / 'summary.json').write_text(json.dumps(
        {'provenance': provenance, 'points': points, 'verdict': verdict,
         'all_sets_equal': all(r['sets_equal'] for r in usable),
         'limits': contract['limits']}, indent=2) + '\n')
    print(json.dumps(verdict, indent=2), flush=True)


if __name__ == '__main__':
    main()
