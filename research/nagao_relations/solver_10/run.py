"""Matched decomposition-solver panel on the real ECC2K-130 curve.

Every solver variant is imported unmodified from its own frozen folder; only
this driver is new.  The toy campaigns (solver_02..09) could lean on an
exhaustive pair oracle over the whole group.  At n = 131 no such oracle exists,
so the two oracle-equality assertions inside the reused SAT drivers are stood
down by an explicitly permissive `expected` object, and correctness is carried
instead by independent re-verification here: every reported abscissa triple is
lifted with pointFromX, re-checked for subspace membership and the exclusions,
and its signed group sum is recomputed and compared with the target.

No relation yield, full-DLP S or rho ratio is computed from these runs; see
contract.json "limits" for why they would be meaningless at d << 45.
"""
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
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CODE = ROOT / 'ecc2k130/codegen'
sys.path.insert(0, str(CODE))
for folder in ('solver_02', 'solver_04', 'solver_06', 'solver_07', 'solver_08'):
    sys.path.insert(0, str(HERE.parent / folder))
import curves
import field
import image_solver
import nagaoannihilator as scalar
import optimized

# Published ECC2K-130 group order: #E = 4r with r a 129-bit prime.
CHALLENGE_R = 680564733841876926932320129493409985129
SEED = 20260916


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


campaign06 = load('campaign06', HERE.parent / 'solver_06/run.py')
campaign04 = load('campaign04', HERE.parent / 'solver_04/run.py')


class AnyExpected:
    """Stands in for the exhaustive oracle that does not exist at n = 131.

    Admits every candidate the reused driver would have checked against the
    oracle, and compares equal to any solution set, so a completed branch is
    reported rather than raising.  It removes checks; it never adds one.  The
    checks it removes are replaced by verify() below.
    """

    def __contains__(self, item):
        return True

    def __len__(self):
        return 0

    def __eq__(self, other):
        return True

    def __ne__(self, other):
        return False

    __hash__ = None


def verify(f, curve, d, target, xs):
    """Re-derive the relation from the abscissas alone, ignoring the solver."""
    if len(set(xs)) != 3:
        raise ArithmeticError('repeated abscissa')
    if any(x == 0 or x >= 1 << d or x == f.toCoords(target[0]) for x in xs):
        raise ArithmeticError('abscissa outside the factor base or excluded')
    total = None
    points = []
    for x in xs:
        p = curve.pointFromX(f.fromCoords(x))
        if p is None:
            raise ArithmeticError('abscissa is not on the curve')
        # Both lifts of x differ by the curve negation; accept the one that
        # makes the signed sum close, which is the relation's own choice.
        points.append(p)
    for signs in range(8):
        total = None
        for i, p in enumerate(points):
            total = curve.add(total, curve.neg(p) if signs >> i & 1 else p)
        if total == target:
            return True
    raise ArithmeticError('signed sum does not reach the target')


def buildPanel(n, d, rng):
    f = scalar.CountedField(n)
    curve = curves.Curve(f)
    order = curves.curveOrder(n)
    base = []
    for i in range(1, 1 << d):
        p = curve.pointFromX(f.fromCoords(i))
        if p is not None:
            if not curve.onCurve(p):
                raise ArithmeticError('factor-base point off curve')
            base.append((i, p))
    uniform = []
    while len(uniform) < 4:
        p = curve.pointFromX(f.randomElement(rng))
        if p is None or f.toCoords(p[0]) < 1 << d:
            continue
        uniform.append((p, None))
    planted = []
    seen = set()
    while len(planted) < 4:
        triple = tuple(sorted(x for x, _ in rng.sample(base, 3)))
        if triple in seen:
            continue
        total = None
        for x in triple:
            total = curve.add(total, dict(base)[x])
        if total is None or f.toCoords(total[0]) < 1 << d:
            continue
        seen.add(triple)
        verify(f, curve, d, total, triple)
        planted.append((total, triple))
    return f, curve, order, base, {'uniform': uniform, 'known_decomposable': planted}


def runVariant(f, curve, d, target, variant, mode, budget):
    coords = [f.toCoords(target[0]), f.toCoords(target[1])]
    if variant == 'quadratic-image':
        return image_solver.cell(f.m, d, coords, mode, budget)
    if variant == 'quadratic-optimized':
        return optimized.cell(f.m, d, coords, mode, budget)
    driver = campaign04 if variant == 'chained-s3' else campaign06
    return driver.runCell(f, curve, d, target, AnyExpected(), variant, mode, budget)


def main():
    raw = HERE / 'raw.jsonl'
    if raw.exists():
        raise SystemExit('Evidence exists; choose a new sibling folder')
    contract = json.loads((HERE / 'contract.json').read_text())
    sources = list(CODE.glob('*.py')) + [HERE / 'run.py', HERE / 'contract.json']
    for folder in ('solver_02', 'solver_04', 'solver_05', 'solver_06', 'solver_07', 'solver_08'):
        sources += list((HERE.parent / folder).glob('*.py'))
    provenance = {
        'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'command': [sys.executable] + sys.argv,
        'python': platform.python_version(),
        'pycryptosat': importlib.metadata.version('pycryptosat'),
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
            rng = random.Random(SEED + n * 1000 + d)
            t = time.perf_counter()
            try:
                f, curve, order, base, strata = buildPanel(n, d, rng)
            except Exception:
                record({'kind': 'validation_failure', 'n': n, 'd': d, 'traceback': traceback.format_exc()})
                raise
            if n == 131 and order != 4 * CHALLENGE_R:
                record({'kind': 'validation_failure', 'n': n, 'd': d, 'reason': 'group order is not 4r'})
                raise SystemExit('group order mismatch: not the ECC2K-130 curve')
            record({'kind': 'panel', 'n': n, 'd': d, 'group_order': order,
                    'challenge_order_checked': n == 131,
                    'factor_base_abscissas': [x for x, _ in base],
                    'factor_base_size': len(base),
                    'targets': {st: [{'target': [f.toCoords(p[0]), f.toCoords(p[1])], 'planted': list(tri) if tri else None}
                                     for p, tri in rowset] for st, rowset in strata.items()},
                    'setup_seconds': time.perf_counter() - t})
            print('panel n=%d d=%d base=%d' % (n, d, len(base)), flush=True)
            for stratum, rowset in strata.items():
                for target, planted in rowset:
                    for mode in contract['modes']:
                        for variant in contract['variants']:
                            try:
                                f.counts = {}
                                f.phase = 'setup'
                                r = runVariant(f, curve, d, target, variant, mode, contract['seconds_per_instance'])
                                found = {tuple(x) for x in r['solutions']}
                                for xs in found:
                                    verify(f, curve, d, target, xs)
                                if planted is not None and r['status'] == 'complete' and tuple(planted) not in found:
                                    raise ArithmeticError('complete run missed its planted triple')
                                r['within_budget'] = r['all_phase_seconds'] <= contract['seconds_per_instance']
                                record({'kind': 'trial', 'stratum': stratum, 'planted': list(planted) if planted else None,
                                        'independently_verified_relations': len(found), **r})
                                print(n, d, stratum, mode, variant, r['status'], len(found),
                                      round(r['all_phase_seconds'], 3), flush=True)
                            except Exception:
                                record({'kind': 'error', 'n': n, 'd': d, 'stratum': stratum, 'mode': mode,
                                        'variant': variant, 'target': [f.toCoords(target[0]), f.toCoords(target[1])],
                                        'traceback': traceback.format_exc()})
                                print(n, d, stratum, mode, variant, 'ERROR', flush=True)
    groups = []
    for panel in contract['panels']:
        for variant in contract['variants']:
            for mode in contract['modes']:
                for st in ('uniform', 'known_decomposable'):
                    rr = [r for r in rows if r['kind'] == 'trial'
                          and (r['n'], r['d'], r['variant'], r['mode'], r['stratum']) == (panel['n'], panel['d'], variant, mode, st)]
                    seconds = sum(r['all_phase_seconds'] for r in rr)
                    relations = sum(r['independently_verified_relations'] for r in rr)
                    ops = [r['field_api_counts']['totals']['fieldOperations'] for r in rr if r.get('field_api_counts')]
                    groups.append({'n': panel['n'], 'd': panel['d'], 'variant': variant, 'mode': mode, 'stratum': st,
                                   'attempted': len(rr),
                                   'resolved': sum(r['status'] != 'timeout' for r in rr),
                                   'resolved_within_budget': sum(r['status'] != 'timeout' and r['within_budget'] for r in rr),
                                   'verified_relations': relations,
                                   'all_phase_seconds': seconds,
                                   'diagnostic_seconds_per_relation': seconds / relations if relations else None,
                                   'field_api_operations': sum(ops) if ops else None,
                                   'full_dlp_S': None, 'rho_ratio': None})
    (HERE / 'summary.json').write_text(json.dumps(
        {'provenance': provenance, 'groups': groups,
         'errors': [r for r in rows if r['kind'] == 'error'],
         'limits': contract['limits']}, indent=2) + '\n')


if __name__ == '__main__':
    main()
