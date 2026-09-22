"""Reuse repository components only inside two fixed miniature fixtures."""
import importlib.util
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / 'research/nagao_relations/coefficient_pullback'
sys.path.insert(0, str(PRIOR))
import experiment as stage
import nagaoannihilator as scalar
import curves
import indexcalc
import image_solver
import pullback
import cnf
import decomp
import nagaocompare
import nagaodecomp

spec = importlib.util.spec_from_file_location(
    'bounded_lab_cached_pullback',
    ROOT / 'research/nagao_relations/pullback_incremental/optimized.py')
cached = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cached)

VARIANTS = ('quadratic-image', 'coefficient-pullback', 'cached-pullback',
            's4-symmetric', 'chained-s3', 'rho')


def guard(n, d):
    if (n, d) not in ((5, 3), (9, 4)):
        raise ValueError('This laboratory accepts only fixed 5/9-bit fixtures')


def sat_cell(n, d, target_coords, mode, budget, variant):
    """SAT encoding -> solving -> lifting, without a ground-truth oracle.

    Based on solver_06/run.py, with expected solutions removed from the
    computation. Ground truth is checked separately after timed runs.
    chained-s3 uses the repository chained system from solver_02 and
    nagaocompare (decomp.buildSystem plus decomp.encode).
    """
    import time
    guard(n, d)
    if variant not in ('s4-symmetric', 'chained-s3'):
        raise ValueError('unsupported SAT encoding')
    start = time.perf_counter()
    deadline = start + budget
    f = scalar.CountedField(n)
    curve = scalar.CountedCurve(f)
    target = tuple(f.fromCoords(x) for x in target_coords)
    c = cnf.Cnf()
    if variant == 'chained-s3':
        prog, roots = decomp.buildSystem(n, f.n, 3, min(12, n))
        xr = target_coords[0]
        pvars = decomp.encode(prog, roots, n, 3, n, xr, c, orderPoints=False)
        nagaodecomp.addRestrictedDomain(c, pvars, xr)
        variables = {'pvars': pvars}
    else:
        _, _, variables = stage.previous.s4.s4.encode(n, d, target, c, variant)
    for xs in variables['pvars']:
        for bit in xs[d:]:
            c.addClause([-bit])
    solver = stage.previous.s4.previous.loadSolver(c)
    built = time.perf_counter()
    phases = {'construct_load': built-start, 'solve': 0., 'extract_verify': 0.}
    solutions = set()
    calls = rejected = 0
    complete = False
    first = None
    while time.perf_counter() < deadline:
        before = time.perf_counter()
        ok, model = solver.solve(time_limit=max(.001, deadline-before))
        solved = time.perf_counter()
        phases['solve'] += solved-before
        calls += 1
        if ok is None:
            break
        if not ok:
            complete = True
            break
        vals = {i: v for i, v in enumerate(model) if i}
        if not nagaocompare.cnfSatisfied(c, vals):
            raise ArithmeticError('CNF model failed independent evaluation')
        xs = sorted(sum(int(vals[v]) << i for i, v in enumerate(bits))
                    for bits in variables['pvars'])
        if (any(x == 0 or x >= 1 << d or x == target_coords[0] for x in xs)
                or len(set(xs)) != 3):
            raise ArithmeticError('SAT result violates the common domain')
        f.phase = 'extract_verify'
        valid = indexcalc.liftAndCheck(f, curve, xs, target) is not None
        if valid:
            solutions.add(tuple(xs))
            if first is None:
                first = time.perf_counter()-start
        else:
            rejected += 1
        solver.add_clause([-v if vals[v] else v
                           for bits in variables['pvars'] for v in bits])
        phases['extract_verify'] += time.perf_counter()-solved
        if valid and mode == 'first':
            break
    elapsed = time.perf_counter()-start
    phases['orchestration'] = max(0., elapsed-sum(phases.values()))
    return {
        'variant': variant, 'n': n, 'd': d, 'target': list(target_coords),
        'mode': mode,
        'status': 'first' if solutions and mode == 'first' else
                  ('complete' if complete else 'timeout'),
        'solutions': [list(x) for x in sorted(solutions)],
        'verified_unique_relations': len(solutions), 'calls': calls,
        'rejected_candidates': rejected, 'first_verified_seconds': first,
        'all_phase_seconds': elapsed, 'phase_seconds': phases,
        'field_api_counts': f.report(),
        'field_api_coverage': 'lifting only; symbolic encoding/search uncounted',
        'solver_operations': None, 'constructed_cnf_totals': c.stats(),
        'full_dlp_S': None, 'rho_ratio': None,
    }


def adapter(variant):
    if variant not in VARIANTS or variant == 'rho':
        raise ValueError('unsupported relation adapter')

    def cell(n, d, target, mode, budget):
        guard(n, d)
        if variant == 'quadratic-image':
            return image_solver.cell(n, d, target, mode, budget)
        if variant == 'coefficient-pullback':
            return pullback.cell(n, d, target, mode, budget)
        if variant == 'cached-pullback':
            return cached.cell(n, d, target, mode, budget, variant)
        return sat_cell(n, d, target, mode, budget, variant)
    return cell
