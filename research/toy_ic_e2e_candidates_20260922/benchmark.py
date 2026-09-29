"""Reproducible, fixed-size end-to-end laboratory. No target input interface."""
import argparse
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
import traceback

IMPORT_STARTED = time.perf_counter()
import bindings as b
import pipeline
IMPORT_SECONDS = time.perf_counter()-IMPORT_STARTED
HERE = Path(__file__).resolve().parent


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def validate():
    oracle = b.stage.previous.PairOracle(5, 3)
    f, curve = oracle.f, oracle.curve
    truth, base_size = b.stage.previous.old.baseline.oracle(f, curve, 3)
    targets = b.nagaocompare.affinePoints(f, curve)
    for target in targets:
        if truth.get(target, set()) != oracle.expected(target):
            raise ArithmeticError('independent pair/triple oracle disagreement')
    checked = 0
    for variant in b.VARIANTS:
        if variant == 'rho':
            continue
        for target in targets:
            result = b.adapter(variant)(5, 3, [f.toCoords(x) for x in target],
                                        'enumerate', 3)
            got = {tuple(x) for x in result['solutions']}
            if result['status'] != 'complete' or got != oracle.expected(target):
                raise ArithmeticError((variant, target, 'exhaustive mismatch'))
            checked += 1
        print('Exhaustive five-bit validation:', variant, '43/43', flush=True)
    for n, d in [(131, 44), (11, 5), (5, 2)]:
        try:
            b.guard(n, d)
        except ValueError:
            continue
        raise ArithmeticError('toy fixture guard failed')
    # Exercise the inconsistent-row rejection independently of successful runs.
    matrix = pipeline.Rows(1, pipeline.ScalarOps(11))
    matrix.add([1, 1])
    try:
        matrix.add([1, 2])
    except ArithmeticError:
        pass
    else:
        raise ArithmeticError('inconsistent matrix accepted')
    return {'exhaustive_targets': len(targets), 'adapter_target_checks': checked,
            'rational_base_abscissae': base_size,
            'pair_and_triple_truth_agree': True,
            'out_of_scope_parameters_rejected': True,
            'inconsistent_relation_rejected': True}


def audit(rows):
    oracles = {}
    seen = set()
    calls = complete_empty = relation_results = 0
    boundaries = []
    for n, d, full in [(5, 3, 44), (9, 4, 508)]:
        oracle = b.stage.previous.PairOracle(n, d)
        oracles[n] = oracle
        points = b.nagaocompare.affinePoints(oracle.f, oracle.curve)
        covered = sum(bool(oracle.expected(p)) for p in points)
        m = len(oracle.base)
        boundaries.append({'n': n, 'd': d, 'rational_abscissae': m,
                           'affine_targets': full-1, 'covered_targets': covered,
                           'exact_coverage': covered/(full-1),
                           'coverage_upper_bound': min(1., 8*math.comb(m, 3)/(full-1))})
    for row in rows:
        if 'target' in row:
            oracle = oracles[row['panel']['n']]
            f, curve = oracle.f, oracle.curve
            G = tuple(f.fromCoords(v) for v in row['generator'])
            Q = tuple(f.fromCoords(v) for v in row['target'])
            if row['status'] == 'verified':
                if curve.mul(G, row['recovered_scalar']) != Q:
                    raise ArithmeticError('independent final verification failed')
        for attempt in row.get('attempts', []):
            result = attempt['result']
            oracle = oracles[result['n']]
            target = tuple(oracle.f.fromCoords(v) for v in result['target'])
            expected = oracle.expected(target)
            got = {tuple(x) for x in result['solutions']}
            if not got <= expected:
                raise ArithmeticError('full-run adapter returned false relation')
            if result['status'] == 'complete' and got != expected:
                raise ArithmeticError('full-run adapter returned false empty set')
            calls += 1
            complete_empty += result['status'] == 'complete' and not got
            relation_results += bool(got)
            seen.add((result['n'], tuple(result['target'])))
    historical_path = b.PRIOR/'e2e_results/raw.jsonl'
    old = [json.loads(x) for x in historical_path.read_text().splitlines()]
    old = {(r['panel']['n'], r['seed'], r['variant']): r for r in old}
    replays = 0
    for row in rows:
        key = (row['panel']['n'], row['seed'], row['variant'])
        if row.get('repetition') != 0 or key not in old:
            continue
        for field in ('status', 'recovered_scalar', 'all_phase_field_api_counts',
                      'scalar_modular_counts'):
            if row.get(field) != old[key].get(field):
                raise ArithmeticError(('historical replay differs', key, field))
        replays += 1
    return {'audited_calls': calls, 'distinct_oracle_targets': len(seen),
            'complete_empty_answers_checked': complete_empty,
            'nonempty_relation_results_checked': relation_results,
            'historical_counter_replays': replays, 'boundaries': boundaries}


def paired_interval(ratios):
    if not ratios:
        return None
    # Ratios are per-seed medians, so repetitions are not independent samples.
    logs = [math.log(x) for x in ratios]
    rng = random.Random(220926)
    boot = sorted(math.exp(statistics.mean(rng.choices(logs, k=len(logs))))
                  for _ in range(2000))
    return {'geomean': math.exp(statistics.mean(logs)),
            'seed_cluster_bootstrap_95pct': [boot[49], boot[1949]],
            'seed_clusters': len(logs)}


def summarize(rows, contract, audit_result):
    groups = []
    for panel in contract['panels']:
        n = panel['n']
        for variant in contract['variants']:
            selected = [r for r in rows if r['panel']['n'] == n and r['variant'] == variant]
            verified = [r for r in selected if r['status'] == 'verified']
            paired = {}
            for reference in ('coefficient-pullback', 'rho'):
                ratios = []
                for seed in sorted({r['seed'] for r in selected}):
                    a = [r for r in selected if r['seed'] == seed]
                    z = [r for r in rows if r['panel']['n'] == n and r['variant'] == reference and r['seed'] == seed]
                    if (not a or len(a) != len(z) or
                            any(r['status'] != 'verified' for r in a+z)):
                        continue
                    ratios.append(statistics.median(r['cold_wall_seconds'] for r in a) /
                                  statistics.median(r['cold_wall_seconds'] for r in z))
                paired[reference] = paired_interval(ratios)
            phases = sorted({k for r in selected for k in r.get('phase_seconds', {})})
            groups.append({
                'n': n, 'subgroup_order': panel['subgroup_order'], 'variant': variant,
                'attempted_runs': len(selected), 'verified_runs': len(verified),
                'failed_runs': len(selected)-len(verified),
                'distinct_final_targets': len({tuple(r['target']) for r in selected if 'target' in r}),
                'median_cold_seconds': statistics.median(r['cold_wall_seconds'] for r in selected) if selected else None,
                'summed_cold_seconds': sum(r.get('cold_wall_seconds', 0) for r in selected),
                'phase_seconds_summed': {k: sum(r.get('phase_seconds', {}).get(k, 0) for r in selected) for k in phases},
                'all_phase_field_api_counts_summed': {k: sum(r.get('all_phase_field_api_counts', {}).get(k, 0) for r in selected)
                                                       for k in ('additions', 'multiplications', 'squarings')},
                'scalar_modular_counts_summed': {k: sum(r.get('scalar_modular_counts', {}).get(k, 0) for r in selected)
                                                  for k in ('additions', 'multiplications', 'inversions')},
                'decomposition_attempts': sum(len(r.get('attempts', [])) for r in selected),
                'decomposition_timeouts': sum(a['result']['status'] == 'timeout' for r in selected for a in r.get('attempts', [])),
                'field_api_coverage': selected[0].get('field_api_coverage') if selected else None,
                'diagnostic_wall_time_over_reference': paired,
                'full_dlp_S': None, 'calibrated_cost_over_rho': None,
                'calibrated_cost_over_floor': None,
                'classification': 'accounting/integration; no common-unit speedup claim',
            })
    return {'contract': contract, 'groups': groups, 'audit': audit_result,
            'any_unverified_run': any(r['status'] != 'verified' for r in rows),
            'process_peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'import_seconds_excluded_from_per_fixture_cold_times': IMPORT_SECONDS,
            'limits': ['Fixed toy fields only', 'No ECC2K-130 runtime estimate',
                       'Shared exhaustive curve setup dominates some small runs',
                       'Plain rho control; no automorphism-aware rho',
                       'SAT search work is not converted to field operations',
                       'No F4 backend measured', 'Fresh seeds may repeat tiny-group targets']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    contract = json.loads((HERE/'contract.json').read_text())
    if args.smoke:
        contract.update(frozen_seeds=[91], fresh_seeds=[], repetitions=1)
    sources = set([HERE/'contract.json', Path(__file__)])
    for module in tuple(sys.modules.values()):
        value = getattr(module, '__file__', None)
        if value:
            path = Path(value).resolve()
            if path.suffix == '.py' and path.is_relative_to(b.ROOT):
                sources.add(path)
    provenance = {
        'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=b.ROOT, text=True).strip(),
        'source_sha256': {str(p.relative_to(b.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(sources)},
        'python': platform.python_version(), 'platform': platform.platform(),
        'cpu_model': next((line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')), 'unknown'),
        'pycryptosat': importlib.metadata.version('pycryptosat'),
        'threads_per_solver': 1, 'command': sys.argv, 'contract': contract,
    }
    write_json(args.output/'provenance.json', provenance)
    checks = validate()
    write_json(args.output/'validation.json', checks)
    if args.validate_only:
        return
    rows = []
    seeds = contract['frozen_seeds']+contract['fresh_seeds']
    with (args.output/'raw.jsonl').open('x') as stream:
        for repetition in range(contract['repetitions']):
            for panel in contract['panels']:
                for seed in seeds:
                    variants = contract['variants'][:]
                    random.Random(seed+repetition*10000+panel['n']).shuffle(variants)
                    for variant in variants:
                        started = time.perf_counter()
                        try:
                            row = pipeline.run(panel, seed, variant, contract,
                                               None if variant == 'rho' else b.adapter(variant))
                        except Exception:
                            row = {'panel': panel, 'seed': seed, 'variant': variant,
                                   'status': 'error', 'traceback': traceback.format_exc(),
                                   'cold_wall_seconds': time.perf_counter()-started}
                        row.update(repetition=repetition,
                                   corpus='frozen' if seed in contract['frozen_seeds'] else 'fresh')
                        rows.append(row)
                        stream.write(json.dumps(row)+'\n')
                        stream.flush()
                        print(panel['n'], seed, repetition, variant, row['status'],
                              round(row['cold_wall_seconds'], 4), flush=True)
    before = time.perf_counter()
    audited = audit(rows)
    audited['offline_audit_seconds'] = time.perf_counter()-before
    result = summarize(rows, contract, audited)
    write_json(args.output/'summary.json', result)
    print(json.dumps({'verified': sum(r['status'] == 'verified' for r in rows),
                      'runs': len(rows), 'audit': audited}), flush=True)
    if result['any_unverified_run']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
