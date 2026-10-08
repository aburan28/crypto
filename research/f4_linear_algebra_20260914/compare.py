#!/usr/bin/env python3
"""Compare frozen pairs; missing runs or changed mathematical answers fail."""
import itertools
import json
import math
import pathlib
import random
import statistics as stats
import sys

folder = pathlib.Path(sys.argv[1])
contract = json.loads((folder / 'contract.json').read_text())
processes = {x['name']: x for x in json.loads((folder / 'processes.json').read_text())}
errors = []
hashes = json.loads((folder / 'hashes.json').read_text())
for binary in ['ic', 'f4_linear_algebra_bench']:
    if hashes['reference-' + binary] == hashes['candidate-' + binary]:
        errors.append('identical reference/candidate binary: ' + binary)


def load(name):
    text = (folder / (name + '.stdout')).read_text() if (folder / (name + '.stdout')).exists() else ''
    try:
        if name.startswith('dlp-'):
            return json.loads(text)
        return [json.loads(line[line.index('{"'):]) for line in text.splitlines() if '{"' in line]
    except (ValueError, TypeError):
        errors.append('invalid output: ' + name)
        return {} if name.startswith('dlp-') else []


def finished(name, allowed=(0,)):
    proc = processes.get(name, {})
    return proc.get('status') == 'finished' and proc.get('returncode') in allowed


def ratio_summary(pairs):
    if not pairs or any(a <= 0 or b <= 0 for a, b in pairs):
        return {'candidate_over_reference': None, 'paired_95pct_ci': None}
    logs = [math.log(b / a) for a, b in pairs]
    rng = random.Random(20260914)
    boot = sorted(math.exp(stats.mean(rng.choices(logs, k=len(logs)))) for _ in range(4000))
    return {'candidate_over_reference': math.exp(stats.mean(logs)),
            'paired_95pct_ci': [boot[100], boot[3899]]}


kernel_raw = load('kernel')
kernel = []
for case, seed in sorted({(r['case'], r['seed']) for r in kernel_raw}):
    samples = [r for r in kernel_raw if (r['case'], r['seed']) == (case, seed)]
    if len(samples) != contract['kernel']['paired_repetitions'] or not all(r['correct'] for r in samples):
        errors.append('incomplete/incorrect kernel case: ' + case)
    first = samples[0]
    if case.startswith('dense-') and first['rank'] < min(first['rows'], first['cols']) - 4:
        errors.append('dense generator produced unexpectedly low rank: ' + case)
    if any(any(r[k] != first[k] for k in ['input_hash', 'rank', 'reference_xors', 'candidate_xors']) for r in samples):
        errors.append('unstable kernel input/result: ' + case)
    kernel.append({'case': case, 'seed': seed, 'split': 'development' if seed == 17 else 'fresh',
        'rows': first['rows'], 'cols': first['cols'], 'rank': first['rank'], 'correct': all(r['correct'] for r in samples),
        'reference_word_xors': first['reference_xors'], 'candidate_word_xors': first['candidate_xors'],
        'word_xor_ratio': first['candidate_xors'] / first['reference_xors'],
        'reference_ns': stats.median(r['reference_ns'] for r in samples),
        'candidate_ns': stats.median(r['candidate_ns'] for r in samples),
        **ratio_summary([(r['reference_ns'], r['candidate_ns']) for r in samples])})
expected_kernel = len(contract['kernel']['seeds']) * (len(contract['kernel']['synthetic']) + len(contract['kernel']['macaulay']))
if len(kernel) != expected_kernel or not finished('kernel'):
    errors.append('kernel suite did not finish all cases')

stages = []
target_checks = 0
for (n, ell, m), split in itertools.product(contract['stages']['cases'], contract['stages']['splits']):
    pairs, rows = [], {'reference': [], 'candidate': []}
    for rep in range(contract['stages']['paired_repetitions']):
        pair = {}
        target_rows = {}
        for variant in rows:
            name = f'stage-n{n}-l{ell}-m{m}-{split}-r{rep}-{variant}'
            raw = load(name)
            targets = [r for r in raw if r['phase'] == 'target']
            expected = (1 << n) if n <= 5 else 32
            if not finished(name) or len(targets) != expected or len({r['target'] for r in targets}) != expected or not all(r['correct'] for r in targets):
                errors.append('incomplete/incorrect stage: ' + name)
                continue
            setup = next(r['field_setup_ns'] for r in raw if r['phase'] == 'setup')
            measured = {'cold_stage_ns': setup + sum(r['encoding_ns'] + r['solver_ns'] for r in targets),
                        'solver_ns': sum(r['solver_ns'] for r in targets),
                        'root_f4_ns': sum(r['root_f4_ns'] for r in targets),
                        'root_word_xors': sum(r['root_word_xors'] for r in targets),
                        'oracle_verification_ns': sum(r['verification_ns'] for r in targets)}
            rows[variant].append(measured)
            pair[variant] = measured['cold_stage_ns']
            target_rows[variant] = targets
            target_checks += len(targets)
        if len(pair) == 2:
            pairs.append((pair['reference'], pair['candidate']))
            for a, b in zip(target_rows['reference'], target_rows['candidate']):
                if any(a[k] != b[k] for k in ['target', 'input_hash', 'output_hash', 'solutions', 'reductions', 'exhausted']):
                    errors.append(f'changed solver answer/work: {n},{ell},{m},{split},{rep}')
    stages.append({'n': n, 'ell': ell, 'm': m, 'split': split,
                   'variants': {v: {k: stats.median(s[k] for s in samples) for k in samples[0]} if samples else {}
                                for v, samples in rows.items()}, **ratio_summary(pairs)})

dlp = []
for n, solver, (seed, log, split), rep in itertools.product(contract['dlp']['degrees'], contract['dlp']['solvers'],
        contract['dlp']['fixtures'], range(contract['dlp']['paired_repetitions'])):
    base = f'dlp-n{n}-{solver}-s{seed}-k{log}-{split}-r{rep}'
    a, b = load(base + '-reference'), load(base + '-candidate')
    matched = bool(a) and bool(b) and all(finished(base + '-' + v, (0, 1)) for v in ['reference', 'candidate'])
    matched &= all(a.get(k) == b.get(k) for k in ['arguments', 'parameters', 'mode', 'factor_base', 'result', 'counts', 'status'])
    if not matched:
        errors.append('changed/incomplete process: ' + base)
    for variant, r in [('reference', a), ('candidate', b)]:
        verified = r.get('result', {}).get('verified', False)
        if finished(base + '-' + variant) and not verified:
            errors.append('successful exit without verified scalar: ' + base + '-' + variant)
        dlp.append({'n': n, 'solver': solver, 'seed': seed, 'known_log': log, 'split': split, 'rep': rep,
                    'variant': variant, 'matches_reference': matched, 'status': r.get('status', 'missing'),
                    'verified': verified, 'elapsed_seconds': r.get('elapsed_seconds'),
                    'phase_costs': r.get('timing_seconds', {}), 'resources': r.get('resources', {})})

dlp_groups = []
for n, solver, split in itertools.product(contract['dlp']['degrees'], contract['dlp']['solvers'], contract['stages']['splits']):
    rows = [r for r in dlp if (r['n'], r['solver'], r['split']) == (n, solver, split)]
    variants = {v: [r for r in rows if r['variant'] == v] for v in ['reference', 'candidate']}
    pairs = [(a['elapsed_seconds'], b['elapsed_seconds']) for a, b in zip(*variants.values())
             if a['elapsed_seconds'] is not None and b['elapsed_seconds'] is not None]
    dlp_groups.append({'n': n, 'solver': solver, 'split': split,
        'variants': {v: {'executions': len(rs), 'verified': sum(r['verified'] for r in rs),
                        'sum_elapsed_seconds': sum(r['elapsed_seconds'] or 0 for r in rs)} for v, rs in variants.items()},
        **ratio_summary(pairs)})

expected_processes = 1 + 2 * len(stages) * contract['stages']['paired_repetitions'] + len(dlp)
if len(processes) != expected_processes:
    errors.append('missing process records')
kernel_geomean = math.exp(stats.mean(math.log(r['candidate_over_reference']) for r in kernel)) if kernel else None
kernel_target_met = kernel_geomean is not None and kernel_geomean <= 1.1 and any(
    r['split'] == 'fresh' and r['cols'] > 64 and r['candidate_ns'] <= .9 * r['reference_ns'] and r['word_xor_ratio'] < 1
    for r in kernel)
summary = {'classification': 'engineering; stage diagnostic', 'errors': errors,
    'correct': not errors, 'kernel_target_met': kernel_target_met,
    'expected_processes': expected_processes, 'recorded_processes': len(processes),
    'stage_target_checks': target_checks, 'kernel_geomean_candidate_over_reference': kernel_geomean,
    'kernel_rows': kernel, 'stage_rows': stages, 'dlp_rows': dlp, 'dlp_groups': dlp_groups,
    'total_calibrated_operations': None, 'S': None, 'rho_ratio': None, 'floor_ratio': None, 'exponent': None,
    'limitations': ['Word XORs cover only elimination, not total attack operations.',
                   'Cold stage includes field setup, encoding and solving; independent exhaustive oracle and extra root reduction are separately measured validation.',
                   'Shared-host runtime and paired bootstrap intervals are diagnostics. Three repetitions per DLP fixture do not establish independent workload sampling.',
                   'Incomplete attempts retained; no full-DLP speedup or asymptotic claim.']}
(folder / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps({k: v for k, v in summary.items() if not k.endswith('_rows') and k != 'dlp_groups'}, indent=2))
if errors:
    raise SystemExit(1)
