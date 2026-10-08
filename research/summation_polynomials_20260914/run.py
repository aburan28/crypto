#!/usr/bin/env python3
"""Frozen ablations of symbolic setup and S4 specialization, with exact replay."""
import argparse
import hashlib
import json
import math
import os
import random
import statistics
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
NAMES = {'micro': 'summation_polynomial_bench', 'stage': 'gaudry_allocation_stage',
         'full_dlp': 'gaudry_cubic_bench'}
TIMING = {'wall_s', 'wall_ms', 'linear_algebra_ms', 'setup_s', 'specialization_s'}

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def scrub(value, varying=()):
    if isinstance(value, dict):
        return {k: scrub(v, varying) for k, v in value.items() if k not in TIMING and k not in varying}
    if isinstance(value, list):
        return [scrub(v, varying) for v in value]
    return value

def interval(pairs, rng):
    groups = sorted({(r['p'], r['seed']) for r in pairs})
    logs = [statistics.mean(math.log(r['ratio']) for r in pairs if (r['p'], r['seed']) == g)
            for g in groups]
    draws = sorted(math.exp(statistics.mean(rng.choices(logs, k=len(logs)))) for _ in range(10000))
    return math.exp(statistics.mean(logs)), [draws[249], draws[9749]]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--build-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    contract = json.loads((HERE/'contract.json').read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    cpu = min(os.sched_getaffinity(0))
    variants = contract['variants']
    bins = {v: args.build_root/'bin'/v for v in variants}
    provenance = {
        'cpu_affinity': [cpu], 'contract_sha256': digest(HERE/'contract.json'),
        'runner_sha256': digest(Path(__file__)),
        'source_sha256': {v: digest(HERE/'sources'/f'{v}.rs') for v in variants},
        'binaries': {v: {n: digest(bins[v]/n) for n in NAMES.values()} for v in variants},
        'harness_sha256': {n: digest(args.build_root/'examples'/f'{n}.rs') for n in NAMES.values()},
        'cargo_lock_sha256': digest(args.build_root/'Cargo.lock'),
        'rustc': subprocess.check_output(['rustc', '--version'], text=True).strip(),
        'cpu': subprocess.check_output(['lscpu'], text=True),
    }
    (args.output/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (args.output/'contract.json').write_text(json.dumps(contract, indent=2)+'\n')
    (args.output/'cargo-lock.txt').write_bytes((args.build_root/'Cargo.lock').read_bytes())
    (args.output/'cargo-manifest.toml').write_bytes((args.build_root/'Cargo.toml').read_bytes())
    schedule = [(kind, p, seed, rep) for kind in NAMES for p, seed in contract[kind+'_cases']
                for rep in range(contract['repetitions'])]
    rng = random.Random(contract['timing_order_seed'])
    rng.shuffle(schedule)
    records = []
    with (args.output/'raw.jsonl').open('w') as raw:
        for kind, p, seed, rep in schedule:
            order = list(variants)
            rng.shuffle(order)
            for variant in order:
                name = f'{kind}-p{p}-seed{seed}-r{rep}-{variant}.json'
                output = args.output/name
                cmd = ['taskset', '-c', str(cpu), str(bins[variant]/NAMES[kind])]
                if kind == 'micro':
                    cmd += [str(p), str(seed), str(contract['micro_targets']),
                            str(contract['setup_repetitions']), str(contract['evaluation_repetitions']),
                            str(output.resolve())]
                elif kind == 'stage':
                    cmd += [str(p), str(seed), str(contract['stage_residuals']), str(output.resolve())]
                else:
                    cmd += ['--p', str(p), '--seed', str(seed), '--groebner', '--cross-check',
                            '--max-residuals', str(contract['max_residuals']), '--json', str(output.resolve())]
                record = dict(kind=kind, p=p, seed=seed, repetition=rep, variant=variant,
                              result=name, command=cmd)
                start = time.perf_counter()
                try:
                    proc = subprocess.run(cmd, capture_output=True, text=True,
                                          timeout=contract['timeout_seconds'])
                    record.update(returncode=proc.returncode, stdout=proc.stdout, stderr=proc.stderr)
                except subprocess.TimeoutExpired as e:
                    record.update(error='timeout', stdout=str(e.stdout), stderr=str(e.stderr))
                record['process_wall_s'] = time.perf_counter()-start
                if output.exists():
                    record['result_sha256'] = digest(output)
                raw.write(json.dumps(record)+'\n')
                raw.flush()
                records.append(record)
    failures, pairs = [], []
    data = {}
    index = {(r['kind'], r['p'], r['seed'], r['repetition'], r['variant']): r for r in records}
    for key, r in index.items():
        if r.get('returncode') != 0 or 'result_sha256' not in r:
            failures.append(r['result'])
            continue
        try:
            data[key] = json.loads((args.output/r['result']).read_text())
        except ValueError:
            failures.append('invalid JSON: '+r['result'])
    # Exact measured specialization charge, constant across targets for these
    # strategies. Reference performs four target-power products and one per term.
    def charge(p, seed, rep, variant):
        d = data['micro', p, seed, rep, variant]
        denominator = d['targets'] * d['evaluation_repetitions']
        assert d['specialization_fp_muls'] % denominator == 0
        actual = d['specialization_fp_muls'] // denominator
        if variant == 'legacy':
            assert actual % 11 == 0
            return 15 * (actual // 11 - 4)
        return actual
    for key, row in index.items():
        kind, p, seed, rep, variant = key
        base_key = kind, p, seed, rep, 'reference'
        if key not in data or base_key not in data:
            continue
        d, ref = data[key], data[base_key]
        try:
            if kind == 'micro':
                assert scrub(d, {'specialization_fp_muls'}) == scrub(ref, {'specialization_fp_muls'})
                if variant in {'legacy', 'setup', 'fixed', 'reference'}:
                    assert d['specialization_fp_muls'] == ref['specialization_fp_muls']
            else:
                varying = {'fp_muls'} if kind == 'stage' else {
                    'fp_muls', 'oracle_fp_muls', 'total_ops', 's', 'fp_muls_per_pair_test'}
                assert scrub(d, varying) == scrub(ref, varying)
                g = d if kind == 'stage' else d[0]['gaudry']
                b = ref if kind == 'stage' else ref[0]['gaudry']
                stats = g['stats'] if kind == 'stage' else g['solve_stats']
                bstats = b['stats'] if kind == 'stage' else b['solve_stats']
                delta = stats['solves'] * (charge(p, seed, rep, variant) - charge(p, seed, rep, 'reference'))
                assert stats['fp_muls'] - bstats['fp_muls'] == delta
                if kind == 'full_dlp':
                    assert g['correct'] and d[0]['rho']['correct']
                    assert g['cross_check_mismatches'] == 0 and g['cross_checked'] == g['residuals']
                    assert g['oracle_fp_muls'] - b['oracle_fp_muls'] == delta
                    assert math.isclose(g['total_ops']-b['total_ops'], delta/g['fp_muls_per_add'], abs_tol=1e-8)
            metrics = ['setup_s', 'specialization_s'] if kind == 'micro' else ['process_wall_s']
            for metric in metrics:
                bt = ref[metric] if kind == 'micro' else index[base_key][metric]
                ct = d[metric] if kind == 'micro' else row[metric]
                pairs.append(dict(kind=kind, metric=metric, p=p, seed=seed, repetition=rep,
                                  variant=variant, baseline_s=bt, candidate_s=ct, ratio=ct/bt))
        except (AssertionError, KeyError) as e:
            failures.append('comparison failed: '+row['result']+' '+str(e))
    summary = dict(complete=not failures, runs=len(records), failures=failures, pairs=pairs,
                   aggregates=[], primary_total_common_operations=None, full_cost_rho_ratio=None,
                   full_cost_floor_ratio=None)
    if not failures:
        for kind in NAMES:
            for metric in (['setup_s', 'specialization_s'] if kind == 'micro' else ['process_wall_s']):
                for variant in variants:
                    ps = [r for r in pairs if r['kind']==kind and r['metric']==metric and r['variant']==variant]
                    ratio, ci = interval(ps, rng)
                    summary['aggregates'].append(dict(kind=kind, metric=metric, variant=variant,
                        pairs=len(ps), sum_s=sum(p['candidate_s'] for p in ps), ratio=ratio, paired_95=ci))
        # Replay all prior frozen completed outputs with the actual legacy
        # implementation. Accounting and optimization rows are kept separate.
        old = REPO/'research/gaudry_allocation_20260914/run-01'
        for kind in ['stage', 'full_dlp']:
            for p, seed in contract['prior_'+kind+'_cases']:
                for rep in range(contract['repetitions']):
                    previous = json.loads((old/f'{kind}-p{p}-seed{seed}-r{rep}-candidate.json').read_text())
                    if scrub(previous) != scrub(data[kind,p,seed,rep,'legacy']):
                        failures.append(f'prior replay failed: {kind} {p} {seed} {rep}')
        summary['complete'] = not failures
    (args.output/'comparison.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary['aggregates'], indent=2))
    print('failures:', failures)
    return bool(failures)

if __name__ == '__main__':
    raise SystemExit(main())
