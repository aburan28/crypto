"""Bounded matched cold-process comparison with independent certificates.

Run only after the query audit passes and model.json is locked. The native
Python adapter does not establish the Rust tournament's common Ir cost unit.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from audit import Checker, require


def verify(result, parameters):
    require(result['fixture']['parameters'] == parameters, 'changed fixed workload')
    checker = Checker(result['fixture'])
    c = checker.curve
    require(result['fixture']['generator_onb'] == [int(parameters['generator'][k], 16) for k in ('x', 'y')],
            'worker changed generator')
    targets = {t['id']: c.decode(checker.convert([int(t[k], 16) for k in ('x', 'y')]))
               for t in parameters['targets']}
    for relation in result['relations']:
        a, b = int(relation['a']), int(relation['b'])
        require(0 < a < c.r and 0 <= b < c.r, 'invalid relation coefficients')
        target = c.mul(c.g, a)
        if relation['stream'] == 'precompute':
            require(b == 0, 'target-dependent precompute')
        else:
            require(b != 0 and relation['stream'][7:] in targets, 'bad target stream')
            target = c.add(target, c.mul(targets[relation['stream'][7:]], b))
        require(target == c.decode(checker.convert(relation['target'])), 'changed relation query')
        checker.relation(target, relation['witness'], relation['row'])
    if result['logs'] is not None:
        require(len(result['logs']) == checker.columns, 'bad log database size')
        for point, scalar in zip(result['fixture']['representatives'], result['logs']):
            require(0 <= scalar < c.r and c.mul(c.g, scalar) == c.decode(checker.convert(point)),
                    'incorrect base logarithm')
    verified = 0
    for name, solution in result['report']['targets'].items():
        if solution['verified']:
            scalar = int(solution['scalar'])
            require(0 <= scalar < c.r and c.mul(c.g, scalar) == targets[name], 'incorrect final scalar')
            verified += 1
    if result['report']['status'] == 'complete':
        require(verified == len(targets), 'missing completed target')
    return {'verified_targets': verified, 'verified_relations': len(result['relations']),
            'status': result['report']['status'], 'independent_checker': 'passed'}


def main():
    root = Path(__file__).resolve().parent
    model = root / 'model.json'
    locked = hashlib.sha256(model.read_bytes()).hexdigest()
    require((root / 'query-audit-v2.json').is_file(), 'query audit must pass before comparison')
    require(json.loads((root / 'aa.json').read_text())['status'] == 'passed', 'baseline A/A must pass')
    out = root / 'full-comparison'
    out.mkdir(exist_ok=False)
    plan = json.loads((root / 'full-plan.json').read_text())
    records = []
    env = dict(os.environ)
    env.pop('ONB_F2_BACKEND', None)
    with (out / 'raw.jsonl').open('w') as stream:
        for cell in plan:
            parameters = json.loads((root / cell['params']).read_text())
            for repetition in range(3):
                variants = ('baseline', 'selector') if repetition % 2 == 0 else ('selector', 'baseline')
                for variant in variants:
                    require(hashlib.sha256(model.read_bytes()).hexdigest() == locked, 'model changed after locking')
                    label = f"{cell['name']}-{repetition}-{variant}"
                    cmd = [sys.executable, str(root / 'worker.py'), '--source',
                           str(root / ('baseline' if variant == 'baseline' else 'candidate-v2') / 'codegen'),
                           '--params', str(root / cell['params']), '--dir', str(out / label)]
                    if variant == 'selector':
                        cmd += ['--model', str(model)]
                    start = time.perf_counter_ns()
                    try:
                        process = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=120)
                    except subprocess.TimeoutExpired as error:
                        def decoded(value):
                            return value.decode(errors='replace') if isinstance(value, bytes) else value or ''
                        process = subprocess.CompletedProcess(cmd, 124, decoded(error.stdout), decoded(error.stderr) + '\nFull worker timeout at 120 seconds')
                    elapsed = time.perf_counter_ns() - start
                    (out / (label + '.stdout.json')).write_text(process.stdout)
                    (out / (label + '.stderr')).write_text(process.stderr)
                    record = {'cell': cell, 'repetition': repetition, 'variant': variant,
                              'process_elapsed_ns': elapsed, 'exit_code': process.returncode,
                              'model_file_sha256': locked if variant == 'selector' else None}
                    try:
                        result = json.loads(process.stdout)
                        record['audit'] = verify(result, parameters)
                        record['relations_sha256'] = hashlib.sha256(json.dumps(result['relations'], sort_keys=True).encode()).hexdigest()
                        record['report'] = result['report']
                    except Exception as error:
                        record['audit'] = {'status': 'error', 'message': str(error)}
                    stream.write(json.dumps(record, sort_keys=True) + '\n')
                    stream.flush()
                    records.append(record)
                print(cell['name'], repetition + 1, '/ 3', flush=True)
    rows = {}
    for variant in ('baseline', 'selector'):
        selected = [r for r in records if r['variant'] == variant]
        rows[variant] = {'runs': len(selected),
                         'complete': sum(r['audit']['status'] == 'complete' for r in selected),
                         'verified_targets': sum(r['audit'].get('verified_targets', 0) for r in selected),
                         'process_elapsed_ns': sum(r['process_elapsed_ns'] for r in selected),
                         'common_ops': None, 'S': None, 'rho_ratio': None, 'floor_ratio': None}
    aa = all(len({r['relations_sha256'] for r in records if r['cell'] == cell and r['variant'] == 'baseline'
                  and 'relations_sha256' in r}) == 1 for cell in plan)
    summary = {'rows': rows, 'model_file_sha256': locked, 'baseline_AA_relation_replay': aa,
               'performance_promoted': False, 'classification': 'accounting; opt-in engineering feature',
               'cost_scope': 'native process wall time including Python startup, cold setup, model load/inference, failed queries, verification, matrix work, database commits, certificate export, output and close; independent audit excluded',
               'limitations': 'no calibrated full-process instruction count; no rho or floor comparison; tiny-field end-to-end checks are functional, not 131-bit extrapolation'}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
