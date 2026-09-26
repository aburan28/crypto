"""Paired legacy/exclusive toy controls; no performance qualification or promotion."""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[4]
HARNESS = ROOT / 'research/ic_candidate_tournament_20260915'
sys.path.insert(0, str(HARNESS))
from generic_phases import verify_native
from generic_queries import verify_queries
from generic_query_law import verify_query_law
from identity import factor_base_inventory
from oracle import require, verify


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')


def panel():
    source = HARNESS / 'goal_20260924/generic-query-law/run_controls.py'
    spec = importlib.util.spec_from_file_location('query_law_controls', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    jobs = module.panel()
    for a in (0, 1):
        job = copy.deepcopy(jobs[0]['job'])
        job.update(mode='rho', curve_a=a)
        jobs.append(dict(name=f'rho-a{a}', job=job, expected='complete'))
    return jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    worker, out = args.worker.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    jobs = panel()
    write(out / 'inputs.json', jobs)
    names = ['src/cryptanalysis/ic_measurement.rs', 'src/cryptanalysis/koblitz_index_calculus.rs',
             'examples/ic_tournament_worker.rs', 'Cargo.lock']
    names += ['research/ic_candidate_tournament_20260915/' + p for p in (
        'generic_phases.py', 'generic_query_law.py', 'generic_queries.py', 'oracle.py', 'identity.py',
        'goal_20260924/generic-exclusive-phases/run_controls.py',
        'goal_20260924/generic-exclusive-phases/PROTOCOL.md',
        'goal_20260924/generic-query-law/run_controls.py')]
    write(out / 'environment.json', dict(
        host=platform.uname()._asdict(), python=sys.version, repetitions=3,
        timeout_seconds=60, rayon_threads=1, memory_limit=None, cpu_affinity=None,
        rustc=subprocess.check_output(['rustc', '--version'], text=True).strip(),
        source_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        worker_sha256=digest(worker), inputs_sha256=digest(out / 'inputs.json'),
        sources={name: digest(ROOT / name) for name in names},
        scope='uncalibrated whole-mode effects; serial collection and tracing both change'))
    env = {**os.environ, 'RAYON_NUM_THREADS': '1'}
    fixtures, outcomes = {}, []
    with (out / 'worker-raw.jsonl').open('x') as raw:
        for item in jobs:
            template = copy.deepcopy(item['job'])
            key = (template['degree'], template['curve_a'])
            if key not in fixtures:
                fixture_job = dict(template, mode='fixture')
                write(out / f'fixture-n{key[0]}a{key[1]}-job.json', fixture_job)
                cp = subprocess.run([str(worker)], input=json.dumps(fixture_job), text=True,
                                    capture_output=True, timeout=60, env=env)
                (out / f'fixture-n{key[0]}a{key[1]}.json').write_text(cp.stdout)
                (out / f'fixture-n{key[0]}a{key[1]}.stderr').write_text(cp.stderr)
                require(cp.returncode == 0, 'fixture generation failed')
                fixtures[key] = json.loads(cp.stdout)['fixture']
            fixture = fixtures[key]
            template['public_targets'] = fixture['targets']
            for repetition in range(3):
                name = f'{item["name"]}-r{repetition}'
                directory = out / name
                directory.mkdir()
                receipt = dict(name=name, status='FAIL', expected=item['expected'], promotion_eligible=False)
                records = {}
                order = ('legacy', 'exclusive') if repetition % 2 == 0 else ('exclusive', 'legacy')
                receipt['order'] = list(order)
                try:
                    for mode in order:
                        job = dict(template, exclusive_phases=(mode == 'exclusive'))
                        run = directory / mode
                        run.mkdir()
                        write(run / 'job.json', job)
                        started = time.monotonic_ns()
                        try:
                            cp = subprocess.run([str(worker)], input=json.dumps(job), text=True,
                                                capture_output=True, timeout=60, env=env)
                        except subprocess.TimeoutExpired as exc:
                            (run / 'timeout-stdout.bin').write_bytes(exc.stdout or b'')
                            (run / 'timeout-stderr.bin').write_bytes(exc.stderr or b'')
                            raise
                        process_ns = time.monotonic_ns() - started
                        (run / 'stdout.json').write_text(cp.stdout)
                        (run / 'stderr.log').write_text(cp.stderr)
                        report = json.loads(cp.stdout)
                        record = dict(name=name, mode=mode, job=job, report=report,
                                      exit_code=cp.returncode, process_wall_ns=process_ns)
                        raw.write(json.dumps(record, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
                        raw.flush()
                        records[mode] = record
                        require(report['status'] in ('complete', 'incomplete'), 'unexpected worker status')
                        require(cp.returncode == (0 if report['status'] == 'complete' else 2), 'wrong worker exit')
                        if item['expected'] != 'complete-or-incomplete':
                            require(report['status'] == item['expected'], 'completion state changed')
                        record['certificate'] = (verify(report, fixture, expected_mode=job['mode'],
                                                        summands=job['config']['summands'])
                                                 if report['status'] == 'complete' else None)
                        if job['mode'] == 'ic':
                            record['queries'] = verify_queries(report, fixture, job['config']['summands'])
                            record['query_law'] = verify_query_law(report, fixture, job)
                            record['inventory'] = factor_base_inventory(report, fixture)
                        if mode == 'exclusive':
                            receipt['phases'] = verify_native(report, job, process_wall_ns=process_ns)
                    legacy, exclusive = records['legacy'], records['exclusive']
                    require(legacy['report']['status'] == exclusive['report']['status'], 'tracing changed completion')
                    require(legacy['certificate'] == exclusive['certificate'], 'tracing changed certificate')
                    if template['mode'] == 'ic':
                        require(legacy['queries'] == exclusive['queries'], 'tracing changed query/counter history')
                        require(legacy['inventory'] == exclusive['inventory'], 'tracing changed base census')
                    a, b = legacy['report']['online_wall_ns'], exclusive['report']['online_wall_ns']
                    receipt.update(status='PASS', worker_status=exclusive['report']['status'],
                                   online_exclusive_over_legacy=None if a is None or b is None else b / a,
                                   process_exclusive_over_legacy=exclusive['process_wall_ns'] / legacy['process_wall_ns'])
                except Exception as exc:
                    receipt['reason'] = f'{type(exc).__name__}: {exc}'
                receipt['artifacts'] = {str(p.relative_to(directory)): digest(p)
                                        for p in sorted(directory.rglob('*')) if p.is_file()}
                write(directory / 'receipt.json', receipt)
                outcomes.append(receipt)
                print(name, receipt['status'], receipt.get('reason', ''), flush=True)
    ratios = [r['online_exclusive_over_legacy'] for r in outcomes
              if r['status'] == 'PASS' and r['online_exclusive_over_legacy'] is not None]
    summary = dict(schema_version=1, pairs=len(outcomes), passed=sum(r['status'] == 'PASS' for r in outcomes),
                   completed_pairs=sum(r.get('worker_status') == 'complete' for r in outcomes),
                   incomplete_pairs=sum(r.get('worker_status') == 'incomplete' for r in outcomes),
                   diagnostic_online_mode_ratio_median=statistics.median(ratios) if ratios else None,
                   diagnostic_online_mode_ratio_range=[min(ratios), max(ratios)] if ratios else None,
                   ratio_limitation='uncalibrated, one host, serial collection and tracing both change; not pure timer overhead',
                   performance_qualified=False, promotion_eligible=False, comparative_online_ms=None,
                   comparative_cold_ms=None, instruction_cost=None, speedup=None, results=outcomes)
    write(out / 'summary.json', summary)
    require(summary['passed'] == summary['pairs'], 'failed controls retained')


if __name__ == '__main__':
    main()
