"""Frozen toy query-law controls. Preserve failures; never rank performance."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[4]
HARNESS = ROOT / 'research/ic_candidate_tournament_20260915'
sys.path.insert(0, str(HARNESS))
from generic_query_law import verify_query_law, verify_rust_vectors
from identity import factor_base_inventory
from oracle import require, verify


def write(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def panel():
    jobs = []
    solvers = ('pair_table', 'enumerate', 'f4', 'f5', 'inherited_f4', 'sat_xor', 'sat_cnf')
    template = dict(mode='ic', degree=9, curve_a=0, target_seeds=[2026092556],
                    algorithm_seed=2026092556, factor_base=dict(kind='factor', index=0),
                    config=dict(solver='pair_table', linear_algebra='dense', summands=2,
                                batch_trials=8, max_trials=256))
    for a in (0, 1):
        for solver in solvers:
            for la in ('dense', 'sparse'):
                job = copy.deepcopy(template)
                job['curve_a'] = a
                job['config'].update(solver=solver, linear_algebra=la)
                jobs.append(dict(name=f'a{a}-{solver}-{la}', job=job, expected='complete'))
    for solver in solvers:
        job = copy.deepcopy(template)
        job['config'].update(solver=solver, batch_trials=1, max_trials=1)
        jobs.append(dict(name=f'incomplete-{solver}', job=job, expected='incomplete'))
    for degree, a in ((9, 0), (13, 0)):
        for window in (0, 1, 1000):
            for la in ('dense', 'sparse'):
                job = copy.deepcopy(template)
                job.update(degree=degree, curve_a=a,
                           factor_base=dict(kind='subgroup_orbits', seed=43,
                                            points=36 if degree == 9 else 78))
                job['config'].update(summands=3, collection_window=window, linear_algebra=la,
                                     batch_trials=128, max_trials=256)
                jobs.append(dict(name=f'n{degree}-window{window}-{la}', job=job,
                                 expected='complete-or-incomplete'))
    return jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', type=Path, required=True)
    parser.add_argument('--vectors', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    worker, out = args.worker.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    jobs = panel()
    write(out / 'inputs.json', jobs)
    names = ('examples/ic_tournament_worker.rs', 'examples/ic_query_law_vectors.rs',
             'src/cryptanalysis/koblitz_index_calculus.rs', 'Cargo.lock',
             'research/ic_candidate_tournament_20260915/generic_query_law.py',
             'research/ic_candidate_tournament_20260915/generic_queries.py',
             'research/ic_candidate_tournament_20260915/oracle.py',
             'research/ic_candidate_tournament_20260915/identity.py',
             'research/ic_candidate_tournament_20260915/goal_20260924/generic-query-law/run_controls.py',
             'research/ic_candidate_tournament_20260915/goal_20260924/generic-query-law/PROTOCOL.md')
    write(out / 'environment.json', dict(
        host=platform.uname()._asdict(), python=sys.version, timeout_seconds=60,
        rayon_threads=1, rustc=subprocess.check_output(['rustc', '--version'], text=True).strip(),
        source_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        worker_sha256=digest(worker), inputs_sha256=digest(out / 'inputs.json'),
        vectors_sha256=digest(args.vectors), sources={name: digest(ROOT / name) for name in names}))
    write(out / 'vector-receipt.json', verify_rust_vectors(json.loads(args.vectors.read_text())))
    historical = []
    for line in (HARNESS / 'goal_20260924/generic-public-inputs/final-worker-raw.jsonl').read_text().splitlines():
        row = json.loads(line)
        if row['job']['mode'] == 'ic':
            historical.append(dict(name=row['name'], receipt=verify_query_law(
                row['report'], row['report']['fixture'], row['job'])))
    write(out / 'historical-replay.json', historical)
    env = {**os.environ, 'RAYON_NUM_THREADS': '1'}
    fixtures, outcomes = {}, []
    with (out / 'worker-raw.jsonl').open('x') as raw:
        for item in jobs:
            name, job = item['name'], copy.deepcopy(item['job'])
            run = out / name
            run.mkdir()
            receipt = dict(name=name, status='FAIL', expected=item['expected'], promotion_eligible=False)
            report = None
            try:
                key = (job['degree'], job['curve_a'])
                if key not in fixtures:
                    fixture_job = dict(job, mode='fixture')
                    write(out / f'fixture-n{key[0]}a{key[1]}-job.json', fixture_job)
                    fixture_process = subprocess.run([str(worker)], input=json.dumps(fixture_job),
                                                     text=True, capture_output=True, timeout=60, env=env)
                    (out / f'fixture-n{key[0]}a{key[1]}.json').write_text(fixture_process.stdout)
                    require(fixture_process.returncode == 0, 'fixture generation failed')
                    fixtures[key] = json.loads(fixture_process.stdout)['fixture']
                fixture = fixtures[key]
                job['public_targets'] = fixture['targets']
                write(run / 'job.json', job)
                cp = subprocess.run([str(worker)], input=json.dumps(job), text=True,
                                    capture_output=True, timeout=60, env=env)
                (run / 'stdout.json').write_text(cp.stdout)
                (run / 'stderr.log').write_text(cp.stderr)
                receipt['exit_code'] = cp.returncode
                report = json.loads(cp.stdout)
                require(report['status'] in ('complete', 'incomplete'), 'unexpected worker status')
                require(cp.returncode == (0 if report['status'] == 'complete' else 2), 'wrong worker exit')
                if item['expected'] != 'complete-or-incomplete':
                    require(report['status'] == item['expected'], 'unexpected completion state')
                receipt['worker_status'] = report['status']
                receipt['query_law'] = verify_query_law(report, fixture, job)
                receipt['inventory'] = factor_base_inventory(report, fixture)
                receipt['certificate'] = (verify(report, fixture, expected_mode='ic',
                                                 summands=job['config']['summands'])
                                          if report['status'] == 'complete' else None)
                receipt['status'] = 'PASS'
            except subprocess.TimeoutExpired as exc:
                (run / 'timeout-stdout.bin').write_bytes(exc.stdout or b'')
                (run / 'timeout-stderr.bin').write_bytes(exc.stderr or b'')
                receipt['reason'] = 'worker exceeded frozen 60-second limit'
            except Exception as exc:
                receipt['reason'] = f'{type(exc).__name__}: {exc}'
            receipt['artifacts'] = {p.name: digest(p) for p in sorted(run.iterdir()) if p.is_file()}
            write(run / 'receipt.json', receipt)
            raw.write(json.dumps(dict(name=name, job=job, report=report, receipt=receipt),
                                 sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')
            raw.flush()
            outcomes.append(receipt)
            print(name, receipt['status'], receipt.get('reason', ''), flush=True)
    summary = dict(schema_version=1, tests=len(outcomes),
                   passed=sum(r['status'] == 'PASS' for r in outcomes),
                   historical_reports=len(historical), promotion_eligible=False,
                   online_ms=None, cold_ms=None, instruction_cost=None, speedup=None,
                   scope='toy query-law controls, no performance comparison', results=outcomes)
    write(out / 'summary.json', summary)
    require(summary['passed'] == summary['tests'], 'query-law failures retained')


if __name__ == '__main__':
    main()
