"""Frozen public toy admission controls, never a performance tournament."""
import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[4]
HARNESS = ROOT / 'research/ic_candidate_tournament_20260915'
sys.path.insert(0, str(HARNESS))
from generic_admission import admit, admit_rho
from execution_ids import allocation, audit_runs
from generic_build import digest, verify_binding, verify_build_record
from generic_stages import verify_stages
from identity import write_immutable
from oracle import require


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')


def panel():
    path = HARNESS / 'goal_20260924/generic-query-law/run_controls.py'
    spec = importlib.util.spec_from_file_location('query_panel', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    jobs = module.panel()
    for row in jobs:
        row['job'].update(target_seeds=[2026092557], algorithm_seed=2026092557, exclusive_phases=True)
    for a in (0, 1):
        job = copy.deepcopy(jobs[0]['job'])
        job.update(mode='rho', curve_a=a)
        jobs.append(dict(name=f'rho-a{a}', job=job, expected='complete'))
    # Independent rank replay of the first panel showed rank 1 with three
    # novel rows at query 4. These controls exercise a paid, failed LA call.
    for la in ('dense', 'sparse'):
        for limit in (4, 256):
            job = copy.deepcopy(jobs[0]['job'])
            job['config'].update(linear_algebra=la, batch_trials=1, max_trials=limit)
            jobs.append(dict(name=f'failed-la-{la}-limit{limit}', job=job,
                             expected='incomplete' if limit == 4 else 'complete'))
    for a in (0, 1):
        job = copy.deepcopy(jobs[0]['job'])
        job['curve_a'] = a
        job['config'].update(linear_algebra='sparse', sparse=dict(filter=dict(
            remove_singletons=False, merge_max_weight=0)))
        jobs.append(dict(name=f'sparse-core-a{a}', job=job, expected='complete'))
    for degree, retained in ((9, 9), (13, 3)):
        union = dict(kind='frobenius_union', seed_masks=[1, 2, 4, 8])
        saturated = dict(kind='two_torsion_saturated', parent=union)
        pruned = dict(kind='pruned', parent=saturated, retained_abscissa_orbits=[retained])
        recipes = [dict(kind='factor', index=0),
                   dict(kind='divisor', indices=[0, 2] if degree == 9 else [1]),
                   union, dict(kind='subgroup_orbits', seed=43, points=2*degree*3),
                   saturated, pruned, dict(kind='two_torsion_saturated', parent=pruned)]
        for i, recipe in enumerate(recipes):
            job = copy.deepcopy(jobs[0]['job'])
            job.update(mode='inventory', degree=degree, factor_base=recipe)
            jobs.append(dict(name=f'inventory-n{degree}-{i}-{recipe["kind"]}', job=job, expected='inventory'))
    return jobs


def reference_readiness_panel():
    """Admission-only vectors on the five registered development cells."""
    jobs = []
    for n, a in ((17, 1), (19, 0), (23, 0), (23, 1), (31, 0)):
        for mode, la in (('ic', 'dense'), ('ic', 'sparse'), ('rho', 'dense')):
            job = dict(mode=mode, degree=n, curve_a=a, target_seeds=[2026092561],
                algorithm_seed=2026092561, exclusive_phases=True,
                factor_base=dict(kind='subgroup_orbits', seed=43, points=16*n),
                config=dict(solver='pair_table', linear_algebra=la, summands=3,
                    batch_trials=1, max_trials=65536, collection_window=0, rho_parallel_walks=4))
            jobs.append(dict(name=f'readiness-n{n}a{a}-{mode}-{la}', job=job,
                             expected='complete-or-incomplete'))
    return jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--run-number-start', type=int, required=True)
    parser.add_argument('--failed-la-only', action='store_true')
    parser.add_argument('--sparse-core-only', action='store_true')
    parser.add_argument('--reference-readiness', action='store_true')
    args = parser.parse_args()
    build_dir, out = args.build.resolve(), args.out.resolve()
    worker = build_dir / 'worker'
    build = json.loads((build_dir / 'build-record.json').read_text())
    source = json.loads((build_dir / 'source-manifest.json').read_text())
    verify_build_record(build, source)
    jobs = panel()
    if args.failed_la_only:
        jobs = [item for item in jobs if item['name'].startswith('failed-la-')]
    if args.sparse_core_only:
        require(not args.failed_la_only, 'select only one supplemental panel')
        jobs = [item for item in jobs if item['name'].startswith('sparse-core-')]
    if args.reference_readiness:
        require(not args.failed_la_only and not args.sparse_core_only, 'select only one panel')
        jobs = reference_readiness_panel()
    plan = allocation(args.run_number_start, [item['name'] for item in jobs])
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'run-number-allocation.json', plan)
    write(out / 'inputs.json', jobs)
    resources = dict(timeout_seconds=60, rayon_threads=1, memory_bytes=None, cpu=None)
    preexec = None
    if platform.system() == 'Linux':
        import resource
        cpu = sorted(os.sched_getaffinity(0))[-1]
        resources.update(cpu=cpu, memory_bytes=8*1024**3)
        def restrict():
            os.sched_setaffinity(0, {cpu})
            resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
        preexec = restrict
    env = dict(PATH=os.environ['PATH'], RAYON_NUM_THREADS='1', LC_ALL='C')
    write(out / 'environment.json', dict(host=platform.uname()._asdict(), resources=resources,
        panel='reference-readiness' if args.reference_readiness else 'scientific-controls',
        build_record_sha256=digest(build_dir / 'build-record.json'), worker_sha256=digest(worker),
        sources={name: digest(HARNESS / name) for name in ('generic_admission.py', 'generic_build.py',
        'generic_bases.py', 'generic_stages.py', 'generic_queries.py', 'generic_query_law.py',
        'generic_phases.py', 'identity.py', 'oracle.py', 'measurement.py', 'execution_ids.py',
        'goal_20260924/generic-query-law/run_controls.py',
        'goal_20260924/generic-scientific-admission/run_controls.py',
        'goal_20260924/generic-scientific-admission/PROTOCOL.md',
        'goal_20260924/generic-reference-readiness/PROTOCOL.md')}))
    fixtures, results, admitted_runs = {}, [], []
    with (out / 'worker-raw.jsonl').open('x') as raw:
        for item, execution in zip(jobs, plan['executions'], strict=True):
            name, job = item['name'], copy.deepcopy(item['job'])
            directory = out / name
            directory.mkdir()
            report = None
            receipt = dict(name=name, status='FAIL', expected=item['expected'], promotion_eligible=False)
            process_ns = None
            try:
                key = (job['degree'], job['curve_a'])
                if key not in fixtures:
                    fixture_job = dict(job, mode='fixture', exclusive_phases=False)
                    write(out / f'fixture-n{key[0]}a{key[1]}-job.json', fixture_job)
                    cp = subprocess.run([str(worker)], input=json.dumps(fixture_job), text=True,
                        capture_output=True, timeout=60, env=env, preexec_fn=preexec)
                    (out / f'fixture-n{key[0]}a{key[1]}.json').write_text(cp.stdout)
                    (out / f'fixture-n{key[0]}a{key[1]}.stderr').write_text(cp.stderr)
                    require(cp.returncode == 0, 'fixture generation failed')
                    fixtures[key] = json.loads(cp.stdout)['fixture']
                fixture = fixtures[key]
                job['public_targets'] = fixture['targets']
                write(directory / 'job.json', job)
                started = time.monotonic_ns()
                cp = subprocess.run([str(worker)], input=json.dumps(job), text=True,
                    capture_output=True, timeout=60, env=env, preexec_fn=preexec)
                process_ns = time.monotonic_ns()-started
                (directory / 'stdout.json').write_text(cp.stdout)
                (directory / 'stderr.log').write_text(cp.stderr)
                receipt['exit_code'] = cp.returncode
                report = json.loads(cp.stdout)
                require(report['status'] in {'complete', 'incomplete', 'inventory'}, 'unexpected worker status')
                require(cp.returncode == (2 if report['status'] == 'incomplete' else 0), 'wrong worker exit')
                if item['expected'] != 'complete-or-incomplete':
                    require(report['status'] == item['expected'], 'unexpected completion state')
                if job['mode'] == 'inventory':
                    receipt['binding'] = verify_binding(report, build, source, executable=worker)
                    receipt['stages'] = verify_stages(report, fixture, job)
                elif job['mode'] == 'rho':
                    receipt['admission'] = admit_rho(report, fixture, job, build, source,
                        executable=worker, process_wall_ns=process_ns)
                else:
                    admitted = admit(report, fixture, job, build, source, executable=worker,
                                     process_wall_ns=process_ns, resources=resources,
                                     number=execution['number'])
                    admitted_runs.append(admitted['run'])
                    receipt['admission'] = admitted
                    for artifact in ('candidate', 'workload', 'run'):
                        write_immutable(directory / f'{artifact}.json', admitted[artifact])
                    if name.startswith('sparse-core-'):
                        require(report['sparse_report']['core_dimension'] > 0
                                and report['sparse_report']['wiedemann'] is not None,
                                'sparse core control did not execute block Wiedemann')
                receipt.update(status='PASS', worker_status=report['status'])
            except subprocess.TimeoutExpired as exc:
                (directory / 'timeout-stdout.bin').write_bytes(exc.stdout or b'')
                (directory / 'timeout-stderr.bin').write_bytes(exc.stderr or b'')
                receipt['reason'] = 'worker exceeded frozen 60-second limit'
            except Exception as exc:
                receipt['reason'] = f'{type(exc).__name__}: {exc}'
            write(directory / 'receipt.json', receipt)
            raw.write(json.dumps(dict(name=name, job=job, report=report,
                process_wall_ns=process_ns, receipt=receipt), sort_keys=True,
                separators=(',', ':'), allow_nan=False) + '\n')
            raw.flush()
            results.append(dict(name=name, status=receipt['status'], reason=receipt.get('reason'),
                                worker_status=report.get('status') if report else None))
            print(name, receipt['status'], receipt.get('reason', ''), flush=True)
    # Refuse inherited method changes before a measured session is begun.
    overrides = ('KIC_F4_INHERIT', 'KIC_F4_KERNEL', 'KIC_GF2_SIMD', 'KIC_SCAN_SIMD',
                 'F4_F2_MAX_ROWS', 'SOLVER_SPLIT_RULE', 'IC_TEMPLATE_MEMO',
                 'IC_CACHE_LOCAL_BYTES', 'IC_REDIS_URL', 'IC_F2_BACKEND', 'IC_ARTIFACT_CACHE')
    rejected = []
    first_job = jobs[0]['job']
    control_job = dict(first_job, public_targets=fixtures[(first_job['degree'], first_job['curve_a'])]['targets'])
    for variable in overrides:
        cp = subprocess.run([str(worker)], input=json.dumps(control_job), text=True,
            capture_output=True, timeout=60, env={**env, variable: 'undeclared-control'}, preexec_fn=preexec)
        write(out / f'environment-control-{variable}.json', dict(variable=variable,
            exit_code=cp.returncode, stdout=cp.stdout, stderr=cp.stderr))
        answer = json.loads(cp.stdout)
        require(cp.returncode == 2 and answer['status'] == 'error'
                and 'undeclared algorithm environment override' in answer['reason'],
                'undeclared environment override accepted')
        rejected.append(variable)
    write(out / 'summary.json', dict(schema_version=1, controls=len(results),
        panel='reference-readiness' if args.reference_readiness else 'scientific-controls',
        readiness_passed=(all(row['status'] == 'PASS' and row['worker_status'] == 'complete'
                              for row in results) if args.reference_readiness else None),
        passed=sum(row['status'] == 'PASS' for row in results), performance_qualified=False,
        environment_overrides_rejected=rejected, run_key_audit=audit_runs(admitted_runs),
        promotion_eligible=False, online_speedup=None, results=results))
    require(all(row['status'] == 'PASS' for row in results), 'admission failures retained')
    if args.reference_readiness:
        require(all(row['worker_status'] == 'complete' for row in results),
                'reference-panel incomplete solves retained')


if __name__ == '__main__':
    main()
