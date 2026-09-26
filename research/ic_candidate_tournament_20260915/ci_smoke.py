#!/usr/bin/env python3
"""Bounded integration checks; deliberately no performance ranking or promotion."""
import argparse
import copy
import json
import platform
import subprocess
from pathlib import Path

from generic_queries import verify_queries
from generic_query_law import verify_query_law
from generic_phases import verify_native as verify_native_phases, parse_profiles as parse_generic_profiles
from generic_admission import admit, admit_rho, scientific_ledger
from generic_build import verify_build_record
from identity import curve_record, factor_base_inventory, write_immutable
from measurement import legacy_ledger
from oracle import require, verify
from tournament import PHASES, child_env, digest, execute, parse_profiles, read, write

SOLVERS = ('pair_table', 'enumerate', 'f4', 'f5', 'inherited_f4', 'sat_xor', 'sat_cnf')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', type=Path, required=True)
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    worker = args.worker.resolve()
    build = read(args.build_dir / 'build-record.json')
    source = read(args.build_dir / 'source-manifest.json')
    verify_build_record(build, source)
    out = args.out.resolve()
    require(platform.system() == 'Linux' and platform.machine() == 'x86_64',
            'CI integration requires Linux amd64')
    version = subprocess.check_output(['valgrind', '--version'], text=True).strip()
    require(version == 'valgrind-3.22.0', 'unexpected Valgrind version')
    out.mkdir(parents=True, exist_ok=False)
    import os
    cpu = sorted(os.sched_getaffinity(0))[-1]
    config = dict(solver='pair_table', linear_algebra='sparse', batch_trials=8,
                  max_trials=256, summands=2, factor_base=dict(kind='factor', index=0))
    jobs = []
    for a in (0, 1):
        template = dict(mode='ic', degree=9, curve_a=a, target_seeds=[2026092425],
                        algorithm_seed=2026092425, factor_base=dict(kind='subgroup_orbits',
                                                                  seed=43, points=54),
                        config=config)
        for solver in SOLVERS:
            for la in ('dense', 'sparse'):
                job = copy.deepcopy(template)
                job['config'].update(solver=solver, linear_algebra=la)
                jobs.append((f'algebra-a{a}-{solver}-{la}', job, False))
        job = copy.deepcopy(template)
        job['mode'] = 'rho'
        jobs.append((f'rho-a{a}', job, True))
        job = copy.deepcopy(template)
        jobs.append((f'profile-ic-a{a}', job, True))
    for name, recipe in (
        ('subgroup', dict(kind='subgroup_orbits', seed=43, points=78)),
        ('frobenius-union', dict(kind='frobenius_union', seed_masks=[1, 2, 4, 8])),
    ):
        job = dict(mode='ic', degree=13, curve_a=0, target_seeds=[2026092425],
                   algorithm_seed=2026092425, factor_base=recipe,
                   config=dict(solver='pair_table', linear_algebra='sparse',
                               batch_trials=1, max_trials=4096, summands=3))
        jobs.append((f'factor-base-{name}', job, False))
    for window in (0, 1, 1000):
        job = dict(mode='ic', degree=9, curve_a=0, target_seeds=[2026092425],
                   algorithm_seed=2026092425,
                   factor_base=dict(kind='subgroup_orbits', seed=43, points=36),
                   config=dict(solver='pair_table', linear_algebra='dense', summands=3,
                               batch_trials=128, max_trials=256, collection_window=window))
        jobs.append((f'query-law-window{window}', job, window == 1))
    for solver in SOLVERS:
        job = copy.deepcopy(jobs[0][1])
        job['config'].update(solver=solver, batch_trials=1, max_trials=1)
        jobs.append((f'incomplete-{solver}', job, False))
    # This is a frozen list of test vectors, not a search over candidate winners.
    write(out / 'inputs.json', [{'test': name, 'job': job, 'profile': profile}
                               for name, job, profile in jobs], exclusive=True)
    write(out / 'environment.json', {
        'scope': 'integration tests; no candidate comparison or promotion',
        'worker_sha256': digest(worker), 'profiler': version,
        'host': platform.uname()._asdict(), 'timeout_seconds': 60,
        'memory_bytes': 8 * 1024**3, 'cpu': cpu,
        'input_sha256': digest(out / 'inputs.json'),
        'evaluator_sha256': {name: digest(Path(__file__).with_name(name))
                             for name in ('ci_smoke.py', 'tournament.py', 'oracle.py', 'portfolio.py',
                                          'identity.py', 'measurement.py', 'generic_queries.py',
                                          'generic_query_law.py', 'generic_phases.py')},
    }, exclusive=True)
    outcomes = []
    fixtures = {}
    for name, job, profile in jobs:
        directory = out / name
        directory.mkdir()
        write(directory / 'job.json', job, exclusive=True)
        expected_failure = name.startswith('incomplete-')
        receipt = dict(test=name, status='FAILED', expected_failure=expected_failure, promotion_eligible=False,
                       end_to_end_speedup=None)
        try:
            key = (job['degree'], job['curve_a'])
            if key not in fixtures:
                fixture_job = dict(job, mode='fixture')
                write(directory / 'fixture-job.json', fixture_job, exclusive=True)
                fixture_process = execute([str(worker)], fixture_job, directory / 'fixture',
                                          60, 8 * 1024**3, cpu)
                write(directory / 'fixture-process.json', fixture_process, exclusive=True)
                require(fixture_process['exit_code'] == 0 and
                        fixture_process['process_status'] == 'EXITED', 'fixture generation failed')
                fixtures[key] = read(directory / 'fixture/stdout.json')['fixture']
            fixture = fixtures[key]
            write(directory / 'expected-fixture.json', fixture, exclusive=True)
            job['public_targets'] = fixture['targets']
            write(directory / 'measured-job.json', job, exclusive=True)
            write_immutable(directory / 'curve-manifest.json', curve_record(fixture))
            process = execute([str(worker)], job, directory / 'native', 60, 8 * 1024**3, cpu)
            receipt['native_process'] = process
            require(process['exit_code'] == (2 if expected_failure else 0)
                    and process['process_status'] == 'EXITED', 'unexpected native exit')
            report = read(directory / 'native/stdout.json')
            require(report['online_timing_schema'] == 1 and
                    report['target_input'] == 'supplied_public_point' and
                    report['reusable_setup_excluded'] is True, 'missing public online boundary')
            if not expected_failure:
                require(type(report['online_wall_ns']) is int and report['online_wall_ns'] > 0
                        and report['scalar_replay_included'] is True, 'online interval missing replay')
            else:
                require(report['online_wall_ns'] is None, 'preparation failure gained online time')
            if expected_failure:
                require(report['status'] == 'incomplete' and report['trials'] == 1,
                        'intentional exhaustion changed')
                require(report['solve_attempts'] == report['log_table_report']['solve_attempts'],
                        'failed matrix attempts lost')
                proof = None
            else:
                proof = verify(report, fixture,
                               expected_mode=job['mode'], summands=job['config']['summands'])
            receipt['certificate'] = proof
            if job['mode'] == 'ic':
                receipt['query_accounting'] = verify_queries(report, fixture, job['config']['summands'])
                receipt['query_law'] = verify_query_law(report, fixture, job)
                inventory = factor_base_inventory(report, fixture)
                write_immutable(directory / 'factor-base-inventory.json', inventory)
                receipt['factor_base_inventory'] = inventory
            if profile:
                profile_dir = directory / 'profile'
                command = ['valgrind', '--tool=callgrind', '--cache-sim=no', '--branch-sim=no',
                           '--separate-threads=no', '--collect-atstart=yes', '--instr-atstart=yes',
                           '--callgrind-out-file=' + str(profile_dir / 'callgrind.out'), str(worker)]
                process = execute(command, job, profile_dir, 60, 8 * 1024**3, cpu)
                receipt['profile_process'] = process
                require(process['exit_code'] == 0 and process['process_status'] == 'EXITED',
                        'profiled worker failed or incomplete')
                profile_report = read(profile_dir / 'stdout.json')
                if job['mode'] == 'rho':
                    require(profile_report['field_kernel'] == report['field_kernel'],
                            'native/profile rho arithmetic dispatch differs')
                profile_proof = verify(profile_report, fixture,
                                       expected_mode=job['mode'], summands=job['config']['summands'])
                if job['mode'] == 'ic':
                    profile_queries = verify_queries(read(profile_dir / 'stdout.json'), fixture,
                                                     job['config']['summands'])
                    require(profile_queries == receipt['query_accounting'],
                            'native/profile query accounting differs')
                    require(verify_query_law(profile_report, fixture, job) == receipt['query_law'],
                            'native/profile query law differs')
                require(profile_proof == proof, 'native/profile certificates differ')
                intervals = parse_profiles(profile_dir)
                expected = {'startup_and_input', 'curve_and_targets', 'final_verification',
                            'reporting_and_cleanup'}
                expected |= {'rho_solve'} if job['mode'] == 'rho' else PHASES - {'rho_solve'}
                require(set(intervals) == expected, 'missing exclusive profiler interval')
                require(sum(intervals.values()) > 0, 'empty instruction trace')
                receipt['profiler_interval_closure_verified'] = True
                # Existing intervals combine some scientific stages. They must not
                # be relabelled as the new exclusive T_* phase decomposition.
                receipt['raw_profiler_intervals'] = intervals
                if job['mode'] == 'ic':
                    ledger = legacy_ledger(intervals, unit='valgrind-3.22-amd64-Ir',
                                           process_operations=sum(intervals.values()))
                    require(ledger['cold_operations'] is None, 'legacy stages became a scientific total')
                    write(directory / 'scientific-phase-ledger.json', ledger, exclusive=True)
            # The opt-in generic tracer has its own namespace and schema. Its
            # intervals must close independently and preserve the untraced job.
            traced_job = dict(job, exclusive_phases=True)
            write(directory / 'exclusive-job.json', traced_job, exclusive=True)
            traced_dir = directory / 'exclusive-native'
            traced_process = execute([str(worker)], traced_job, traced_dir, 60, 8 * 1024**3, cpu)
            receipt['exclusive_native_process'] = traced_process
            require(traced_process['exit_code'] == (2 if expected_failure else 0)
                    and traced_process['process_status'] == 'EXITED', 'exclusive native exit changed')
            traced_report = read(traced_dir / 'stdout.json')
            native_ns = round(traced_process['process_wall_seconds'] * 1e9)
            if job['mode'] == 'ic':
                admitted = admit(traced_report, fixture, traced_job, build, source,
                    executable=worker, process_wall_ns=native_ns,
                    resources=dict(cpu=cpu, rayon_threads=1, memory_bytes=8*1024**3, timeout_seconds=60))
                receipt['scientific_admission'] = admitted
                for artifact in ('candidate', 'workload', 'run'):
                    write_immutable(directory / f'generic-{artifact}.json', admitted[artifact])
            else:
                receipt['scientific_admission'] = admit_rho(traced_report, fixture, traced_job,
                    build, source, executable=worker, process_wall_ns=native_ns)
            receipt['exclusive_native_phases'] = verify_native_phases(
                traced_report, traced_job,
                process_wall_ns=round(traced_process['process_wall_seconds'] * 1e9))
            traced_proof = None if expected_failure else verify(
                traced_report, fixture, expected_mode=job['mode'], summands=job['config']['summands'])
            require(traced_proof == proof, 'tracing changed certificates')
            if job['mode'] == 'ic':
                require(verify_queries(traced_report, fixture, job['config']['summands'])
                        == receipt['query_accounting'], 'tracing changed queries/counters')
                receipt['exclusive_query_law'] = verify_query_law(traced_report, fixture, traced_job)
            # Cover each backend, sparse LA, both rho curves, a collection walk,
            # and an unstarted target under Callgrind without profiling every
            # duplicate LA/curve combination in this integration panel.
            profile_exclusive = (profile or name == 'incomplete-pair_table' or
                                 (name.startswith('algebra-a0-') and job['config']['linear_algebra'] == 'dense'))
            if profile_exclusive:
                traced_profile = directory / 'exclusive-profile'
                command = ['valgrind', '--tool=callgrind', '--cache-sim=no', '--branch-sim=no',
                           '--separate-threads=no', '--collect-atstart=yes', '--instr-atstart=yes',
                           '--callgrind-out-file=' + str(traced_profile / 'callgrind.out'), str(worker)]
                profile_process = execute(command, traced_job, traced_profile, 60, 8 * 1024**3, cpu)
                receipt['exclusive_profile_process'] = profile_process
                require(profile_process['exit_code'] == (2 if expected_failure else 0)
                        and profile_process['process_status'] == 'EXITED', 'exclusive profiler exit changed')
                profiled = read(traced_profile / 'stdout.json')
                receipt['exclusive_instructions'] = parse_generic_profiles(traced_profile, profiled, traced_job)
                if job['mode'] == 'ic':
                    profiled_admission = admit(profiled, fixture, traced_job, build, source,
                        executable=worker, process_wall_ns=round(profile_process['process_wall_seconds'] * 1e9),
                        resources=dict(cpu=cpu, rayon_threads=1, memory_bytes=8*1024**3, timeout_seconds=60))
                    require(profiled_admission['candidate'] == receipt['scientific_admission']['candidate'],
                            'native/profile scientific method mismatch')
                    instructions = receipt['exclusive_instructions']
                    receipt['scientific_instruction_ledger'] = scientific_ledger(
                        instructions['phases'], unit=instructions['unit'],
                        process_total=instructions['process_instructions'])
                else:
                    admit_rho(profiled, fixture, traced_job, build, source, executable=worker,
                              process_wall_ns=round(profile_process['process_wall_seconds'] * 1e9))
                profile_proof = None if expected_failure else verify(
                    profiled, fixture, expected_mode=job['mode'], summands=job['config']['summands'])
                require(profile_proof == traced_proof, 'exclusive native/profile certificate mismatch')
                if job['mode'] == 'ic':
                    require(verify_query_law(profiled, fixture, traced_job) == receipt['exclusive_query_law'],
                            'exclusive native/profile query law mismatch')
                else:
                    require(profiled['field_kernel'] == traced_report['field_kernel'],
                            'exclusive native/profile rho dispatch mismatch')
            receipt['status'] = 'VERIFIED'
        except Exception as exc:
            receipt['reason'] = f'{type(exc).__name__}: {exc}'
        receipt['artifacts'] = {str(path.relative_to(directory)): digest(path)
                                for path in sorted(directory.rglob('*')) if path.is_file()}
        write(directory / 'receipt.json', receipt, exclusive=True)
        outcomes.append(dict(test=name, status=receipt['status'], reason=receipt.get('reason')))
        print(json.dumps(outcomes[-1]), flush=True)
    summary = dict(scope='integration, no performance claim', promotion_eligible=False,
                   tests=len(outcomes), verified=sum(row['status'] == 'VERIFIED' for row in outcomes),
                   outcomes=outcomes)
    write(out / 'summary.json', summary, exclusive=True)
    require(summary['verified'] == summary['tests'], 'integration failures retained in evidence')


if __name__ == '__main__':
    main()
