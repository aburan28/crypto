#!/usr/bin/env python3
"""Execute the fixed optimized-producer admission protocol; never select a winner."""
import argparse
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from identity import candidate_manifest, sha256, workload_manifest, write_immutable
from measurement import PHASES, exclusive_ledger, measured_run, report_sha256
from oracle import require, verify
from tournament import digest, execute, parse_profiles, read, write
from producer.evidence import audit_stages, check_build_identity, method_record, scientific_ledger
from producer.timing import native_intervals


def process_ok(process):
    require(process['exit_code'] == 0 and process['process_status'] == 'EXITED',
            'worker failed, timed out or did not complete')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--worker', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    prepared, worker, out = args.prepared.resolve(), args.worker.resolve(), args.out.resolve()
    require(platform.system() == 'Linux' and platform.machine() == 'x86_64', 'requires Linux amd64')
    require(subprocess.check_output(['valgrind', '--version'], text=True).strip() == 'valgrind-3.22.0',
            'unqualified profiler')
    out.mkdir(parents=True, exist_ok=False)
    preparation = read(prepared/'preparation.json')
    manifest = read(prepared/'source-manifest.json')
    require(sha256(manifest) == preparation['source_manifest_sha256'], 'changed prepared manifest')
    for relative, expected in manifest.items():
        require(digest(prepared/'source'/relative) == expected, 'changed prepared source')
    reference = preparation['reference']
    cells = [(13, 0), (17, 1), (23, 1)] if reference == 'both' else [(13, 0), (23, 1), (37, 0), (43, 1), (61, 1)]
    jobs = [dict(mode='ic', degree=n, curve_a=a, target_seeds=[2026092526], algorithm_seed=2026092527,
                 factor_base={'kind': 'subgroup_orbits', 'seed': 43, 'points': 6*n},
                 config=dict(solver='pair_table', linear_algebra='tiny_gauss', batch_trials=1,
                             max_trials=65536, summands=3)) for n, a in cells]
    cpu = sorted(os.sched_getaffinity(0))[-1]
    resources = dict(cpu=cpu, worker_threads=1, memory_bytes=8*1024**3, timeout_seconds=180)
    compiler = subprocess.check_output(['rustc', '--version'], text=True).strip()
    require(compiler.startswith('rustc 1.94.1 '), 'unqualified compiler')
    require(not any(k in os.environ for k in ('RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS'))
            and not any(k.startswith('CARGO_PROFILE_RELEASE_') for k in os.environ),
            'undeclared override of the archived build flags')
    build = dict(compiler=compiler, target='x86_64-unknown-linux-musl',
                 cargo_config_sha256=manifest['.cargo/config.toml'])
    host = platform.uname()._asdict()
    # Fixture construction is a separate process and is outside every measured
    # IC/rho job. Freeze the resulting points in the protocol before admission.
    for job in jobs:
        directory = out/'fixtures'/f"n{job['degree']}a{job['curve_a']}"
        fixture_process = execute([str(worker)], dict(job, mode='fixture'), directory,
            resources['timeout_seconds'], resources['memory_bytes'], cpu)
        write(directory/'process.json', fixture_process, exclusive=True)
        process_ok(fixture_process)
        fixture_report = read(directory/'stdout.json')
        require(fixture_report['status'] == 'fixture', 'fixture preparation failed')
        job['public_targets'] = fixture_report['fixture']['targets']
    write(out/'protocol.json', dict(scope='fixed integration; no performance comparison',
        jobs=jobs, repetitions=3, rho_controls_per_cell=1,
        resources=resources, build=build, host=host,
        worker_sha256=digest(worker), source_manifest_sha256=sha256(manifest),
        preparation=preparation, evaluator_sha256={str(p.relative_to(HERE.parent)): digest(p)
            for p in [HERE/'integration.py', HERE/'evidence.py', HERE/'timing.py', HERE/'audit.py', HERE.parent/'oracle.py',
                      HERE.parent/'identity.py', HERE.parent/'measurement.py', HERE.parent/'tournament.py']},
        protocol_sha256=digest(HERE/'PROTOCOL.md')), exclusive=True)
    (out/'cpuinfo.txt').write_text(Path('/proc/cpuinfo').read_text())
    # Retain the actual executable and the evaluator with the reports. A source
    # digest alone cannot reproduce the bytes that were measured after the CI
    # runner and its build directory disappear.
    shutil.copy2(worker, out/'worker')
    for relative in ('identity.py', 'measurement.py', 'oracle.py', 'tournament.py',
                     'portfolio.py', 'producer/evidence.py', 'producer/integration.py', 'producer/timing.py', 'producer/audit.py',
                     'producer/PROTOCOL.md'):
        target = out/'evaluator'/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(HERE.parent/relative, target)
    outcomes = []
    rho_outcomes = []
    for job in jobs:
        cell = f"n{job['degree']}a{job['curve_a']}"
        cell_dir = out/cell
        cell_dir.mkdir()
        admission_process = execute([str(worker)], dict(job, mode='inventory'), cell_dir/'admission',
            resources['timeout_seconds'], resources['memory_bytes'], cpu)
        write(cell_dir/'admission/process.json', admission_process, exclusive=True)
        try:
            process_ok(admission_process)
            admission = read(cell_dir/'admission/stdout.json')
            require(admission['status'] == 'inventory', 'inventory admission failed')
            check_build_identity(admission, sha256(manifest))
            require(admission.get('field_kernel') in ('portable', 'pclmulqdq'), 'unknown field backend')
            fixture = admission['fixture']
            require(fixture == read(out/'fixtures'/cell/'stdout.json')['fixture'], 'admission changed prepared public point')
            method = method_record(job, fixture, manifest, sha256(manifest), reference,
                                   dict(build, field_kernel=admission['field_kernel']))
            candidate = candidate_manifest(fixture, admission, method)
            workload = workload_manifest(fixture, input_law='ic-workflow-public-target-v1; public point supplied after separate hash-to-curve-cofactor preparation',
                algorithm_seed=job['algorithm_seed'], resource_envelope=resources)
            write_immutable(cell_dir/'candidate.json', candidate)
            write_immutable(cell_dir/'workload.json', workload)
            write_immutable(cell_dir/'method.json', method)
        except Exception as exc:
            outcomes.append(dict(cell=cell, status='ADMISSION_FAILED', reason=f'{type(exc).__name__}: {exc}'))
            write(cell_dir/'admission/failure.json', outcomes[-1], exclusive=True)
            continue
        for repetition in range(3):
            directory = cell_dir/f'rep-{repetition}'
            directory.mkdir()
            write(directory/'job.json', job, exclusive=True)
            result = dict(cell=cell, repetition=repetition, status='FAILED', promotion_eligible=False)
            native = profile = None
            try:
                native = execute([str(worker)], job, directory/'native', 180, 8*1024**3, cpu)
                write(directory/'native/process.json', native, exclusive=True)
                process_ok(native)
                native_report = read(directory/'native/stdout.json')
                check_build_identity(native_report, sha256(manifest), admission['field_kernel'])
                proof = verify(native_report, fixture, summands=3)
                native_audit = audit_stages(native_report, fixture, job['algorithm_seed'])
                command = ['valgrind', '--tool=callgrind', '--cache-sim=no', '--branch-sim=no',
                           '--separate-threads=no', '--collect-atstart=yes', '--instr-atstart=yes',
                           '--callgrind-out-file='+str(directory/'profile/callgrind.out'), str(worker)]
                profile = execute(command, job, directory/'profile', 180, 8*1024**3, cpu)
                write(directory/'profile/process.json', profile, exclusive=True)
                process_ok(profile)
                report = read(directory/'profile/stdout.json')
                check_build_identity(report, sha256(manifest), admission['field_kernel'])
                require(verify(report, fixture, summands=3) == proof, 'native/profile proof differs')
                audit = audit_stages(report, fixture, job['algorithm_seed'])
                require(audit == native_audit, 'native/profile matrix or query diagnostics differ')
                costs = parse_profiles(directory/'profile', phase_schema=3)
                ledger = scientific_ledger(report, costs)
                record = measured_run(candidate=candidate, workload=workload, report=report, fixture=fixture,
                    method=method, number=repetition, ledger=ledger, native_wall_ns=native['process_wall_ns'],
                    status='complete', diagnostics=audit['queries'], admission_report=admission,
                    provenance=dict(source_manifest_sha256=sha256(manifest), worker_sha256=digest(worker),
                        host_id=sha256(host), resource_envelope_id=sha256(resources),
                        calibration_id='valgrind-3.22-amd64-Ir', report_sha256=report_sha256(report)))
                record['stage_audit'] = audit
                record['peak_rss_bytes'] = native['peak_rss_bytes']
                record['native_timing'] = native_intervals(native_report, native['process_wall_ns'])
                write_immutable(directory/'run.json', record)
                result.update(status='VERIFIED', candidate_id=candidate['candidate_id'], run_id=record['run_id'])
            except Exception as exc:
                result['reason'] = f'{type(exc).__name__}: {exc}'
                processes = [p for p in (native, profile) if p is not None]
                status = 'timeout' if any(p['process_status'] == 'TIMEOUT' for p in processes) else 'error'
                errors = '\n'.join(p.read_text(errors='replace') for p in directory.rglob('stderr.txt'))
                if status != 'timeout' and 'memory allocation' in errors:
                    status = 'oom'
                # A rejected/raw report remains an artifact; it cannot certify
                # a failed run or change its already frozen candidate identity.
                record = measured_run(candidate=candidate, workload=workload, report=None, fixture=fixture,
                    method=method, number=repetition, admission_report=admission,
                    ledger=exclusive_ledger(dict.fromkeys(PHASES), unit='valgrind-3.22-amd64-Ir',
                                             process_operations=None, zero_reasons={}),
                    native_wall_ns=native['process_wall_ns'] if native else None, status=status,
                    provenance=dict(source_manifest_sha256=sha256(manifest), worker_sha256=digest(worker),
                        host_id=sha256(host), resource_envelope_id=sha256(resources),
                        calibration_id='valgrind-3.22-amd64-Ir', report_sha256=None))
                record['failure_reason'] = result['reason']
                write_immutable(directory/'run.json', record)
                result.update(candidate_id=candidate['candidate_id'], run_id=record['run_id'])
            result['artifacts'] = {str(p.relative_to(directory)): digest(p)
                                   for p in sorted(directory.rglob('*')) if p.is_file()}
            write(directory/'receipt.json', result, exclusive=True)
            outcomes.append({k: v for k, v in result.items() if k != 'artifacts'})
            print(json.dumps(outcomes[-1]), flush=True)
        # This validates the matched point and reference timing boundary. It
        # selects no rho width and does not compare performance across sources.
        directory = cell_dir/'rho-control'
        directory.mkdir()
        rho_job = dict(job, mode='rho')
        write(directory/'job.json', rho_job, exclusive=True)
        rho_result = dict(cell=cell, status='FAILED', promotion_eligible=False)
        try:
            native = execute([str(worker)], rho_job, directory/'native', 180, 8*1024**3, cpu)
            write(directory/'native/process.json', native, exclusive=True)
            process_ok(native)
            native_report = read(directory/'native/stdout.json')
            check_build_identity(native_report, sha256(manifest), admission['field_kernel'])
            proof = verify(native_report, fixture, expected_mode='rho')
            timing = native_intervals(native_report, native['process_wall_ns'])
            command = ['valgrind', '--tool=callgrind', '--cache-sim=no', '--branch-sim=no',
                       '--separate-threads=no', '--collect-atstart=yes', '--instr-atstart=yes',
                       '--callgrind-out-file='+str(directory/'profile/callgrind.out'), str(worker)]
            profile = execute(command, rho_job, directory/'profile', 180, 8*1024**3, cpu)
            write(directory/'profile/process.json', profile, exclusive=True)
            process_ok(profile)
            report = read(directory/'profile/stdout.json')
            check_build_identity(report, sha256(manifest), native_report['field_kernel'])
            require(report.get('rho_reusable_setup_excluded') is True, 'profile rho preparation boundary differs')
            require(verify(report, fixture, expected_mode='rho') == proof, 'rho native/profile proof differs')
            costs = parse_profiles(directory/'profile', phase_schema=3)
            require(set(costs) == {'setup', 'reference_solve', 'recovery_check'}, 'rho phase closure failed')
            rho_result.update(status='VERIFIED', native_timing=timing, instruction_phases=costs,
                              certificate=proof, workload_id=workload['workload_id'])
        except Exception as exc:
            rho_result['reason'] = f'{type(exc).__name__}: {exc}'
        rho_result['artifacts'] = {str(p.relative_to(directory)): digest(p)
                                  for p in sorted(directory.rglob('*')) if p.is_file()}
        write(directory/'receipt.json', rho_result, exclusive=True)
        rho_outcomes.append({k: v for k, v in rho_result.items() if k != 'artifacts'})
        print(json.dumps(dict(cell=cell, rho_status=rho_result['status'])), flush=True)
    summary = dict(scope='optimized producer integration; no performance comparison', promotion_eligible=False,
                   scheduled=len(jobs)*3, verified=sum(r['status'] == 'VERIFIED' for r in outcomes), outcomes=outcomes,
                   rho_scheduled=len(jobs), rho_verified=sum(r['status'] == 'VERIFIED' for r in rho_outcomes),
                   rho_outcomes=rho_outcomes)
    write(out/'summary.json', summary, exclusive=True)
    require(summary['verified'] == summary['scheduled'], 'retained producer admission failures')
    require(summary['rho_verified'] == summary['rho_scheduled'], 'retained rho timing-control failures')


if __name__ == '__main__':
    main()
