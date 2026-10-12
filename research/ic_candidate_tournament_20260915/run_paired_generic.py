#!/usr/bin/env python3
"""Run one registered source-bound F5 or rho arm; retain every outcome."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time

from generic_admission import admit, admit_rho
from generic_build import verify_build_record
from generic_solver_feasibility import check_source_checkout
from identity import sha256, workload_manifest
from oracle import require
from register_paired_generic import (HERE, INPUT_LAW, REGISTRATION,
                                     RESOURCES, digest, identities)
from tournament import read, write


def sealed_inputs(arm, source_root):
    require(arm in ('f5', 'rho'), 'unknown frozen arm')
    data = identities()
    seal = read(REGISTRATION/'seal.json')
    require(data['panel'] == read(REGISTRATION/'panel.json')
            and data['method'] == read(REGISTRATION/'method.json')
            and data['candidate'] == read(REGISTRATION/'candidate.json')
            and data['f5_workload'] == read(REGISTRATION/'f5-workload.json')
            and data['rho_workload'] == read(REGISTRATION/'rho-workload.json')
            and data['f5_job'] == read(REGISTRATION/'f5-job.json')
            and data['rho_job'] == read(REGISTRATION/'rho-job.json')
            and seal['files'] == {name: digest(REGISTRATION/name)
                for name in seal['files']}
            and seal['controller_sha256'] == digest(Path(__file__))
            and seal['registration_sha256']
                == digest(HERE/'register_paired_generic.py'),
            'generic paired registration or controller changed')
    check_source_checkout(source_root)
    require(sha256(data['build']) == data['panel']['build_record_sha256']
            and data['build']['source_manifest_sha256']
                == data['panel']['source_manifest_sha256']
            and hashlib.sha256(data['files']['build/worker']).hexdigest()
                == data['panel']['worker_sha256'],
            'registered worker source/build identity changed')
    verify_build_record(data['build'], data['source'])
    return data


def run_worker(worker, job, out):
    environment = {name: os.environ[name] for name in
                   ('PATH', 'HOME', 'TMPDIR', 'LANG', 'LC_ALL')
                   if name in os.environ}
    environment['RAYON_NUM_THREADS'] = '1'
    started = time.monotonic_ns()
    with (out/'stdout.json').open('x') as stdout, \
         (out/'stderr.txt').open('x') as stderr:
        process = subprocess.Popen([str(worker)], stdin=subprocess.PIPE,
                                   stdout=stdout, stderr=stderr, text=True,
                                   env=environment)
        process.stdin.write(json.dumps(job, sort_keys=True))
        process.stdin.close()
        if hasattr(os, 'wait4'):
            _, status, usage = os.wait4(process.pid, 0)
            exit_code = os.waitstatus_to_exitcode(status)
            process.returncode = exit_code
            peak_rss_bytes = (usage.ru_maxrss if platform.system() == 'Darwin'
                              else usage.ru_maxrss * 1024
                              if platform.system() == 'Linux' else None)
            peak_scope = 'OS child-process high-water RSS via wait4'
        else:
            exit_code = process.wait()
            peak_rss_bytes = None
            peak_scope = 'wait4 unavailable; peak RSS unknown'
    wall_ns = time.monotonic_ns()-started
    result = dict(exit_code=exit_code, process_wall_ns=wall_ns,
                  process_start_monotonic_ns=started,
                  memory_peak_bytes=peak_rss_bytes,
                  memory_scope=peak_scope+'; resource limit not enforced',
                  environment_policy='PATH/HOME/TMPDIR/LANG/LC_ALL and RAYON_NUM_THREADS=1',
                  hard_process_wall_limit_seconds=None,
                  hard_memory_limit_bytes=None)
    write(out/'process.json', result, exclusive=True)
    return result


def audit_f5(report, job, process, data, worker):
    require(report['status'] in ('complete', 'incomplete'),
            'F5 worker returned no IC report')
    require((process['exit_code'] == 0 and report['status'] == 'complete')
            or (process['exit_code'] == 2 and report['status'] == 'incomplete'),
            'F5 process exit contradicts report status')
    expected = dict(data['old']['fixture'],
                    targets=job['public_targets'],
                    target_seeds=job['target_seeds'],
                    target_scalar_constructed=False)
    require(report['fixture'] == expected,
            'F5 worker substituted curve, target or seed')
    admitted = admit(report, expected, job, data['build'], data['source'],
                     executable=worker,
                     process_wall_ns=process['process_wall_ns'],
                     resources=RESOURCES, number=0)
    require(admitted['candidate'] == data['candidate']
            and admitted['workload'] == data['f5_workload']
            and admitted['run']['run_id'] == data['panel']['f5_run_id']
            and (admitted['run']['certificate'] is not None)
                == (report['status'] == 'complete'),
            'F5 candidate/workload identity or scalar certificate differs')
    return admitted


def audit_rho(report, job, process, data, worker):
    require(process['exit_code'] == 0 and report['status'] == 'complete',
            'rho process did not complete one target')
    expected = dict(data['old']['fixture'],
                    targets=job['public_targets'],
                    target_seeds=job['target_seeds'],
                    target_scalar_constructed=False)
    require(report['fixture'] == expected,
            'rho worker substituted curve, target or seed')
    admitted = admit_rho(report, expected, job, data['build'], data['source'],
                         executable=worker,
                         process_wall_ns=process['process_wall_ns'])
    workload = workload_manifest(
        expected, input_law=INPUT_LAW,
        algorithm_seed=job['algorithm_seed'],
        resource_envelope=RESOURCES, cache_policy='cold')
    require(workload == data['rho_workload']
            and admitted['certificate']['verified_targets'] == 1
            and admitted['phases']['online_wall_ns'] is not None,
            'rho public workload or verified online interval differs')
    return admitted


def run(arm, source_root, out):
    data = sealed_inputs(arm, source_root)
    require(not out.exists(), 'one-shot paired-arm output already exists')
    out.mkdir(parents=True)
    for name in ('panel.json', 'seal.json', 'method.json',
                 'candidate.json', 'f5-workload.json', 'rho-workload.json',
                 'f5-job.json', 'rho-job.json', 'PROTOCOL.md'):
        shutil.copy2(REGISTRATION/name, out/name)
    shutil.copy2(Path(__file__), out/'registered-runner.py')
    shutil.copy2(HERE/'register_paired_generic.py',
                 out/'registered-registration.py')
    write(out/'source-manifest.json', data['source'], exclusive=True)
    write(out/'build-record.json', data['build'], exclusive=True)
    (out/'root-source.tar.gz').write_bytes(
        data['files']['build/root-source.tar.gz'])
    worker = out/'worker'
    worker.write_bytes(data['files']['build/worker'])
    worker.chmod(0o755)
    require(digest(worker) == data['panel']['worker_sha256'],
            'copied worker binary changed')
    identity = json.loads(subprocess.check_output(
        [str(worker), '--build-identity'], text=True))
    require(identity == data['build']['identity'],
            'worker binary embeds a different build identity')
    write(out/'host.json', dict(system=platform.system(),
         machine=platform.machine(), processor=platform.processor(),
         cpu_count=os.cpu_count(), scope='physical local one-thread diagnostic'),
         exclusive=True)
    job = data[arm+'_job']
    write(out/'executed-job.json', job, exclusive=True)
    process = run_worker(worker, job, out)
    result = dict(schema_version=1, arm=arm,
                  status='PRODUCER_FAILURE',
                  process=process,
                  candidate_id=(data['panel']['candidate_id'] if arm == 'f5'
                                else None),
                  reference_id=(data['panel']['rho_reference_id'] if arm == 'rho'
                                else None),
                  workload_id=data[arm+'_workload']['workload_id'],
                  run_id=data['panel'][arm+'_run_id'],
                  target_point=data['target_panel']['target_input']['point'],
                  verified=False, online_wall_ns=None,
                  paired_rho_speedup=None, promotion_eligible=False)
    try:
        report = json.loads((out/'stdout.json').read_text())
        admitted = (audit_f5(report, job, process, data, worker)
                    if arm == 'f5' else
                    audit_rho(report, job, process, data, worker))
        write(out/'admission.json', admitted, exclusive=True)
        complete = report['status'] == 'complete'
        result.update(status=('AUDITED_COMPLETE' if complete
                              else 'AUDITED_BOUNDED_INCOMPLETE'),
                      verified=complete,
                      online_wall_ns=(admitted['run']['online_wall_ns']
                                      if arm == 'f5' and complete else
                                      admitted['phases']['online_wall_ns']
                                      if arm == 'rho' and complete else None))
    except Exception as error:
        result.update(audit_error=f'{type(error).__name__}: {error}')
    write(out/'result.json', result, exclusive=True)
    print(json.dumps(result, sort_keys=True), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arm', choices=('f5', 'rho'), required=True)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.arm, args.source_root.resolve(), args.out.resolve())


if __name__ == '__main__':
    main()
