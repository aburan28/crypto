#!/usr/bin/env python3
"""Run one presealed local pairinv incumbent or signed rho reference arm."""
import argparse
import json
import os
from pathlib import Path
import platform
import shutil

from campaign_rules import IC_SOURCE
from identity import candidate_manifest, sha256
from oracle import require, verify
from producer.evidence import audit_stages, check_build_identity
from producer.timing import native_intervals
from register_local_pairinv import (REGISTRATION, TARGET, build_inputs,
                                    jobs)
from register_paired_generic import digest
from run_paired_generic import run_worker
from tournament import read, write


def sealed_inputs(prepared, built, arm):
    require(arm in ('ic', 'rho'), 'unknown optimized local arm')
    source, build = build_inputs(prepared, built)
    panel, seal = read(REGISTRATION/'panel.json'), read(REGISTRATION/'seal.json')
    require(panel['status'] == 'REGISTERED_BEFORE_ARM_EXECUTION'
            and panel['accepted_source_manifest_sha256'] == IC_SOURCE
            and panel['build_record_sha256'] == sha256(build)
            and panel['worker_sha256'] == digest(Path(built)/'worker')
            and panel['target_panel_sha256'] == digest(TARGET)
            and seal['controller_sha256'] == digest(__file__)
            and seal['registration_sha256']
                == digest(Path(__file__).with_name('register_local_pairinv.py')),
            'optimized local registration, target or binary changed')
    filenames = {key: key.replace('_','-')+'.json' for key in (
        'panel', 'method', 'candidate', 'ic_workload', 'rho_workload',
        'ic_job', 'rho_job', 'build_record', 'dependency_manifest')}
    require(seal['files'] == {key: digest(REGISTRATION/name)
                              for key, name in filenames.items()}
            and seal['inventory_files'] == {name: digest(REGISTRATION/'inventory'/name)
                for name in ('job.json', 'stdout.json', 'stderr.txt',
                             'process.json')},
            'optimized local manifest or inventory bytes differ from seal')
    ic_job, rho_job = jobs()
    require(read(REGISTRATION/'ic-job.json') == ic_job
            and read(REGISTRATION/'rho-job.json') == rho_job
            and read(REGISTRATION/'build-record.json') == build,
            'optimized local job or build changed')
    inventory = read(REGISTRATION/'inventory/stdout.json')
    require(inventory['status'] == 'inventory'
            and inventory['fixture']['targets'] == ic_job['public_targets']
            and inventory['fixture']['target_seeds'] == ic_job['target_seeds']
            and read(REGISTRATION/'inventory/job.json')
                == (ic_job | {'mode':'inventory'})
            and panel['inventory_job_sha256']
                == digest(REGISTRATION/'inventory/job.json')
            and panel['inventory_report_sha256']
                == digest(REGISTRATION/'inventory/stdout.json')
            and panel['inventory_process_sha256']
                == digest(REGISTRATION/'inventory/process.json'),
            'optimized local factor-base inventory changed')
    check_build_identity(inventory, IC_SOURCE, panel['field_kernel'])
    return dict(source=source, build=build, panel=panel,
                inventory=inventory, job=ic_job if arm == 'ic' else rho_job,
                method=read(REGISTRATION/'method.json'),
                candidate=read(REGISTRATION/'candidate.json'))


def audit_ic(report, process, data):
    certificate = verify(report, data['inventory']['fixture'], summands=3)
    stages = audit_stages(report, data['inventory']['fixture'],
                          data['job']['algorithm_seed'])
    require(candidate_manifest(report['fixture'], report, data['method'])
                == data['candidate'],
            'measured pairinv candidate differs from frozen inventory/method')
    timing = native_intervals(report, process['process_wall_ns'])
    require(certificate['verified_targets'] == 1
            and stages['queries']['verified_relations'] > 0,
            'optimized incumbent lacks complete IC certificate')
    return dict(certificate=certificate, stages=stages,
                timing=timing, online_wall_ns=timing['online']['wall_ns'])


def audit_rho(report, process, data):
    certificate = verify(report, data['inventory']['fixture'],
                         expected_mode='rho')
    timing = native_intervals(report, process['process_wall_ns'])
    require(report['executed_method'] == {
                'reference':'signed_frobenius_rho', 'requested_walks':1}
            and len(report['solutions']) == 1
            and report['solutions'][0]['effective_walks'] == 1
            and report['automorphism_order']
                == 2 * int(report['fixture']['degree'])
            and certificate['verified_targets'] == 1,
            'rho walk policy or one-target certificate changed')
    return dict(certificate=certificate,
                timing=timing, online_wall_ns=timing['online']['wall_ns'],
                rho_dispatch=report['executed_method'],
                iterations=report['solutions'][0]['iterations'],
                restarts=report['solutions'][0]['restarts'])


def run(prepared, built, arm, out):
    data = sealed_inputs(prepared, built, arm)
    require(not out.exists(), 'one-shot optimized arm output already exists')
    out.mkdir(parents=True)
    for name in ('panel.json', 'seal.json', 'method.json',
                 'candidate.json', 'ic-workload.json', 'rho-workload.json',
                 'ic-job.json', 'rho-job.json', 'build-record.json',
                 'dependency-manifest.json',
                 'PROTOCOL.md'):
        shutil.copy2(REGISTRATION/name, out/name)
    shutil.copy2(Path(built)/'build.log', out/'build.log')
    require(digest(out/'build.log') == data['build']['build_log_sha256'],
            'copied optimized reference build log changed')
    shutil.copytree(REGISTRATION/'inventory', out/'inventory')
    shutil.copy2(Path(__file__), out/'registered-runner.py')
    shutil.copy2(Path(__file__).with_name('register_local_pairinv.py'),
                 out/'registered-registration.py')
    worker = out/'worker'
    shutil.copy2(Path(built)/'worker', worker)
    require(digest(worker) == data['panel']['worker_sha256'],
            'copied optimized reference executable changed')
    write(out/'host.json', dict(system=platform.system(),
          machine=platform.machine(), processor=platform.processor(),
          cpu_count=os.cpu_count(),
          scope='local single-target diagnostic; hardware calibration not claimed'),
          exclusive=True)
    job = data['job']
    write(out/'executed-job.json', job, exclusive=True)
    process = run_worker(worker, job, out)
    result = dict(schema_version=1, arm=arm,
                  status='PRODUCER_FAILURE', process=process,
                  candidate_id=(data['panel']['candidate_id'] if arm == 'ic'
                                else None),
                  reference_id=(data['panel']['rho_reference_id']
                                if arm == 'rho' else None),
                  workload_id=data['panel'][arm+'_workload_id'],
                  run_id=data['panel'][arm+'_run_id'],
                  target_point=read(TARGET)['target_input']['point'],
                  verified=False, online_wall_ns=None,
                  paired_rho_speedup=None, promotion_eligible=False)
    try:
        report = json.loads((out/'stdout.json').read_text())
        result['raw_report_status'] = report.get('status')
        require(report['fixture'] == data['inventory']['fixture'],
                'worker substituted public target or subgroup')
        check_build_identity(report, IC_SOURCE, data['panel']['field_kernel'])
        if report['status'] == 'incomplete' and process['exit_code'] == 2:
            result['status'] = 'BOUNDED_INCOMPLETE'
        else:
            require(process['exit_code'] == 0
                    and report['status'] == 'complete',
                    'worker exit or report is not complete')
            admitted = (audit_ic(report, process, data) if arm == 'ic'
                        else audit_rho(report, process, data))
            write(out/'admission.json', admitted, exclusive=True)
            result.update(status='AUDITED_COMPLETE', verified=True,
                          online_wall_ns=admitted['online_wall_ns'])
    except Exception as error:
        result.update(status='AUDIT_OR_PRODUCER_FAILURE',
                      audit_error=f'{type(error).__name__}: {error}')
    write(out/'result.json', result, exclusive=True)
    print(json.dumps(result, sort_keys=True), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--arm', choices=('ic', 'rho'), required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(args.prepared, args.build, args.arm, args.out.resolve())


if __name__ == '__main__':
    main()
