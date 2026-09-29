"""Seal the accepted pairinv source as local incumbent and same-point rho."""
import argparse
import hashlib
import json
from pathlib import Path

from campaign_rules import IC_SOURCE
from identity import (candidate_manifest, run_id, sha256,
                      workload_manifest, write_immutable)
from oracle import Curve, require
from producer.evidence import (check_build_identity, method_record)
from producer.prepare import verify_source
from register_paired_generic import RESOURCES, digest
from run_paired_generic import run_worker
from tournament import read

HERE = Path(__file__).resolve().parent
REGISTRATION = HERE/'goal_20260924/paired-fresh-n17a1/pairinv-local'
TARGET = REGISTRATION.parent/'target-panel.json'
INPUT_LAW = 'one-supplied-public-point; seed-is-provenance'


def jobs():
    target = read(TARGET)['target_input']
    point = [str(value) for value in target['point']]
    common = dict(degree=17, curve_a=1, public_targets=[point],
                  target_seeds=[target['seed']],
                  factor_base=dict(kind='subgroup_orbits', seed=43, points=102))
    ic = dict(common, mode='ic', algorithm_seed=2026092960,
              config=dict(solver='pair_table', linear_algebra='tiny_gauss',
                          summands=3, batch_trials=1, max_trials=65536))
    rho = dict(common, mode='rho', algorithm_seed=2026092961,
               config=dict(solver='pair_table', linear_algebra='tiny_gauss',
                           summands=3, batch_trials=1, max_trials=65536,
                           rho_parallel_walks=1))
    return ic, rho


def build_inputs(prepared, built):
    prepared, built = Path(prepared).resolve(), Path(built).resolve()
    receipt = read(prepared/'preparation.json')
    require(receipt['reference'] == 'pairinv'
            and receipt['source_manifest_sha256'] == IC_SOURCE,
            'prepared source is not accepted pairinv')
    source = verify_source(prepared, IC_SOURCE)
    record = read(built/'build-record.json')
    dependencies = read(built/'dependency-manifest.json')
    require(record['accepted_source_manifest_sha256'] == IC_SOURCE
            and record['source_manifest'] == source
            and record['prepared_source_receipt'] == receipt
            and read(built/'build-policy.json') == record['build_policy']
            and read(built/'build-exit.json') == {'exit_code': 0}
            and record['dependency_manifest_sha256'] == sha256(dependencies)
            and record['build_policy']['dependency_manifest_sha256']
                == sha256(dependencies)
            and record['worker_sha256'] == digest(built/'worker')
            and record['build_policy_sha256']
                == sha256(record['build_policy'])
            and record['build_log_sha256'] == digest(built/'build.log')
            and record['builder_sha256']
                == digest(HERE/'build_local_pairinv.py'),
            'local pairinv build receipt differs from retained source/binary')
    return source, record


def register(prepared, built):
    source, build = build_inputs(prepared, built)
    require(not (REGISTRATION/'panel.json').exists(),
            'local incumbent registration already exists')
    REGISTRATION.mkdir(parents=True, exist_ok=True)
    ic_job, rho_job = jobs()
    inventory = REGISTRATION/'inventory'
    inventory.mkdir(exist_ok=False)
    inventory_job = ic_job | {'mode':'inventory'}
    write_immutable(inventory/'job.json', inventory_job)
    process = run_worker(Path(built).resolve()/'worker', inventory_job,
                         inventory)
    require(process['exit_code'] == 0, 'pre-execution incumbent inventory failed')
    report = json.loads((inventory/'stdout.json').read_text())
    require(report['status'] == 'inventory'
            and report['fixture']['targets'] == ic_job['public_targets']
            and report['fixture']['target_seeds'] == ic_job['target_seeds']
            and report['fixture']['target_scalar_constructed'] is False,
            'inventory did not use the frozen public point')
    check_build_identity(report, IC_SOURCE)
    curve = Curve(report['fixture'])
    target = curve.decode(ic_job['public_targets'][0])
    require(target is not None and curve.mul(target, curve.r) is None,
            'local incumbent target leaves subgroup')
    field = report['field_kernel']
    require(field == 'portable',
            'local ARM64 reference selected an unexpected field backend')
    build_flags = dict(compiler=build['build_policy']['compiler'],
                       cargo=build['build_policy']['cargo'],
                       target=build['build_policy']['target'],
                       cargo_config_sha256=build['build_policy']['cargo_config_sha256'],
                       cargo_overrides=build['build_policy']['cargo_overrides'],
                       dependency_manifest_sha256=build['dependency_manifest_sha256'],
                       build_policy_sha256=build['build_policy_sha256'],
                       worker_sha256=build['worker_sha256'],
                       field_kernel=field)
    method = method_record(ic_job, report['fixture'], source, IC_SOURCE,
                           'pairinv', build_flags)
    candidate = candidate_manifest(report['fixture'], report, method)
    ic_workload = workload_manifest(
        report['fixture'], input_law=INPUT_LAW,
        algorithm_seed=ic_job['algorithm_seed'],
        resource_envelope=RESOURCES, cache_policy='cold')
    rho_workload = workload_manifest(
        report['fixture'], input_law=INPUT_LAW,
        algorithm_seed=rho_job['algorithm_seed'],
        resource_envelope=RESOURCES, cache_policy='cold')
    rho_id = 'RHO1N17Ckb1h'+sha256(dict(
        curve_id=candidate['record']['curve']['curve_id'],
        source_manifest_sha256=IC_SOURCE,
        worker_sha256=build['worker_sha256'],
        config=rho_job['config']))[:12]
    panel = dict(schema_version=1,
                 status='REGISTERED_BEFORE_ARM_EXECUTION',
                 target_panel_sha256=digest(TARGET),
                 inventory_job_sha256=digest(inventory/'job.json'),
                 inventory_report_sha256=digest(inventory/'stdout.json'),
                 inventory_process_sha256=digest(inventory/'process.json'),
                 accepted_source_manifest_sha256=IC_SOURCE,
                 source_preparation_sha256=sha256(build['prepared_source_receipt']),
                 build_record_sha256=sha256(build),
                 worker_sha256=build['worker_sha256'],
                 field_kernel=field,
                 candidate_id=candidate['candidate_id'],
                 ic_workload_id=ic_workload['workload_id'],
                 ic_run_id=run_id(candidate['candidate_id'],
                                  ic_workload['workload_id'], 0),
                 rho_reference_id=rho_id,
                 rho_workload_id=rho_workload['workload_id'],
                 rho_run_id=f'{rho_id}W{rho_workload["workload_id"]}R0',
                 resource_envelope=RESOURCES,
                 incumbent_factor_base_policy='sampled signed Frobenius orbits; '
                                              'different from standard subspace',
                 comparison_kind='complete-pipeline-and-factor-base-policy',
                 headline_speedup=None,
                 claim_boundary='same-point single-arm diagnostics until all '
                                'four arms, rho quality and host timing are audited')
    values = dict(panel=panel, method=method, candidate=candidate,
                  ic_workload=ic_workload, rho_workload=rho_workload,
                  ic_job=ic_job, rho_job=rho_job,
                  build_record=build,
                  dependency_manifest=read(Path(built)/'dependency-manifest.json'))
    for key, value in values.items():
        write_immutable(REGISTRATION/(key.replace('_','-')+'.json'), value)
    seal = dict(schema_version=1,
                files={key: digest(REGISTRATION/(key.replace('_','-')+'.json'))
                       for key in values},
                inventory_files={name: digest(inventory/name) for name in
                                 ('job.json', 'stdout.json', 'stderr.txt',
                                  'process.json')},
                controller_sha256=digest(HERE/'run_local_pairinv.py'),
                registration_sha256=digest(Path(__file__)))
    write_immutable(REGISTRATION/'seal.json', seal)
    return panel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--build', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(register(args.prepared, args.build), sort_keys=True))


if __name__ == '__main__':
    main()
