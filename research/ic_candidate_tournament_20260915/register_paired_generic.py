"""Seal the archived source-bound F5 and rho arms on one fresh public point."""
import copy
import hashlib
import io
from pathlib import Path
import tarfile

from generic_admission import method_record
from generic_bases import verify_base
from generic_build import verify_build_record
from identity import (candidate_manifest, run_id, sha256,
                      workload_manifest, write_immutable)
from oracle import Curve, require
from run_generic_exact_yield_audit import (PANEL as PARENT_PANEL,
                                           load_evidence, record)
from tournament import read

HERE = Path(__file__).resolve().parent
REGISTRATION = HERE/'goal_20260924/paired-fresh-n17a1/generic-arms'
TARGET = REGISTRATION.parent/'target-panel.json'
PARENT_ARCHIVE = (HERE/'goal_20260924/generic-backend-recovery-pilot/'
                  'evidence.tar.gz')
PARENT_ARCHIVE_SHA256 = 'e0ce19fc28c58e2dc1cae9649a16af74099016fff8183e5ce69f65b39c804f02'
SOURCE_COMMIT = '765c3c5f19032bd852163805f257c56babef2040'
RESOURCES = dict(host_class='physical-macos-arm64', cpu_workers=1,
                 target_count=1, memory_limit_bytes=None,
                 total_wall_limit_seconds=None)
INPUT_LAW = 'one-supplied-public-point; seed-is-provenance'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inputs():
    require(digest(PARENT_ARCHIVE) == PARENT_ARCHIVE_SHA256,
            'source-bound F5 archive changed')
    files = load_evidence(read(PARENT_PANEL))
    build = record(files, 'build/build-record.json')
    source = record(files, 'build/source-manifest.json')
    verify_build_record(build, source)
    require(hashlib.sha256(files['build/worker']).hexdigest()
                == build['worker_sha256'],
            'retained source-bound worker changed')
    with tarfile.open(fileobj=io.BytesIO(files['build/root-source.tar.gz']),
                      mode='r:gz') as archive:
        require(set(archive.getnames()) == set(source['root_files']),
                'retained Rust source snapshot member set changed')
        for name, expected in source['root_files'].items():
            require(hashlib.sha256(archive.extractfile(name).read()).hexdigest()
                        == expected,
                    'retained Rust source snapshot differs: '+name)
    old = record(files, 'jobs/n17a1/f5/stdout.json')
    curve = Curve(old['fixture'])
    target_panel = read(TARGET)
    point = curve.decode(target_panel['target_input']['point'])
    require(target_panel['status'] == 'TARGET_FROZEN_BEFORE_ARM_OUTCOMES'
            and point is not None and curve.mul(point, curve.r) is None,
            'fresh paired point is not in the intended subgroup')
    require(old['fixture']['degree'] == 17
            and old['fixture']['curve_a'] == 1
            and len(old['factor_base']) == 63
            and old['columns'] == 29,
            'archived F5 reference differs from registered n17a1 geometry')
    return files, build, source, old, target_panel


def jobs(old, target_panel):
    target = [str(value) for value in target_panel['target_input']['point']]
    common = dict(degree=17, curve_a=1,
                  public_targets=[target],
                  target_seeds=[target_panel['target_input']['seed']],
                  factor_base=dict(kind='standard_subspace', dimension=6),
                  exclusive_phases=True)
    f5 = dict(common, mode='ic', algorithm_seed=2026092955,
              config=dict(solver='f5', linear_algebra='dense', summands=3,
                          groebner_degree=3, node_budget=4096,
                          conflict_budget=100000, batch_trials=8,
                          max_trials=256))
    rho = dict(common, mode='rho', algorithm_seed=2026092958,
               config=dict(solver='pair_table', linear_algebra='dense',
                           summands=3, rho_parallel_walks=1,
                           max_trials=65536))
    return f5, rho


def identities():
    files, build, source, old, target_panel = inputs()
    f5_job, rho_job = jobs(old, target_panel)
    # Inventory and observed F5 dispatch come from the audited complete
    # source-bound pilot. Only registered limit fields change; the real arm
    # must independently report and match this manifest after execution.
    pre_run = copy.deepcopy(old)
    pre_run['effective_config']['max_trials'] = 256
    base = verify_base(old, old['fixture'], f5_job)
    method = method_record(f5_job, pre_run, {'base':base}, build)
    candidate = candidate_manifest(old['fixture'], old, method)
    fixture = copy.deepcopy(old['fixture'])
    fixture.update(targets=[list(map(int, f5_job['public_targets'][0]))],
                   target_seeds=[target_panel['target_input']['seed']],
                   target_scalar_constructed=False)
    f5_workload = workload_manifest(
        fixture, input_law=INPUT_LAW,
        algorithm_seed=f5_job['algorithm_seed'],
        resource_envelope=RESOURCES, cache_policy='cold')
    rho_workload = workload_manifest(
        fixture, input_law=INPUT_LAW,
        algorithm_seed=rho_job['algorithm_seed'],
        resource_envelope=RESOURCES, cache_policy='cold')
    rho_id = 'RHO1N17Ckb1h'+sha256(dict(
        curve_id=candidate['record']['curve']['curve_id'],
        source_manifest_sha256=build['source_manifest_sha256'],
        worker_sha256=build['worker_sha256'],
        config=rho_job['config']))[:12]
    panel = dict(schema_version=1,
                 status='REGISTERED_BEFORE_ARM_EXECUTION',
                 purpose='same-fresh-point-F5-and-rho-source-bound-arms',
                 target_panel_sha256=digest(TARGET),
                 parent_evidence_sha256=PARENT_ARCHIVE_SHA256,
                 source_commit=SOURCE_COMMIT,
                 source_manifest_sha256=build['source_manifest_sha256'],
                 build_record_sha256=sha256(build),
                 worker_sha256=build['worker_sha256'],
                 candidate_id=candidate['candidate_id'],
                 f5_workload_id=f5_workload['workload_id'],
                 f5_run_id=run_id(candidate['candidate_id'],
                                  f5_workload['workload_id'], 0),
                 rho_reference_id=rho_id,
                 rho_workload_id=rho_workload['workload_id'],
                 rho_run_id=f'{rho_id}W{rho_workload["workload_id"]}R0',
                 resource_envelope=RESOURCES,
                 f5_attempt_cap=256, rho_restart_iteration_cap=65536,
                 worker_threads=1,
                 scheduling='execute each sealed arm once; preserve failures; '
                            'no contingent target substitution',
                 claim_boundary='single-arm admission only until the full '
                                'target panel, same-point reference and '
                                'host-calibration audit are complete')
    return dict(files=files, build=build, source=source, old=old,
                target_panel=target_panel, f5_job=f5_job, rho_job=rho_job,
                method=method, candidate=candidate,
                f5_workload=f5_workload, rho_workload=rho_workload,
                panel=panel)


def register():
    data = identities()
    require(not (REGISTRATION/'panel.json').exists(),
            'generic paired arms already registered')
    REGISTRATION.mkdir(parents=True, exist_ok=True)
    for name, value in (
            ('panel.json', data['panel']),
            ('method.json', data['method']),
            ('candidate.json', data['candidate']),
            ('f5-workload.json', data['f5_workload']),
            ('rho-workload.json', data['rho_workload']),
            ('f5-job.json', data['f5_job']),
            ('rho-job.json', data['rho_job'])):
        write_immutable(REGISTRATION/name, value)
    seal = dict(schema_version=1,
                files={name: digest(REGISTRATION/name) for name in (
                    'panel.json', 'method.json', 'candidate.json',
                    'f5-workload.json', 'rho-workload.json',
                    'f5-job.json', 'rho-job.json')},
                controller_sha256=digest(HERE/'run_paired_generic.py'),
                registration_sha256=digest(Path(__file__)))
    write_immutable(REGISTRATION/'seal.json', seal)
    return data['panel']


if __name__ == '__main__':
    print(register())
