"""Generic scientific admission for the existing qualification driver.

Inventory declares a source-bound method without solving the supplied target.
Measured runs must independently verify that declaration and all executed stages.
No qualification or promotion decision is made by this adapter.
"""
import copy
import json
from pathlib import Path
import shutil
import tarfile

from generic_admission import admit_rho, method_record, scientific_ledger
from generic_build import digest, verify_binding, verify_build_record
from generic_phases import verify_native
from generic_stages import STRATEGIES, effective_config, kernel, verify_stages
from identity import candidate_manifest, natural, run_id, sha256, workload_manifest
from measurement import PHASES, exclusive_ledger, measured_run, report_sha256
from oracle import require

ADAPTER = 'generic-v1'
INVENTORY_POLICY = 'configuration-only-no-query-or-target-solve-v1'


def snapshot(build_dir, destination, compiler):
    """Transport a controlled build, its exact source, and its retained receipts."""
    from tournament import write
    build_dir = Path(build_dir).resolve()
    source = json.loads((build_dir/'source-manifest.json').read_text())
    build = json.loads((build_dir/'build-record.json').read_text())
    verify_build_record(build, source)
    require(build['build']['target_arch'] == 'x86_64' and build['build']['target_os'] == 'linux'
            and build['build']['rustc'].splitlines()[0] == compiler,
            'generic qualification needs the current Linux amd64 compiler build')
    require(digest(build_dir/'worker') == build['worker_sha256'], 'changed generic executable')
    for name in ('worker', 'source-manifest.json', 'build-record.json', 'root-source.tar.gz', 'build.log'):
        shutil.copy2(build_dir/name, destination/name)
    root = destination/'source'
    root.mkdir()
    with tarfile.open(destination/'root-source.tar.gz') as archive:
        archive.extractall(root, filter='data')
    actual = {p.relative_to(root).as_posix(): digest(p) for p in root.rglob('*') if p.is_file()}
    require(sha256(actual) == sha256(source['root_files']), 'generic root source archive differs')
    metadata = dict(adapter=ADAPTER, build_record=build)
    write(destination/'producer.json', metadata, exclusive=True)
    write(destination/'preparation.json', dict(adapter=ADAPTER,
        status='SOURCE_BOUND_NOT_PERFORMANCE_QUALIFIED', source_manifest_sha256=sha256(source)), exclusive=True)
    return destination/'worker', source


def make_admission(*, job, fixture, report, manifest, metadata, resources,
                   worker_sha256, executable):
    from driver_admission import INPUT_LAW
    require(metadata.get('adapter') == ADAPTER and job.get('exclusive_phases') is True,
            'generic adapter requires declared exclusive measurements')
    build = metadata['build_record']
    verify_binding(report, build, manifest, executable=executable)
    require(worker_sha256 == build['worker_sha256'], 'changed declared generic binary')
    require(report['status'] == 'inventory' and report['mode'] == job['mode']
            and report['fixture'] == fixture and job['public_targets'] == fixture['targets']
            and len(fixture['targets']) == 1, 'generic inventory changed mode or fixture')
    require(report.get('inventory_policy') == INVENTORY_POLICY
            and report.get('generic_runtime_policy') == 'default-environment-one-rayon-v1'
            and type(report.get('generic_admission_schema')) is int
            and report['generic_admission_schema'] == 1,
            'missing generic inventory policy')
    require(report['online_wall_ns'] is None and report['scalar_replay_included'] is False
            and not report.get('solutions') and not report.get('collection_reports'),
            'inventory performed target work')
    cfg = effective_config(job, report)
    field = kernel(report['field_kernel'])
    workload = workload_manifest(fixture, input_law=INPUT_LAW,
        algorithm_seed=job['algorithm_seed'], resource_envelope=resources)
    common = dict(adapter=ADAPTER, job=copy.deepcopy(job), fixture=fixture, report=report,
        workload=workload, source_manifest_sha256=sha256(manifest), worker_sha256=worker_sha256,
        mode=job['mode'], metadata=metadata, field_kernel=field)
    if job['mode'] == 'ic':
        stages = verify_stages(report, fixture, dict(job, mode='inventory'))
        pair = cfg['solver'] == 'pair_table'
        window = cfg['collection_window']
        active = pair and cfg['summands'] == 3 and window is not None and 0 < window < len(report['factor_base'])
        plan = dict(strategy=STRATEGIES[cfg['solver']], field_kernel=field, pair_table=pair,
            query_rule='windowed-walk-64' if active else 'trial-keyed-sample',
            collection_window=window if active else None)
        method = method_record(job, report, stages, build, collector_plan=plan)
        return dict(common, method=method, candidate=candidate_manifest(fixture, report, method),
                    collector_plan=plan)
    require(job['mode'] == 'rho', 'unknown generic algorithm mode')
    reference = dict(curve_id=workload['record']['curve_id'], algorithm='signed-Frobenius-rho',
        configuration=dict(requested_walks=cfg['rho_parallel_walks'], max_iterations=cfg['max_trials'],
                           jump_count=16, max_restarts=64, progress_interval=256),
        source_manifest_sha256=sha256(manifest), build_sha256=build['build_sha256'], field_kernel=field)
    return dict(common, reference=dict(reference_id='RHO1h'+sha256(reference)[:12],
        record_sha256=sha256(reference), record=reference))


def native_timing(report, job, process_wall_ns):
    clocks = verify_native(report, job, process_wall_ns=process_wall_ns)
    require(clocks['complete_phase_coverage'], 'incomplete generic native interval')
    cold = (scientific_ledger(clocks['process_phases_ns'], unit='native_wall_ns',
                             process_total=process_wall_ns)['operations']
            if job['mode'] == 'ic' else
            {k:v for k,v in clocks['process_phases_ns'].items() if v is not None})
    return dict(schema_version=1, unit='native_monotonic_ns', target_count=1,
        native_report_sha256=report_sha256(report),
        online=dict(phase_wall_ns={k:v for k,v in clocks['online_phases_ns'].items() if v is not None},
            wall_ns=clocks['online_wall_ns'], target_generation_included=False, scalar_replay_included=True,
            boundary='one supplied point; after reusable preparation through scalar replay'),
        cold=dict(phase_wall_ns=cold, wall_ns=process_wall_ns,
            worker_snapshot_wall_ns=clocks['observed_wall_ns'],
            external_setup_remainder_ns=clocks['external_setup_remainder_ns'],
            remainder_policy='setup includes process/input and report/exit tail', unattributed_wall_ns=0))


def run_record(admitted, *, manifest, executable, number, host_id, status,
               native=None, process_wall_ns=None, profile=None, costs=None, profile_wall_ns=None,
               native_status=None, profile_status=None):
    natural(number, 'generic native run number')
    require(status in ('complete', 'timeout', 'oom', 'error'), 'unknown generic run status')
    job, fixture = admitted['job'], admitted['fixture']
    build = admitted['metadata']['build_record']
    proof = timing = stages = None
    ledger = exclusive_ledger(dict.fromkeys(PHASES), unit='valgrind-3.22-amd64-Ir',
                             process_operations=None, zero_reasons={})
    if status == 'complete':
        require(native is not None and profile is not None and costs is not None,
                'generic qualification requires a native/profile pair')
        natural(profile_wall_ns, 'profiled process nanoseconds', positive=True)
        observations = []
        for report in (native, profile):
            require(report['fixture'] == fixture and report['status'] == 'complete',
                    'generic measured fixture or completion differs')
            verify_binding(report, build, manifest, executable=executable)
            verify_native(report, job, process_wall_ns=process_wall_ns if report is native else profile_wall_ns)
            if job['mode'] == 'ic':
                audit = verify_stages(report, fixture, job)
                method = method_record(job, report, audit, build)
                require(sha256(method) == sha256(admitted['method']) and
                        candidate_manifest(fixture, report, method) == admitted['candidate'],
                        'generic measured method differs from inventory')
                observations.append(audit)
                proof = audit['certificate']
                stages = audit
            else:
                audit = admit_rho(report, fixture, job, build, manifest, executable=executable,
                    process_wall_ns=process_wall_ns if report is native else profile_wall_ns)
                require(report['field_kernel'] == admitted['field_kernel'], 'generic rho kernel differs')
                observations.append(dict(certificate=audit['certificate'], dispatch=audit['dispatch']))
                proof = audit['certificate']
        require(report_sha256(observations[0]) == report_sha256(observations[1]),
                'generic native/profile stage evidence differs')
        timing = native_timing(native, job, process_wall_ns)
        total = sum(v for v in costs.values() if v is not None)
        if job['mode'] == 'ic':
            ledger = scientific_ledger(costs, unit='valgrind-3.22-amd64-Ir', process_total=total)
            require(ledger['complete'], 'complete generic result has unpriced scientific phases')
        else:
            require(set(k for k,v in costs.items() if v is not None)
                    == {'setup', 'precompute', 'rho_solve', 'recovery_check'}, 'invalid generic rho phases')
            require(all(type(v) is int and v >= 0 for v in costs.values() if v is not None)
                    and total > 0, 'invalid generic rho instruction costs')
    else:
        require(native is None and profile is None and costs is None,
                'failed generic execution cannot certify a measured result')
        total = None
    workload = admitted['workload']
    provenance = dict(source_manifest_sha256=admitted['source_manifest_sha256'],
        worker_sha256=admitted['worker_sha256'], host_id=host_id,
        resource_envelope_id=sha256(workload['record']['resource_envelope']),
        calibration_id='valgrind-3.22-amd64-Ir',
        report_sha256=report_sha256(profile) if status == 'complete' else None)
    if job['mode'] == 'ic':
        record = measured_run(candidate=admitted['candidate'], workload=workload, report=profile,
            fixture=fixture, method=admitted['method'], number=number, ledger=ledger,
            native_wall_ns=process_wall_ns, status=status, provenance=provenance,
            admission_report=admitted['report'])
        identity = dict(candidate_id=admitted['candidate']['candidate_id'])
        profile_id = run_id(identity['candidate_id'], workload['workload_id'], number+1)
        record['stage_audit'] = stages
    else:
        reference_id = admitted['reference']['reference_id']
        identity = dict(reference_id=reference_id)
        record = dict(schema_version=2, **identity, workload_id=workload['workload_id'],
            run_id=f"{reference_id}W{workload['workload_id']}R{number}", status=status,
            instruction_phases=costs, total_operations=total, native_wall_ns=process_wall_ns,
            certificate=proof, provenance=provenance, promotion_eligible=False)
        profile_id = f"{reference_id}W{workload['workload_id']}R{number+1}"
    record['native_timing'] = timing
    record['adapter'] = ADAPTER
    record['native_process_status'] = native_status
    record['profile_execution'] = dict(**identity, workload_id=workload['workload_id'],
        run_id=profile_id, status=status, process_status=profile_status, process_wall_ns=profile_wall_ns,
        total_operations=total, certificate=proof, provenance=provenance,
        unit='valgrind-3.22-amd64-Ir', timing_scope='profiled process; never native timing',
        promotion_eligible=False)
    return record
