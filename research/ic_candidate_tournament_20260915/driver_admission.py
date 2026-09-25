"""Shared admission and run records for the existing native and instruction drivers.

Only the audited optimized schema-3 producers are admitted here. Adding another
backend requires a method adapter and independent stage audit, not a solver tag.
"""
import copy
import json
from pathlib import Path
import statistics
import tomllib

from identity import candidate_manifest, sha256, workload_manifest, write_immutable
from measurement import PHASES, exclusive_ledger, measured_run, report_sha256
from oracle import require, verify
from producer.evidence import audit_stages, check_build_identity, method_record, scientific_ledger, executed_policy
from producer.timing import native_intervals

EVALUATOR = ('autolab.py', 'tournament.py', 'portfolio.py', 'oracle.py', 'identity.py',
             'campaign_rules.py', 'target_history.py',
             'measurement.py', 'driver_admission.py', 'qualification.py', 'producer/evidence.py', 'producer/timing.py')
INPUT_LAW = ('public-hash-to-curve-cofactor-v1; independently generated fixture, '
             'one supplied public point, no planted scalar')


def load(path):
    return json.loads(Path(path).read_text())


def producer_metadata(source, manifest, compiler):
    """Bind the prepared source policy; arbitrary aliases cannot name a backend."""
    preparation = load(source.parent/'preparation.json')
    require(preparation['instrumented'] is True and
            preparation['source_manifest_sha256'] == sha256(manifest),
            'prepared producer source manifest differs')
    require(preparation['reference'] in ('both', 'scaled', 'pairinv'), 'unknown prepared source policy')
    require(manifest.get('src/cryptanalysis/ic_phase.rs') == preparation['phase_module_sha256'],
            'prepared timing module differs')
    config = tomllib.loads((source/'.cargo/config.toml').read_text())
    target = config.get('build', {}).get('target')
    require(isinstance(target, str) and target, 'prepared build must declare one target')
    return {'preparation': preparation, 'build': {'compiler': compiler, 'target': target,
            'cargo_config_sha256': manifest['.cargo/config.toml']}}


def make_admission(*, job, fixture, report, manifest, metadata, resources, worker_sha256):
    require(len(fixture['targets']) == 1 and job.get('public_targets') == fixture['targets'],
            'admission requires one frozen supplied public point')
    require(report.get('status') == 'inventory' and report['fixture'] == fixture,
            'inventory changed the prepared fixture')
    require(report.get('phase_schema') == 3 and report.get('target_input') == 'supplied_public_point',
            'inventory lacks public-point scientific admission')
    source_sha = sha256(manifest)
    require(metadata['preparation']['source_manifest_sha256'] == source_sha, 'changed prepared source')
    check_build_identity(report, source_sha)
    kernel = report['field_kernel']
    require(kernel in ('portable', 'pclmulqdq'), 'unknown admitted arithmetic kernel')
    build = dict(metadata['build'], field_kernel=kernel)
    workload = workload_manifest(fixture, input_law=INPUT_LAW,
        algorithm_seed=job['algorithm_seed'], resource_envelope=resources)
    common = dict(job=copy.deepcopy(job), fixture=fixture, report=report, workload=workload,
                  source_manifest_sha256=source_sha, worker_sha256=worker_sha256,
                  build=build, mode=job['mode'])
    if job['mode'] == 'ic':
        policy = executed_policy(job['config'], metadata['preparation'].get('candidate_panel'))
        require(report.get('implementation_policy') == policy, 'inventory executed a different candidate policy')
        method = method_record(job, fixture, manifest, source_sha,
                               metadata['preparation']['reference'], build,
                               metadata['preparation'].get('candidate_panel'))
        candidate = candidate_manifest(fixture, report, method)
        return dict(common, candidate=candidate, method=method, implementation_policy=policy)
    require(job['mode'] == 'rho', 'unknown admitted algorithm')
    # Rho is a reference, never an IC candidate with fictitious PDP/LA stages.
    reference = {'curve_id': workload['record']['curve_id'], 'algorithm': 'signed-Frobenius-rho',
                 'configuration': job['config'], 'source_manifest_sha256': source_sha,
                 'build': build}
    return dict(common, reference={'reference_id': 'RHO1h'+sha256(reference)[:12],
                                  'record_sha256': sha256(reference), 'record': reference})


def freeze_admission(directory, *, binary, job, fixture, manifest, metadata, resources,
                     worker_sha256, execute):
    inventory_job = dict(job, mode='inventory')
    process = execute(binary, inventory_job, directory)
    if (directory/'job.json').exists():
        require(load(directory/'job.json') == inventory_job, 'changed inventory job')
    else:
        write_immutable(directory/'job.json', inventory_job)
    # Process durations are measurements and may be floats. They never enter
    # the canonical method hash, whose JSON representation forbids floats.
    with (directory/'process.json').open('x') as stream:
        json.dump(process, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')
    require(process['exit_code'] == 0 and
            process.get('process_status', process.get('status')) == 'EXITED',
            'inventory admission failed; raw failure retained')
    admitted = make_admission(job=job, fixture=fixture, report=load(directory/'stdout.json'),
        manifest=manifest, metadata=metadata, resources=resources, worker_sha256=worker_sha256)
    for name in ('workload', 'candidate', 'method', 'reference'):
        if name in admitted:
            write_immutable(directory/(name+'.json'), admitted[name])
    write_immutable(directory/'admission.json', admitted)
    return admitted


def check_admission(admitted, *, job, fixture, manifest, metadata, resources, worker_sha256):
    expected = make_admission(job=job, fixture=fixture, report=admitted['report'], manifest=manifest,
        metadata=metadata, resources=resources, worker_sha256=worker_sha256)
    require(admitted == expected, 'changed canonical admission')


def distinct_candidates(named_admissions):
    """Ignored configuration flags cannot create another measured IC method."""
    seen = set()
    for alias, admitted in named_admissions:
        if alias == 'aa_control' or admitted['mode'] != 'ic':
            continue
        identity = admitted['candidate']['record_sha256']
        require(identity not in seen, 'duplicate admitted IC method under different aliases')
        seen.add(identity)


def run_record(admitted, *, number, host_id, status, native=None, process_wall_ns=None,
               profile=None, costs=None):
    """Reconstructible success/failure record; never rank a failed or partial solve."""
    fixture, job = admitted['fixture'], admitted['job']
    proof = timing = stage_audit = None
    ledger = exclusive_ledger(dict.fromkeys(PHASES), unit='valgrind-3.22-amd64-Ir',
                             process_operations=None, zero_reasons={})
    require(status in ('complete', 'timeout', 'oom', 'error'), 'unknown driver outcome')
    if status == 'complete':
        require(native is not None, 'complete run requires native replay')
        reports = [native] + ([profile] if profile is not None else [])
        for report in reports:
            check_build_identity(report, admitted['source_manifest_sha256'], admitted['build']['field_kernel'])
            checked = verify(report, fixture, expected_mode=job['mode'], summands=job['config']['summands'])
            require(not checked.get('degenerate_descents', 0), 'direct collision cannot enter IC')
            require(proof is None or checked == proof, 'native/profile certificate differs')
            proof = checked
            if job['mode'] == 'ic':
                require(report.get('diagnostics', {}).get('implementation_policy') == admitted.get('implementation_policy'),
                        'native/profile candidate policy differs from admission')
                audit = audit_stages(report, fixture, job['algorithm_seed'])
                require(stage_audit is None or audit == stage_audit, 'native/profile stage diagnostics differ')
                stage_audit = audit
            else:
                require(report.get('rho_reusable_setup_excluded') is True, 'rho preparation interval differs')
                require(report.get('executed_method') == {'reference': 'signed_frobenius_rho',
                    'requested_walks': job['config'].get('rho_parallel_walks', 32)},
                    'rho executed a different requested walk width')
        timing = native_intervals(native, process_wall_ns)
        if profile is not None:
            require(costs is not None, 'missing instruction profile')
            if job['mode'] == 'ic':
                ledger = scientific_ledger(profile, costs)
            else:
                require(set(costs) == {'setup', 'reference_solve', 'recovery_check'} and
                        all(type(v) is int and v >= 0 for v in costs.values()) and sum(costs.values()) > 0,
                        'invalid rho instruction phases')
        else:
            require(costs is None, 'unbound instruction profile')
    else:
        require(native is None and profile is None and costs is None,
                'failed raw reports cannot certify a measured record')
    workload = admitted['workload']
    provenance = dict(source_manifest_sha256=admitted['source_manifest_sha256'],
        worker_sha256=admitted['worker_sha256'], host_id=host_id,
        resource_envelope_id=sha256(workload['record']['resource_envelope']),
        calibration_id='valgrind-3.22-amd64-Ir' if profile is not None else 'native-only-no-operation-count',
        report_sha256=report_sha256(profile or native) if status == 'complete' else None)
    if job['mode'] == 'ic':
        record = measured_run(candidate=admitted['candidate'], workload=workload,
            report=profile or native, fixture=fixture, method=admitted['method'], number=number,
            ledger=ledger, native_wall_ns=process_wall_ns, status=status, provenance=provenance,
            diagnostics=stage_audit['queries'] if stage_audit else None, admission_report=admitted['report'])
        record['stage_audit'] = stage_audit
    else:
        reference_id = admitted['reference']['reference_id']
        record = dict(schema_version=2, reference_id=reference_id, workload_id=workload['workload_id'],
            run_id=f"{reference_id}W{workload['workload_id']}R{number}", status=status,
            instruction_phases=costs, total_operations=sum(costs.values()) if costs else None,
            native_wall_ns=process_wall_ns, certificate=proof, provenance=provenance,
            promotion_eligible=False)
    record['native_timing'] = timing
    return record


def online_table(rows, cases, arms, repetitions, rho_aliases):
    """One paired row per target/IC/reference; failures never acquire a speedup."""
    table = []
    for case in cases:
        for arm in arms:
            if arm['id'] in rho_aliases or arm['id'] == 'aa_control':
                continue
            ic = [r for r in rows if r['case'] == case['id'] and r['arm'] == arm['id']]
            for rho_alias in rho_aliases:
                rho = [r for r in rows if r['case'] == case['id'] and r['arm'] == rho_alias]
                complete = all(len(group) == repetitions and
                    {r['repetition'] for r in group} == set(range(repetitions)) and
                    all(r['status'] == 'VERIFIED' and r['measurement']['native_timing'] is not None for r in group)
                    for group in (ic, rho))
                ic_ns = statistics.median(r['measurement']['native_timing']['online']['wall_ns'] for r in ic) if complete else None
                rho_ns = statistics.median(r['measurement']['native_timing']['online']['wall_ns'] for r in rho) if complete else None
                table.append(dict(case=case['id'], public_target=case['fixture']['targets'][0],
                    arm=arm['id'], rho_alias=rho_alias, verified=complete,
                    candidate_ids=sorted({r['measurement']['candidate_id'] for r in ic}),
                    rho_reference_ids=sorted({r['measurement']['reference_id'] for r in rho}),
                    workload_ids=sorted({r['measurement']['workload_id'] for r in ic+rho}),
                    run_ids=[r['measurement']['run_id'] for r in ic+rho],
                    IC_online_ms=ic_ns/1e6 if complete else None,
                    rho_online_ms=rho_ns/1e6 if complete else None,
                    online_speedup=rho_ns/ic_ns if complete else None,
                    boundary='one supplied point, after reusable preparation through scalar replay',
                    aggregation='median of process repetitions of this same point; no target amortization'))
    return table
