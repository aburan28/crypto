"""Versioned complete IC admission with bounded exact group negatives.

Source/build, observed dispatch, base/matrix/log/scalar and exclusive accounting
checks reuse the accepted implementations. No native execution occurs here.
"""
import time
from generic_admission import scientific_ledger, method_record
from generic_build import verify_binding
from generic_phases import verify_native
from generic_stages_exact_v1 import verify_stages
from identity import candidate_manifest, natural, run_id, workload_manifest
from measurement import report_sha256
from oracle import require


def admit(report, fixture, job, build, source, *, executable, process_wall_ns,
          resources, number):
    natural(number, 'run number')
    started = time.monotonic_ns()
    require(report.get('generic_runtime_policy') == 'default-environment-one-rayon-v1',
            'missing runtime override policy')
    binding = verify_binding(report, build, source, executable=executable)
    require(job['mode'] == 'ic', 'IC admission requires an IC job')
    stages = verify_stages(report, fixture, job)
    clocks = verify_native(report, job, process_wall_ns=process_wall_ns)
    method = method_record(job, report, stages, build)
    candidate = candidate_manifest(fixture, report, method)
    workload = workload_manifest(fixture, input_law='one-supplied-public-point; seed-is-provenance',
        algorithm_seed=job['algorithm_seed'], resource_envelope=resources, cache_policy='cold')
    ledger = scientific_ledger(clocks['process_phases_ns'], unit='native_wall_ns',
                               process_total=process_wall_ns)
    complete = report['status'] == 'complete'
    require(not complete or ledger['complete'], 'complete candidate has missing scientific costs')
    run = dict(schema_version=1, candidate_id=candidate['candidate_id'], workload_id=workload['workload_id'],
        run_id=run_id(candidate['candidate_id'], workload['workload_id'], number), status=report['status'],
        headline_metric='one-target-online-native-wall-ns', online_wall_ns=clocks['online_wall_ns'],
        online_phases_ns=clocks['online_phases_ns'], cold_phase_ledger=ledger,
        cold_wall_ns=process_wall_ns if complete else None,
        certificate=stages['certificate'], report_sha256=report_sha256(report), binding=binding,
        independent_audit_wall_ns=time.monotonic_ns()-started,
        independent_audit_timing='external Python audit excluded from worker online interval; worker scalar replay included',
        qualification=None, performance_qualified=False, promotion_eligible=False, online_speedup=None,
        instruction_cost=None, normalized_S=None,
        limitation='scientific admission only; scheduling/observer and matched reference qualification pending')
    return dict(schema_version=1, status='PASS', candidate=candidate, workload=workload,
                run=run, stages=stages, phases=clocks, promotion_eligible=False)
