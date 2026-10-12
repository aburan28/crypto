"""Fail-closed scientific accounting, separate from historical profiler labels."""
import copy
import hashlib
import json

from identity import (candidate_manifest, digest, natural, run_id, sha256,
                      workload_manifest)
from oracle import require, verify


PHASES = ('setup', 'isogeny', 'factor_base', 'precompute', 'queries', 'pdp',
          'relation_check', 'matrix_build', 'relation_la', 'target_descent', 'recovery_check')
PDP_OUTCOMES = ('verified', 'proved_unsat', 'timeout', 'budget', 'error', 'lift_rejected',
                'unresolved')


def report_sha256(report):
    # Worker reports contain diagnostic elapsed seconds. Floats are permitted in
    # run artifacts, but never in canonical mathematical/method identities.
    return hashlib.sha256(json.dumps(report, sort_keys=True, separators=(',', ':'),
        ensure_ascii=False, allow_nan=False).encode('utf-8')).hexdigest()


def exclusive_ledger(operations, *, unit, process_operations, zero_reasons):
    require(type(operations) is dict and set(operations) == set(PHASES),
            'ledger must explicitly contain all eleven phases')
    require(type(unit) is str and unit, 'missing operation unit')
    require(type(zero_reasons) is dict and set(zero_reasons) <= set(PHASES),
            'unknown zero-cost phase')
    if process_operations is not None:
        natural(process_operations, 'whole-process operations', positive=True)
    known = 0
    missing = []
    for phase, value in operations.items():
        if value is None:
            missing.append(phase)
        else:
            natural(value, phase + ' operations')
            known += value
            if value == 0:
                require(type(zero_reasons.get(phase)) is str and zero_reasons[phase],
                        'zero cost needs an explicit absence/accounting reason: ' + phase)
        require(phase not in zero_reasons or value == 0, 'zero reason for a nonzero or unknown phase')
    if process_operations is not None:
        require(known <= process_operations, 'phase costs exceed whole process')
    closed = not missing and process_operations is not None
    if closed:
        require(known == process_operations, 'exclusive phases do not close against the whole process')
    return {'unit': unit, 'operations': copy.deepcopy(operations),
            'zero_reasons': copy.deepcopy(zero_reasons), 'missing_phases': sorted(missing),
            'whole_process_operations': process_operations,
            'unattributed_operations': None if process_operations is None else process_operations-known,
            'complete': closed, 'cold_operations': known if closed else None}


def legacy_ledger(costs, *, unit, process_operations):
    """Combined historical labels cannot be reverse-engineered into stage costs.

    Preserve their independently closed process measurement as a diagnostic.
    Even wholly named intervals exclude some setup/reporting work, so do not
    silently claim those are complete scientific phase measurements either.
    """
    require(type(costs) is dict and costs, 'missing historical intervals')
    for label, cost in costs.items():
        require(type(label) is str and label, 'invalid historical interval label')
        natural(cost, 'historical interval cost')
    require(sum(costs.values()) == process_operations, 'historical intervals do not close')
    ledger = exclusive_ledger(dict.fromkeys(PHASES), unit=unit,
                              process_operations=process_operations, zero_reasons={})
    ledger['legacy_intervals'] = copy.deepcopy(costs)
    ledger['limitation'] = 'Combined historical intervals; eleven scientific phases were not measured.'
    return ledger


def query_diagnostics(*, attempts, outcomes, ordinary_queries, verified_relations,
                      novel_rows, final_rank, effective_columns):
    """No-result from a bounded solver is unresolved, never a proof of UNSAT."""
    natural(attempts, 'PDP attempts')
    require(type(outcomes) is dict and set(outcomes) == set(PDP_OUTCOMES),
            'missing PDP outcome categories')
    for status, count in outcomes.items():
        natural(count, 'PDP ' + status)
    require(sum(outcomes.values()) == attempts, 'PDP outcomes do not sum to attempts')
    for name, count in [('ordinary queries', ordinary_queries), ('verified relations', verified_relations),
                        ('novel rows', novel_rows), ('final rank', final_rank),
                        ('effective columns', effective_columns)]:
        natural(count, name)
    require(final_rank <= min(novel_rows, effective_columns), 'impossible final rank')
    require(novel_rows <= verified_relations, 'novel rows exceed verified relations')
    # A query may yield multiple decompositions; verified outcome counts queries,
    # while verified_relations counts rows. Do not conflate their denominators.
    require(outcomes['verified'] <= ordinary_queries <= attempts, 'inconsistent ordinary-query counts')
    return dict(attempts=attempts, outcomes=copy.deepcopy(outcomes), ordinary_queries=ordinary_queries,
                verified_relations=verified_relations, novel_rows=novel_rows, final_rank=final_rank,
                effective_columns=effective_columns)


def measured_run(*, candidate, workload, report, fixture, method, number, ledger,
                 native_wall_ns, status, provenance, diagnostics=None, admission_report=None):
    """Re-derive identities and certificates; supplied labels cannot bless a run.

    An incomplete ledger may have a verified answer, but has no complete cost
    or speedup. Failed executions still retain the frozen candidate/workload key.
    Their base was checked during admission, before this measured execution.
    """
    require(status in ('complete', 'timeout', 'budget', 'oom', 'error', 'insufficient_relations'),
            'unknown run status')
    record = workload['record']
    expected_workload = workload_manifest(fixture, input_law=record['input_law'],
        algorithm_seed=record['algorithm_seed'], resource_envelope=record['resource_envelope'],
        cache_policy=record['cache_policy'])
    require(workload == expected_workload, 'changed workload manifest')
    admission = admission_report if admission_report is not None else report
    require(admission is not None, 'failed run requires the frozen base admission')
    require(candidate_manifest(fixture, admission, method) == candidate, 'changed candidate admission')
    require(candidate['record_sha256'] == sha256(candidate['record']), 'changed candidate record')
    require(candidate['record']['curve']['curve_id'] == record['curve_id'], 'candidate/workload curve mismatch')
    require(candidate['record']['point_decomposition'] == method['point_decomposition'], 'changed PDP method')
    # Re-derive the full candidate even for failed jobs from its admitted base.
    # Failed jobs may have no report, so retain the immutable admission separately.
    for section in method:
        actual = candidate['record'][section]
        if section == 'factor_base':
            actual = {k: v for k, v in actual.items() if k != 'inventory'}
        require(actual == method[section], 'changed method section ' + section)
    checked = exclusive_ledger(ledger['operations'], unit=ledger['unit'],
        process_operations=ledger['whole_process_operations'], zero_reasons=ledger['zero_reasons'])
    require(all(ledger[k] == v for k, v in checked.items()), 'changed phase ledger')
    if native_wall_ns is not None:
        natural(native_wall_ns, 'native cold wall nanoseconds', positive=True)
    required_provenance = {'source_manifest_sha256', 'worker_sha256', 'host_id',
                           'resource_envelope_id', 'calibration_id', 'report_sha256'}
    require(type(provenance) is dict and set(provenance) == required_provenance, 'missing run provenance')
    for key in ('source_manifest_sha256', 'worker_sha256'):
        digest(provenance[key], key)
    for key in ('host_id', 'resource_envelope_id', 'calibration_id'):
        require(type(provenance[key]) is str and provenance[key], 'missing ' + key)
    require(provenance['source_manifest_sha256'] == candidate['record']['implementation']['source_manifest_sha256'],
            'run source differs from candidate')
    require(provenance['resource_envelope_id'] == sha256(record['resource_envelope']),
            'changed run resources')
    require(provenance['report_sha256'] == (report_sha256(report) if report is not None else None),
            'changed run report')
    certificate = None
    if status == 'complete':
        require(report is not None, 'missing complete report')
        require(candidate_manifest(fixture, report, method) == candidate, 'changed candidate/base identity')
        certificate = verify(report, fixture, summands=method['point_decomposition']['summands'])
    elif report is not None and ('factor_base' in report or 'factor_base_orbits' in report):
        require(candidate_manifest(fixture, report, method) == candidate, 'changed failed-run base')
    if diagnostics is not None:
        require(query_diagnostics(**diagnostics) == diagnostics, 'changed query diagnostics')
        if certificate is not None:
            require(diagnostics['verified_relations'] == certificate['verified_relations'] and
                    diagnostics['novel_rows'] == certificate['fresh_rows'] and
                    diagnostics['final_rank'] == certificate['rank'], 'diagnostics disagree with certificate')
    complete = certificate is not None and ledger['complete'] and record['cache_policy'] == 'cold'
    total = ledger['cold_operations'] if complete else None
    return {'schema_version': 2, 'candidate_id': candidate['candidate_id'],
            'workload_id': workload['workload_id'], 'run_id': run_id(candidate['candidate_id'], workload['workload_id'], number),
            'status': status, 'unit': ledger['unit'], 'phase_ledger': copy.deepcopy(ledger),
            'total_operations': total, 'native_wall_ns': native_wall_ns,
            'complete_cold_cost': complete, 'promotion_eligible': False,
            'promotion_note': 'A run record does not establish paired confirmation or a promotion decision.',
            'certificate': certificate, 'diagnostics': diagnostics, 'provenance': copy.deepcopy(provenance)}
