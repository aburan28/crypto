"""Compose independent generic IC admission and immutable experiment records.

Admitted correctness/provenance is not reference qualification. In particular,
this adapter never produces a speedup, candidate ranking or promotion decision.
"""
import math
import time

from generic_build import verify_binding
from generic_phases import PHASES as RAW_PHASES, verify_native
from generic_stages import effective_config, kernel, verify_stages
from identity import candidate_manifest, natural, run_id, sha256, workload_manifest
from measurement import PHASES, exclusive_ledger, report_sha256
from oracle import require, verify


def scientific_ledger(values, *, unit, process_total):
    """The four target-work phases fold into cold target_descent exactly once."""
    require(set(values) == set(RAW_PHASES) and values['rho_solve'] is None,
            'incorrect IC scientific phase domain')
    costs = {name: values[name] for name in PHASES if name != 'isogeny'}
    children = [values[name] for name in ('target_query', 'target_pdp',
                                         'target_relation_check', 'target_descent')]
    costs['target_descent'] = None if None in children else sum(children)
    costs['isogeny'] = 0
    reasons = {'isogeny': 'isogeny:none; no transport executed'}
    for name, value in costs.items():
        if value == 0 and name != 'isogeny':
            reasons[name] = 'observed interval rounded to zero in the declared unit'
    return exclusive_ledger(costs, unit=unit, process_operations=process_total, zero_reasons=reasons)


def method_record(job, report, stages, build):
    cfg = effective_config(job, report)
    source = build['source_manifest_sha256']
    solver, m = cfg['solver'], cfg['summands']
    algebra = solver not in {'pair_table', 'enumerate'}
    sat = solver.startswith('sat_')
    inherited = solver == 'inherited_f4'
    recipe = stages['base']['recipe']
    base_dimension = stages['base']['construction']['nominal_dimension']
    ambient_sat = sat and m == 3 and base_dimension == job['degree']
    collector = report['collector_dispatch']
    code = dict(pair_table='pair', enumerate='enum', f4='f4', f5='f5',
                inherited_f4='if4', sat_xor='satxor', sat_cnf='satcnf')[solver]
    encoding = ('ambient-S4-circuit-with-finite-base-membership' if ambient_sat else
                'polynomial-basis-Weil-descent-Boolean-quotient' if algebra else
                'full-pair-sum-table' if solver == 'pair_table' else 'group-enumeration')
    limits = {'summands': m}
    if sat:
        limits.update(conflict_budget=cfg['conflict_budget'], max_models=64,
                      macaulay_degree='none' if ambient_sat else 2)
    elif algebra:
        limits.update(max_degree=cfg['groebner_degree'], node_budget=cfg['node_budget'],
                      max_variables=64, split_rule='highest-free' if inherited else 'lowest-free')
    else:
        limits['search'] = 'complete-enumeration' if solver == 'enumerate' else 'source-defined-pair-scan'
    method = dict(
        isogeny='none',
        endomorphism=dict(order_conductor=None, frobenius_order_conductor=None, volcano_levels=[]),
        factor_base=dict(construction=dict(recipe=recipe,
            ordering='constructor-abscissa-order-then-half-trace-root-and-negative',
            subgroup_map='cofactor-multiplication-remove-identity-deduplicate',
            factor_enumeration='degree-at-most-24', sampled_orbit_batch=8),
            nominal_bound=base_dimension),
        point_decomposition=dict(summands=m, solver=code,
            summation_polynomial=('S4-direct' if ambient_sat else 'S3' if m == 2 else 'S3-chain')
            if algebra else 'none', encoding=encoding,
            equation_order=('interleaved-highest-free' if inherited and m >= 3 else
                            'source-defined-ascending-layout') if algebra else 'none',
            monomial_order='degrevlex' if algebra and not sat else 'none',
            internal_matrix_kernel=(('native-xor-CDCL' if solver == 'sat_xor' else 'Tseitin-CNF-CDCL')
                                    + ('' if ambient_sat else '-plus-degree2-Macaulay')) if sat else
                                   'source-default-packed-GF2-dispatch' if algebra else 'none',
            limits=limits, cache_policy='process-local-template-and-layout-caches-source-defaults',
            source_sha256=source),
        relation_collection=dict(collector='walk' if collector['collection_window'] is not None else 'sample',
            query_distribution='pinned-StdRng08-nonzero-u64-samples; walked-queries-correlated',
            query_rule=dict(policy=collector['query_rule'],
                            window=collector['collection_window'] or 'none'),
            filtering='first-group-verified-witness', verification='group-readd-every-relation',
            duplicates='scalar-and-sorted-base-indices', dependencies='retain-unique-even-if-dependent',
            stop_rule=dict(batch_trials=cfg['batch_trials'], max_trials=cfg['max_trials'],
                           success='certified-column-logs; independent-admission-requires-full-rank'),
            source_sha256=source),
        relation_linear_algebra=dict(solver='bw' if cfg['linear_algebra'] == 'sparse' else 'gauss',
            modulus=int(report['fixture']['subgroup_order']),
            matrix_construction='cofactor-project; sorted-sign-Frobenius-columns; rhs=h*a; sum-signed-lambda-powers',
            orbit_quotient='sign-and-Frobenius', rank_criterion='independent-full-column-rank',
            block_parameters=cfg['sparse'] if cfg['linear_algebra'] == 'sparse' else 'none',
            preconditioner='filter-fold-reconstruct' if cfg['linear_algebra'] == 'sparse' else 'none',
            source_sha256=source),
        target_descent=dict(method='pdp',
            policy='walked-aG+bQ-64' if solver == 'pair_table' else 'sampled-aG+bQ',
            recursive_solvers='same-PDP; no-recursion; direct-collision-disabled',
            success_rule='relation-derived-scalar-and-general-group-replay',
            stop_rule=dict(max_trials=cfg['max_trials']), source_sha256=source),
        implementation=dict(source_manifest_sha256=source,
            components=[dict(role='complete-worker-library-and-dependency-manifest', sha256=source)],
            flags=dict(build_sha256=build['build_sha256'],
                actual_field_kernel=collector['field_kernel'],
                runtime_policy='default-environment-one-rayon-v1',
                observer_policy='exclusive-owner-thread-v1; queries-matrix-and-dispatch-charged')))
    return method


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


def admit_rho(report, fixture, job, build, source, *, executable, process_wall_ns):
    require(job['mode'] == 'rho' and type(report.get('generic_admission_schema')) is int
            and report['generic_admission_schema'] == 1
            and report.get('generic_runtime_policy') == 'default-environment-one-rayon-v1',
            'missing rho admission policy')
    binding = verify_binding(report, build, source, executable=executable)
    cfg = effective_config(job, report)
    observed = report['rho_dispatch']
    field = kernel(report['field_kernel'])
    require(sha256(observed) == sha256(dict(algorithm='signed-frobenius-batched-affine', jump_count=16,
        max_restarts=64, max_iterations_per_restart=cfg['max_trials'], progress_interval=256,
        requested_walks=cfg['rho_parallel_walks'], seed=job['algorithm_seed'], field_kernel=field)),
        'rho dispatch mismatch')
    expected_walks = max(1, min(cfg['rho_parallel_walks'],
        int(math.sqrt(math.pi*int(fixture['subgroup_order'])/2)/math.sqrt(2*fixture['degree'])/64)))
    require(len(report['solutions']) == 1 and type(report['solutions'][0]['effective_walks']) is int
            and report['solutions'][0]['effective_walks'] == expected_walks,
            'rho effective width mismatch')
    clocks = verify_native(report, job, process_wall_ns=process_wall_ns)
    certificate = verify(report, fixture, expected_mode='rho') if report['status'] == 'complete' else None
    return dict(schema_version=1, status='PASS', binding=binding, phases=clocks, certificate=certificate,
                dispatch=observed, performance_qualified=False, promotion_eligible=False)
