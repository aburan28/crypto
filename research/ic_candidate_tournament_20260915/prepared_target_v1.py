"""Target-only IC adapter cores for the accepted public-synthetic n17 state.

These APIs do not sample targets, dispatch a campaign or certify a speedup.
prepared_f5_runtime_v1 registers the native job against a new build. SAT's
query callback must be the isolated source-bound query runner in its registration.
Historical complete-solve invocations remain consumed and unchanged.
"""
from collections import Counter
import copy
import time

from generic_query_law import descent_coefficients
from generic_queries_exact_v1 import verify_queries
from generic_stages import effective_config, kernel
from identity import curve_record, natural, sha256
from oracle import Curve, require
from prepared_ic_state_v1 import STATE_SHA256, verify
from static_sat_matrix import RelationMatrix

CERTIFICATE_SEALS = {
    'f5': 'd8c5d6679fe89561606154785fdfba763835bf261dceb08fae193eff9f5636ae',
    'sat': '91856ab78550436d3f668367f9aebd9e2c0604bd64b1472d9d19ec318e2b144e',
}
ONLINE_PHASES = ('target_query', 'target_pdp', 'target_relation_check',
                 'target_descent', 'target_recovery_check')


def load_preparation(document, expected_sha256, *, family):
    require(family in CERTIFICATE_SEALS and expected_sha256 == CERTIFICATE_SEALS[family],
            'target adapter requires the exact external family certificate seal')
    receipt = verify(document, expected_sha256)
    require(document['provenance']['family'] == family,
            'target adapter preparation family differs from its run evidence')
    inputs = document['certificate']['inputs']
    curve = Curve(inputs['fixture'])
    base = tuple(curve.decode(p) for p in inputs['base'])
    matrix = RelationMatrix(curve, base)
    for attempt in inputs['attempts']:
        if attempt['outcome'] in ('witness', 'VALID_POINT_WITNESS'):
            matrix.push(attempt['scalar'], attempt['indices'])
    logs = tuple(matrix.solve())
    require(matrix.snapshot() == document['certificate']['proof']['matrix']
            and logs == tuple(item['log'] for item in document['record']['column_logs']),
            'loaded preparation differs from independently replayed matrix or logs')
    return curve, base, matrix, logs, receipt


def supplied_target(curve, encoded):
    require(type(encoded) is list and len(encoded) == 2
            and all(type(value) is int and 0 <= value < 2**curve.n for value in encoded),
            'target adapter requires canonical integer public coordinates')
    target = curve.decode(encoded)
    require(target is not None and curve.mul(target, curve.r) is None,
            'target adapter public input is outside the nonidentity subgroup')
    return target


def native_job(document, expected_sha256, *, point, algorithm_seed, max_attempts=8):
    """Prepare stdin for the new native mode, never execute the old invocation.

    Only mathematical state enters native reusable inputs. Certificate source
    history, original target/scalar, run seeds and costs are not copied into it.
    A later registrar binds the complete certificate in the run manifest.
    """
    curve, base, matrix, logs, _ = load_preparation(document, expected_sha256, family='f5')
    target = supplied_target(curve, point)
    require(type(algorithm_seed) is int and 0 <= algorithm_seed < 2**64,
            'native target seed outside pinned RNG domain')
    require(type(max_attempts) is int and 1 <= max_attempts <= 8,
            'native target control cap outside frozen eight-attempt envelope')
    return dict(mode='ic', degree=17, curve_a=1, public_targets=[[str(v) for v in target]],
        target_seeds=[], algorithm_seed=algorithm_seed, exclusive_phases=True,
        factor_base=dict(kind='standard_subspace', dimension=6),
        config=dict(solver='f5', linear_algebra='dense', summands=3,
                    groebner_degree=3, node_budget=8192, conflict_budget=100000,
                    batch_trials=1, max_trials=max_attempts),
        prepared=dict(mathematical_state_sha256=STATE_SHA256,
            factor_base=[[str(v) for v in p] for p in base],
            columns=[dict(point=[str(v) for v in p], log=str(log))
                     for p, log in zip(matrix.columns, logs)]))


def audit_native_target(report, job, document, expected_sha256):
    """Independently audit a new warm report; source/build transport is separate."""
    curve, base, matrix, logs, preparation = load_preparation(
        document, expected_sha256, family='f5')
    require(job == native_job(document, expected_sha256,
                point=[int(v) for v in job['public_targets'][0]],
                algorithm_seed=job['algorithm_seed'], max_attempts=job['config']['max_trials']),
            'native target job differs from the exact prepared control method')
    target = supplied_target(curve, [int(v) for v in job['public_targets'][0]])
    require(report['status'] in ('complete', 'incomplete')
            and report['preparation_mode'] == 'imported-certified-log-table-v1'
            and report['preparation_mathematical_state_sha256'] == STATE_SHA256
            and report['reusable_symbolic_template_prepared'] is True
            and report['factor_base'] == job['prepared']['factor_base']
            and report['column_logs'] == job['prepared']['columns']
            and report['trials'] == 0 and report['relations'] == []
            and report['collection_reports'] == [] and report['solve_attempts'] == 0,
            'target-only native report ran ordinary queries or changed preparation')
    # This exact query-law adapter proves all recorded negative m=3 queries and
    # group-readds every witness, including repeated geometry and killed torsion.
    require(curve_record(report['fixture']) == document['record']['curve']
            and report['fixture']['targets'] == job['public_targets']
            and report['fixture']['target_scalar_constructed'] is False,
            'native target curve or supplied point differs from preparation')
    query_audit = verify_queries(report, report['fixture'], 3)
    effective_config(job, report)
    dispatch = report['descent_dispatch']
    require(dispatch == dict(strategy='Groebner', field_kernel=kernel(dispatch['field_kernel']),
                pair_table=False, query_rule='seeded-sample', summands=3, direct_collision=False),
            'prepared native target dispatch substituted a different mechanism')
    solutions = report['solutions']
    require(len(solutions) == 1, 'target-only native report has multiple targets')
    solution = solutions[0]
    stream = descent_coefficients(job['algorithm_seed'], curve.r, walked=False)
    require(type(solution['trials']) is int and 0 < solution['trials'] <= job['config']['max_trials'],
            'native target attempt count exceeds its frozen cap')
    for attempt in solution['attempts']:
        require((attempt['a'], attempt['b']) == next(stream), 'native target query stream changed')
        stats = attempt['pdp']['stats']
        if attempt['pdp']['outcome'] != 'identity':
            require(set(stats) == {'family', 'engine', 'stats'} and stats['family'] == 'groebner'
                    and stats['engine'] == {'MatrixF5':{'max_degree':3}},
                    'prepared native target did not execute the declared MatrixF5 engine')
    complete = report['status'] == 'complete'
    if complete:
        relation = solution['relation']
        require(type(relation) is dict and relation.get('points') is not None,
                'target-only answer lacks a factor-base decomposition')
        scalar = matrix.recover_target(logs, target, relation['a'], relation['b'], relation['points'])
        require(str(scalar) == solution['recovered'] and report['scalar_verified'] is True,
                'target-only native scalar differs from independent relation recovery')
    else:
        require(solution['recovered'] is None and report['scalar_verified'] is False,
                'incomplete native target claims a scalar')
        require(solution['trials'] == job['config']['max_trials'], 'incomplete target stopped before cap')
        scalar = None
    trace = report['generic_phase_timing']
    require(report['generic_phase_policy'] == 'exclusive-owner-thread-v1'
            and report['online_timing_schema'] == 2
            and report['reusable_setup_excluded'] is True,
            'target-only native interval policy differs')
    observed = trace['online_phases_ns']
    require(set(observed) == {'target_query', 'target_pdp', 'target_relation_check',
                             'target_descent', 'recovery_check', 'rho_solve'}
            and observed['rho_solve'] is None,
            'target-only native schema gained rho work or omitted phase slots')
    raw = {k:v for k,v in observed.items() if k != 'rho_solve'}
    require(all((v is None and not complete) or (type(v) is int and v >= 0) for v in raw.values()),
            'target-only native exclusive phases missing or malformed')
    attempted = natural(trace['online_wall_ns'], 'native attempted online ns', positive=True)
    require(sum(v for v in raw.values() if v is not None) == attempted == report['online_wall_ns']
            and report['outer_online_wall_ns'] <= attempted <= trace['observed_wall_ns'],
            'target-only native phase clocks fail closure')
    phases = {('target_recovery_check' if k == 'recovery_check' else k):v for k,v in raw.items()}
    return dict(status='COMPLETE' if complete else 'INCOMPLETE_TARGET',
        preparation=preparation, ordinary_queries_executed=0, query_audit=query_audit,
        recovered_scalar=scalar, scalar_verified=complete,
        online_attempt_wall_ns=attempted, online_wall_ns=attempted if complete else None,
        online_phases_ns=phases, source_bound_execution_admitted=False,
        scope='independent mathematical and online-clock audit; source transport not supplied',
        fresh_paired_qualification=False, promotion_eligible=False, online_speedup=None)


def solve_sat_target(document, expected_sha256, *, point, panel, query,
                     progress=None, clock=time.monotonic_ns):
    """Core used by a future sealed SAT entrypoint; callback runs one bound query.

    Fake callbacks in tests are mathematical controls, not solver yield. Native
    use must pass static_sat_query_v3.one_query through the frozen meter context.
    """
    started = clock()
    curve, base, matrix, logs, preparation = load_preparation(document, expected_sha256, family='sat')
    target = supplied_target(curve, point)
    require(type(panel['descent_query_seed']) is int
            and 0 <= panel['descent_query_seed'] < 2**64
            and type(panel['max_descent_queries']) is int
            and 1 <= panel['max_descent_queries'] <= 8,
            'SAT target control stream or cap outside frozen envelope')
    require(callable(query) and (progress is None or callable(progress)), 'invalid SAT query adapter')
    stream = descent_coefficients(panel['descent_query_seed'], curve.r, walked=False)
    phases = {name:0 for name in ONLINE_PHASES}
    attempts, recovered = [], None
    preparation_ns = clock()-started
    online_start = clock()
    for trial in range(panel['max_descent_queries']):
        start = clock()
        a,b = next(stream)
        public = curve.add(curve.mul(curve.g, a), curve.mul(target, b))
        phases['target_query'] += clock()-start
        item = dict(trial=trial, probe_scalar=a, point=None if public is None else list(public))
        start = clock()
        row = copy.deepcopy(query(panel, item, curve, base))
        query_ns = clock()-start
        verification_ns = natural(row['verification_wall_ns'], 'SAT query verification ns')
        require(verification_ns <= query_ns, 'SAT verification exceeds exclusive query interval')
        require(row['trial'] == trial and row['probe_scalar'] == a
                and row['public_point'] == item['point'], 'SAT query receipt input differs')
        phases['target_pdp'] += query_ns-verification_ns
        phases['target_relation_check'] += verification_ns
        row.update(a=a, b=b, target_query_index=trial, query_wall_ns=query_ns,
                   pdp_wall_ns=query_ns-verification_ns)
        if row['status'] == 'VALID_POINT_WITNESS':
            require(row['source_model_valid'] is True and type(row['point_witness']) is dict
                    and row['point_witness']['group_replay'] is True,
                    'SAT witness lacks source-model and lifting verification')
            indices = row['point_witness']['point_indices']
            require(type(indices) is list and len(indices) == 3
                    and all(type(i) is int and 0 <= i < len(base) for i in indices),
                    'SAT target witness indices malformed')
            start = clock()
            total = None
            for index in indices:
                total = curve.add(total, base[index])
            require(total == public, 'SAT target witness does not group-readd to aG+bQ')
            phases['target_relation_check'] += clock()-start
            start = clock()
            candidate = matrix.descent_scalar(logs, a,b,indices)
            phases['target_descent'] += clock()-start
            start = clock()
            require(curve.mul(curve.g,candidate) == target, 'SAT target scalar replay failed')
            phases['target_recovery_check'] += clock()-start
            recovered = candidate
            online_stop = clock()
            row.update(candidate_scalar=str(candidate), scalar_replay_verified=True)
        else:
            require(row['status'] in ('IDENTITY_QUERY', 'EXPORT_FAILURE', 'INVALID_EXPORT',
                    'TIMEOUT', 'INVALID_SOURCE_MODEL', 'SOURCE_MODEL_NONLIFTING',
                    'SOURCE_UNSAT', 'CONFLICT_BUDGET_INCONCLUSIVE', 'UNKNOWN_INCONCLUSIVE',
                    'SOLVER_ERROR'), 'unknown SAT attempt status')
        attempts.append(row)
        if recovered is not None:
            break
        if progress is not None:
            progress(row)  # All intermediate reporting is charged to this target.
    if recovered is None:
        online_stop = clock()
    online_ns = online_stop-online_start
    residual = online_ns-sum(phases.values())
    require(residual >= 0 and online_ns > 0, 'SAT exclusive online clocks overlap or vanish')
    phases['target_query'] += residual
    return dict(status='COMPLETE' if recovered is not None else 'INCOMPLETE_TARGET',
        preparation=preparation, preparation_wall_ns=preparation_ns,
        ordinary_queries_executed=0, target_input=list(target), target_attempts=attempts,
        target_status_mix=dict(sorted(Counter(row['status'] for row in attempts).items())),
        recovered_scalar=recovered, scalar_verified=recovered is not None,
        online_wall_ns=online_ns if recovered is not None else None,
        online_attempt_wall_ns=online_ns, online_phases_ns=phases,
        online_start_monotonic_ns=online_start, online_stop_monotonic_ns=online_stop,
        online_stop_event='independent-scalar-replay' if recovered is not None else 'frozen-target-attempt-cap',
        online_bookkeeping_assigned_to_target_query_ns=residual,
        source_bound_execution_admitted=False, fresh_paired_qualification=False,
        scope='target adapter core; source-bound callback and transport require separate admission',
        promotion_eligible=False, online_speedup=None)


def preparation_exposures(documents):
    """Required preparation exclusion union; not the complete historical census."""
    require(set(documents) == set(CERTIFICATE_SEALS), 'preparation exposure union omits an arm')
    points, categories = set(), {}
    for family, document in documents.items():
        curve, _, matrix, _, _ = load_preparation(document, CERTIFICATE_SEALS[family], family=family)
        ordinary = {curve.mul(curve.g, row['scalar']) for row in document['certificate']['inputs']['attempts']}
        known = set()
        for column in matrix.columns:
            current = column
            for _ in range(curve.n):
                known.update((current, curve.neg(current)))
                current = curve.frob(current)
        require(None not in ordinary | known, 'identity leaked into preparation exposure set')
        points.update(ordinary | known)
        categories[family] = dict(ordinary_query_points=[list(p) for p in sorted(ordinary)],
            known_log_orbit_points=[list(p) for p in sorted(known)])
    record = dict(schema_version=1, curve_id=documents['f5']['record']['curve']['curve']['curve_id'],
        source_certificate_seals=CERTIFICATE_SEALS, categories=categories,
        points=[list(p) for p in sorted(points)],
        scope='preparation exclusions only; all historical and current exposures still required')
    return dict(record=record, record_sha256=sha256(record), point_count=len(points),
                fresh_sampling_authorized=False)
