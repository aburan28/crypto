"""Admission adapter for the archived optimized producer's actual execution."""
import hashlib
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from identity import base_points, natural
from measurement import PHASES, exclusive_ledger, query_diagnostics
from oracle import Curve, require, verify


def check_build_identity(report, expected, kernel=None):
    source = report.get('diagnostics', report)
    require(source.get('source_manifest_sha256') == expected, 'stale or untagged compiled source')
    if kernel is not None:
        require(source.get('field_kernel') == kernel, 'native/profile/admission field dispatch differs')


def executed_policy(config, panel=None):
    names = {'orbit_batch', 'orbit_target', 'row_kernel', 'full_pair_table'}
    if panel is None:
        require(not (names & config.keys()), 'candidate policy requested from an archived worker')
        return None
    require(panel == 'round1-v1', 'unknown candidate policy adapter')
    policy = dict(orbit_batch=config.get('orbit_batch', 8), orbit_target=config.get('orbit_target'),
                  row_kernel=config.get('row_kernel', 'full'), full_pair_table=config.get('full_pair_table', False))
    require(type(policy['orbit_batch']) is int and policy['orbit_batch'] in (1, 2, 4, 8), 'invalid orbit batch')
    require(policy['orbit_target'] is None or (type(policy['orbit_target']) is int and
            1 <= policy['orbit_target'] <= 64), 'invalid orbit target')
    require(policy['row_kernel'] in ('full', 'suffix', 'bounded', 'word') and
            type(policy['full_pair_table']) is bool, 'invalid row/table policy')
    return policy


def walk_parameters(seed, r):
    """Independent replay of the declared nonzero SplitMix64 scalar sampler."""
    mask = (1 << 64)-1
    state = seed ^ 0x54494E5949430001
    def scalar():
        nonlocal state
        span = r-1
        threshold = (1 << 64) % span
        while True:
            state = (state + 0x9E3779B97F4A7C15) & mask
            z = state
            z = ((z ^ (z >> 30))*0xBF58476D1CE4E5B9) & mask
            z = ((z ^ (z >> 27))*0x94D049BB133111EB) & mask
            word = z ^ (z >> 31)
            if word >= threshold:
                return 1 + word % span
    return scalar(), scalar()


def audit_stages(report, fixture, algorithm_seed):
    """Replay relation scalars, matrix, rank trajectory and all outcome counts.

    Hash the exact matrix in the external certificate audit. Instrumentation
    need not run another cryptographic hash in the measured solver hot path.
    """
    require(report.get('phase_schema') in (2, 3), 'missing scientific producer schema')
    policy = report.get('diagnostics', {}).get('implementation_policy')
    if policy is not None:
        require(policy == executed_policy(policy, 'round1-v1'), 'invalid executed candidate policy')
    full = policy is not None and policy['full_pair_table']
    require(report.get('executed_method') == {
        'pdp': 'full_folded_pair_table' if full else 'partial_folded_pair_table', 'collection': 'additive_walk',
        'relation_la': 'incremental_gauss', 'descent': 'walked_pdp', 'direct_collision': False},
        'unexpected executed backend')
    proof = verify(report, fixture, summands=3)
    c = Curve(fixture)
    require(report['column_convention'] == 'representative', 'wrong column convention')
    base = base_points(report, c)
    mapping = {}
    columns = report['columns']
    if full:
        require(report['pair_table']['rows'] == report['pair_table']['rows_possible'] == columns,
                'full-table policy omitted a row')
    for column, entry in enumerate(report['column_logs']):
        p, coefficient = c.decode(entry['point']), 1
        for _ in range(c.n):
            mapping[p] = column, coefficient
            mapping[c.neg(p)] = column, -coefficient % c.r
            p, coefficient = c.frob(p), coefficient*c.lam % c.r
    initial, stride = walk_parameters(algorithm_seed, c.r)
    trials = natural(report['trials'], 'collection trials')
    first_identity = (-initial*pow(stride, -1, c.r)) % c.r
    identities = 0 if first_identity >= trials else 1+(trials-1-first_identity)//c.r
    d = report['diagnostics']
    for name in ('identity_queries', 'pdp_attempts', 'matrix_nonzeros'):
        natural(d[name], name)
    require(type(d['rank_events']) is list, 'invalid rank events')
    for event in d['rank_events']:
        require(type(event) is list and len(event) == 3, 'invalid rank event')
        for value in event:
            natural(value, 'rank event component')
    require(d['identity_queries'] == identities, 'incorrect identity-query count')
    require(d['ordinary_queries'] == trials-identities, 'incorrect ordinary-query count')
    require(d['pdp_attempts'] == d['ordinary_queries'], 'incorrect PDP attempt count')
    require(d['pdp_outcomes']['verified'] == len(report['relations']), 'incorrect successful-query count')
    require(d['pdp_outcomes']['unresolved'] == d['ordinary_queries']-len(report['relations']),
            'incorrect unresolved-query count')
    checked = query_diagnostics(attempts=d['pdp_attempts'], outcomes=d['pdp_outcomes'],
        ordinary_queries=d['ordinary_queries'], verified_relations=proof['verified_relations'],
        novel_rows=proof['fresh_rows'], final_rank=d['final_rank'], effective_columns=columns)
    seen, trials_seen, pivots, events = set(), set(), {}, []
    matrix = hashlib.sha256()
    nonzeros = 0
    previous_trial = -1
    for rel in report['relations']:
        trial = natural(rel['trial'], 'relation trial')
        require(previous_trial < trial < trials and trial not in trials_seen, 'invalid relation chronology')
        previous_trial = trial
        trials_seen.add(trial)
        require(int(rel['a']) == (initial+trial*stride) % c.r, 'relation was not an ordinary walk query')
        key = int(rel['a']), tuple(sorted(rel['points']))
        if key in seen:
            continue
        seen.add(key)
        row = [0]*columns
        for index in rel['points']:
            column, coefficient = mapping[base[index]]
            row[column] = (row[column]+coefficient) % c.r
        nonzeros += sum(v != 0 for v in row)
        for value in row+[int(rel['a'])]:
            matrix.update(value.to_bytes(8, 'little'))
        before = len(pivots)
        for column in range(columns):
            # Reduction replaces row; read the current coefficient rather than
            # continuing an iterator over the row from before elimination.
            value = row[column]
            if value == 0:
                continue
            if column in pivots:
                row = [(a-value*b) % c.r for a, b in zip(row, pivots[column])]
            else:
                inv = pow(value, -1, c.r)
                pivots[column] = [a*inv % c.r for a in row]
                break
        events.append([trial, before, len(pivots)])
    require(d['rank_events'] == events, 'incorrect rank trajectory')
    require(d['final_rank'] == len(pivots) == proof['rank'], 'incorrect final rank')
    require(d['matrix_nonzeros'] == nonzeros, 'incorrect matrix nonzero count')
    require(d['matrix_encoding'] == 'u64-le-row-major-including-rhs', 'unknown matrix encoding')
    return {'queries': checked, 'rank_events': events, 'matrix_nonzeros': nonzeros,
            'matrix_sha256': matrix.hexdigest(), 'matrix_encoding': d['matrix_encoding'],
            'query_dependence': 'one additive walk; process repetitions are not new query samples'}


def scientific_ledger(report, costs):
    schema = report.get('phase_schema')
    require(schema in (2, 3) and report['mode'] == 'ic', 'not an instrumented IC report')
    target_children = {'target_query', 'target_pdp', 'target_relation_check'} if schema == 3 else set()
    require(set(costs) == (set(PHASES)-{'isogeny'}) | target_children,
            'missing or extra scientific phase intervals')
    cold = {k: v for k, v in costs.items() if k not in target_children}
    cold['target_descent'] += sum(costs[p] for p in target_children)
    return exclusive_ledger(dict(cold, isogeny=0), unit='valgrind-3.22-amd64-Ir',
        process_operations=sum(costs.values()), zero_reasons={'isogeny': 'isogeny:none; no transport executed'})


def method_record(job, fixture, source_manifest, source_manifest_sha256, reference, build, panel=None):
    """Resolve only this audited backend; raw flags cannot name another LA solver."""
    cfg = job['config']
    require(cfg['solver'] == 'pair_table' and cfg['linear_algebra'] == 'tiny_gauss'
            and cfg['summands'] == 3 and cfg.get('collection_window') is None,
            'configuration does not dispatch the instrumented tiny backend')
    require(job['factor_base']['kind'] == 'subgroup_orbits', 'unsupported optimized base')
    require(reference in ('both', 'scaled', 'pairinv'), 'unknown source policy')
    source = source_manifest_sha256
    variant = executed_policy(cfg, panel)
    policy = {'rule': 'retain requested point bound'} if reference == 'both' else {
        'rule': '8*max(1,round(max(1,(r/4196903)^0.157))); retain point bound when rule<=8',
        'arithmetic': 'binary64-powf-round; exact implementation and build pinned'}
    method = {
        'isogeny': 'none',
        'endomorphism': {'order_conductor': None, 'frobenius_order_conductor': None, 'volcano_levels': []},
        'factor_base': {'construction': {'recipe': job['factor_base'], 'orbit_policy': policy,
            'ordering': 'signed-Frobenius-orbit-major', 'subgroup_filter': 'cofactor-project-reject-identity'},
            'nominal_bound': job['factor_base']['points']},
        'point_decomposition': {'summands': 3, 'solver': 'pair',
            'summation_polynomial': 'none; direct three-point decomposition',
            'encoding': 'normal-basis-signed-Frobenius-keyed-partial-pair-table',
            'equation_order': 'cyclic-base-scan-with-frozen-cost-model-blocks',
            'monomial_order': 'none', 'internal_matrix_kernel': 'none',
            'limits': {'scan': 'at-most-one-base-pass'}, 'cache_policy': 'cold-per-job', 'source_sha256': source},
        'relation_collection': {'collector': 'walk',
            'query_distribution': 'nonzero-start-and-stride-SplitMix64-rejection-sampling; correlated-additive-walk',
            'query_rule': 'R[t+1]=R[t]+S; identity queries charged and skipped',
            'filtering': 'first-verified-three-point-witness', 'verification': 'group-readd-every-witness',
            'duplicates': 'scalar-and-sorted-indices', 'dependencies': 'incremental-rank-over-r',
            'stop_rule': {'max_trials': cfg['max_trials'], 'batch_trials': cfg['batch_trials'],
                          'success': 'full-rank-and-certified-column-logs'}, 'source_sha256': source},
        'relation_linear_algebra': {'solver': 'gauss', 'modulus': int(fixture['subgroup_order']),
            'matrix_construction': 'sum-signed-lambda-powers; representative-columns; rhs=a-mod-r',
            'orbit_quotient': 'sign-and-Frobenius', 'rank_criterion': 'full-column-rank',
            'block_parameters': 'none', 'preconditioner': 'none', 'source_sha256': source},
        'target_descent': {'method': 'pdp', 'policy': 'walk-aG+bQ-with-nonzero-b; three-point-witness-required',
            'recursive_solvers': 'same-partial-pair-table; no recursion',
            'success_rule': 'witness-derived-scalar-and-group-replay',
            'stop_rule': {'max_trials': cfg['max_trials']}, 'source_sha256': source},
        'implementation': {'source_manifest_sha256': source, 'components': [
            {'role': role, 'sha256': source_manifest[path]} for role, path in (
                ('worker', 'examples/ic_tournament_worker.rs'),
                ('optimized-ic', 'src/cryptanalysis/koblitz_tiny_ic.rs'),
                ('phase-markers', 'src/cryptanalysis/ic_phase.rs'))],
            'flags': {'build': build, 'reference_policy': reference,
                      'observer_cost': 'fully-charged-markers-native-clocks-and-diagnostics'}},
    }
    if variant is not None:
        construction = method['factor_base']['construction']
        construction['orbit_policy'] = dict(policy, sampler_batch=variant['orbit_batch'],
            explicit_orbit_target=variant['orbit_target'],
            stop='after each declared batch once signed closure reaches the point target')
        if variant['orbit_target'] is not None:
            method['factor_base']['nominal_bound'] = 2 * fixture['degree'] * variant['orbit_target']
        method['implementation']['flags']['candidate_policy'] = variant
        if variant['full_pair_table']:
            method['point_decomposition']['encoding'] = 'normal-basis-signed-Frobenius-keyed-full-pair-table'
            method['target_descent']['recursive_solvers'] = 'same-full-pair-table; no recursion'
    return method
