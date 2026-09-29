"""Independent bounded dispatch, stored-matrix and batch-stop admission.

Successful group equations alone do not prove that a named pipeline ran. This
adapter binds them to observed dispatch and the actual LA input, including every
unsuccessful batch. It does not establish performance qualification.
"""
import copy

from generic_bases import verify_base
from generic_query_law import SOLVERS, verify_query_law
from identity import sha256
from oracle import Curve, require, verify

SPARSE = dict(filter=dict(remove_duplicates=True, remove_singletons=True,
                         target_excess=32, merge_max_weight=8, max_row_weight=32),
              wiedemann=dict(block_m=4, block_n=4, margin=8),
              fold=3, attempts=3, seed=0x5350415253454c41)
DEFAULTS = dict(solver='pair_table', linear_algebra='sparse', batch_trials=64,
                max_trials=4096, summands=3, collection_window=None,
                sparse=SPARSE, factor_base_orbits=None, factor_base_cube_root=False,
                factor_base=None, groebner_degree=3, node_budget=4096,
                conflict_budget=100_000, rho_parallel_walks=32)
STRATEGIES = dict(pair_table='PairTable', enumerate='Enumerate', f4='Groebner',
                  f5='Groebner', inherited_f4='Groebner', sat_xor='Sat', sat_cnf='Sat')


def effective_config(job, report):
    requested = job['config']
    require(type(requested) is dict and set(requested) <= set(DEFAULTS), 'unknown configuration')
    cfg = copy.deepcopy(DEFAULTS)
    cfg.update(copy.deepcopy(requested))
    def merge_sparse(defaults, supplied):
        require(type(supplied) is dict and set(supplied) <= set(defaults), 'unknown sparse option')
        return {key: (merge_sparse(value, supplied.get(key, {})) if type(value) is dict
                      else supplied.get(key, value)) for key, value in defaults.items()}
    cfg['sparse'] = merge_sparse(SPARSE, requested.get('sparse', {}))
    require(sha256(cfg) == sha256(report.get('effective_config')), 'effective configuration differs from job')
    require(cfg['solver'] in SOLVERS and cfg['linear_algebra'] in {'dense', 'sparse'},
            'unsupported observed configuration')
    return cfg


def kernel(value):
    require(value in {'portable', 'pclmulqdq', 'pmull'}, 'missing/unknown actual field kernel')
    return value


def dispatch(report, cfg):
    strategy, pair = STRATEGIES[cfg['solver']], cfg['solver'] == 'pair_table'
    window = cfg['collection_window']
    active = (pair and cfg['summands'] == 3 and window is not None
              and 0 < window < len(report['factor_base']))
    collector = report['collector_dispatch']
    field = kernel(collector['field_kernel'])
    require(sha256(collector) == sha256(dict(strategy=strategy, field_kernel=field, pair_table=pair,
                              query_rule='windowed-walk-64' if active else 'trial-keyed-sample',
                              collection_window=window if active else None)),
            'collector dispatch mismatch')
    solutions = report.get('solutions', [])
    descent = report['descent_dispatch']
    if solutions:
        require(sha256(descent) == sha256(dict(strategy=strategy, field_kernel=field, pair_table=pair,
                                query_rule='parallel-walk-64' if pair else 'seeded-sample',
                                summands=cfg['summands'], direct_collision=False)),
                'descent dispatch mismatch')
    else:
        require(descent is None, 'unstarted descent has dispatch evidence')
    attempts = [a for b in report['collection_reports'] for a in b['attempts']]
    attempts += [a for s in solutions for a in s['attempts']]
    families = {}
    for attempt in attempts:
        pdp, solver = attempt['pdp'], cfg['solver']
        stats = pdp['stats']
        if pdp['outcome'] == 'identity' or solver in {'pair_table', 'enumerate'}:
            require(stats == {'family': 'none'}, 'unexpected algebra dispatch')
        elif solver in {'f4', 'f5', 'inherited_f4'}:
            engine = dict(f4='MatrixF4', f5='MatrixF5', inherited_f4='InheritedF4')[solver]
            require(set(stats) == {'family', 'engine', 'stats'} and stats['family'] == 'groebner'
                    and stats['engine'] == {engine: {'max_degree': cfg['groebner_degree']}},
                    'Groebner engine substitution')
        else:
            require(set(stats) == {'family', 'backend', 'stats'} and stats['family'] == 'sat'
                    and stats['backend'] == ('cnf' if solver == 'sat_cnf' else 'native_xor'),
                    'SAT backend substitution')
        families[stats['family']] = families.get(stats['family'], 0) + 1
    return dict(collector=collector, descent=descent, frontend_calls=families)


def projected_columns(curve, base):
    projected = [curve.mul(p, curve.h) for p in base]

    def orbit(p):
        points = []
        for _ in range(curve.n):
            points.extend((p, curve.neg(p)))
            p = curve.frob(p)
        return points

    representatives = sorted({min(orbit(p)) for p in projected if p is not None})
    mapping = {}
    for j, p in enumerate(representatives):
        coefficient = 1
        for _ in range(curve.n):
            for q, value in ((p, coefficient), (curve.neg(p), -coefficient % curve.r)):
                require(q not in mapping or mapping[q] == (j, value), 'ambiguous orbit coefficient')
                mapping[q] = j, value
            p, coefficient = curve.frob(p), coefficient * curve.lam % curve.r
    return representatives, [mapping.get(p) for p in projected]


def matrix_audit(report, fixture, cfg, *, inventory=False):
    c = Curve(fixture)
    base = [c.decode(p) for p in report['factor_base']]
    representatives, projected = projected_columns(c, base)
    columns = [[str(x), str(y)] for x, y in representatives]
    require(len(columns) == report['columns'] > 0, 'matrix column count mismatch')
    snapshot = report['relation_matrix']
    sparse = cfg['linear_algebra'] == 'sparse'
    # Every admitted subgroup is below 2^31, inside the sparse u64 modulus gate.
    require(snapshot['modulus'] == str(c.r) and snapshot['column_points'] == columns,
            'matrix modulus or column order mismatch')
    require(snapshot['solver'] == ('sparse-filter-block-wiedemann' if sparse else 'dense-gauss')
            and sha256(snapshot['sparse_options']) == sha256(cfg['sparse'] if sparse else None),
            'relation LA dispatch mismatch')
    expected, seen, pivots, trajectory = [], set(), {}, []
    queries = duplicates = attempts = 0
    batches = report.get('matrix_batches', [])
    collection = report.get('collection_reports', [])
    require(len(batches) == len(collection), 'missing matrix batch observation')
    terminal = False
    for batch, observed in zip(collection, batches):
        require(not terminal, 'collection continued after certified LA success')
        for query in batch['attempts']:
            points = query['pdp']['points']
            novel = False
            if points is not None:
                key = query['a'], tuple(sorted(points))
                if key in seen:
                    duplicates += 1
                else:
                    seen.add(key)
                    novel = True
                    row = [0] * len(columns)
                    for index in points:
                        if projected[index] is not None:
                            j, coefficient = projected[index]
                            row[j] = (row[j] + coefficient) % c.r
                    expected.append(dict(entries=[[j, str(v)] for j, v in enumerate(row) if v],
                                         rhs=str(c.h * query['a'] % c.r)))
                    for j in range(len(columns)):
                        value = row[j]
                        if value == 0:
                            continue
                        if j in pivots:
                            row = [(a - value * b) % c.r for a, b in zip(row, pivots[j])]
                        else:
                            inverse = pow(value, -1, c.r)
                            pivots[j] = [v * inverse % c.r for v in row]
                            break
            queries += 1
            trajectory.append(dict(query=queries, novel_row=novel, rank=len(pivots)))
        attempts += len(expected) >= len(columns)
        for name, value in dict(queries=queries, accepted_rows=len(expected),
                                duplicate_relations=duplicates, rejected_relations=0,
                                solve_attempts=attempts).items():
            require(type(observed[name]) is int and observed[name] == value,
                    'matrix batch accounting mismatch: ' + name)
        require(type(observed['verified']) is bool, 'invalid batch certificate verdict')
        terminal = observed['verified']
        require(not terminal or len(pivots) == len(columns), 'certified batch lacks full rank')
        require((observed['sparse_report'] is not None) == (sparse and attempts > 0),
                'missing/unexpected sparse solve report')
        if sparse and attempts:
            diagnostic = observed['sparse_report']
            filtering = diagnostic['filter']
            require(all(type(diagnostic[name]) is int and diagnostic[name] >= 0 for name in
                        ('attempts', 'core_dimension', 'core_nonzeros', 'reconstructed_columns')),
                    'invalid sparse solve counters')
            require(filtering['rows_in'] == len(expected) and filtering['columns_in'] == len(columns),
                    'sparse input dimensions mismatch')
            require(all(type(value) is int and value >= 0 for value in filtering.values()),
                    'invalid sparse filter counters')
            require(0 <= diagnostic['attempts'] <= max(cfg['sparse']['attempts'], 1)
                    and 0 <= diagnostic['core_dimension'] <= len(columns), 'invalid sparse core counters')
            core = diagnostic['wiedemann']
            if core is not None:
                dimension = diagnostic['core_dimension'] + 1  # Homogenised RHS column.
                left = max(cfg['sparse']['wiedemann']['block_m'], 1)
                right = max(cfg['sparse']['wiedemann']['block_n'], 1)
                length = ((dimension+left-1)//left + (dimension+right-1)//right
                          + max(cfg['sparse']['wiedemann']['margin'], 2))
                require(diagnostic['attempts'] > 0 and core['dimension'] == dimension
                        and core['block_m'] == left and core['block_n'] == right
                        and core['sequence_length'] == length, 'block Wiedemann dispatch mismatch')
                require(all(type(value) is int and value >= 0 for value in core.values()),
                        'invalid block Wiedemann counters')
            covered = {entry[0] for row in expected for entry in row['entries']}
            if len(covered) < len(columns):
                require(filtering['uncovered_columns'] == len(columns)-len(covered)
                        and diagnostic['attempts'] == 0 and diagnostic['wiedemann'] is None,
                        'uncovered-column rejection ran an unreported core')
            else:
                require(filtering['nonzeros_in'] == sum(len(row['entries']) for row in expected),
                        'sparse input nonzero count mismatch')
            if terminal:
                require(diagnostic['reconstructed_columns'] + diagnostic['core_dimension'] == len(columns),
                        'sparse reconstruction column accounting mismatch')
        if not sparse:
            require(terminal == (len(pivots) == len(columns)), 'dense full-rank solve state mismatch')
    require(sha256(snapshot['rows']) == sha256(expected), 'stored matrix rows differ from query witnesses')
    if inventory:
        require(not expected and not batches and report['status'] == 'inventory', 'inventory ran collection')
    else:
        require(batches, 'missing matrix chronology')
        logs = report.get('column_logs')
        require(terminal == (logs is not None), 'terminal matrix certificate mismatch')
        final = report['log_table_report']
        for name, value in dict(relations=len(expected), columns=len(columns),
                                duplicate_relations=duplicates, rejected_relations=0,
                                solve_attempts=attempts, verified=terminal, sparse=sparse).items():
            require(type(final[name]) is type(value) and final[name] == value,
                    'final matrix accounting mismatch: ' + name)
        require(final['sparse_report'] == batches[-1]['sparse_report'], 'final sparse report changed')
        require(type(report['solve_attempts']) is int and report['solve_attempts'] == attempts,
                'top-level solve attempts mismatch')
        if terminal:
            require([entry['point'] for entry in logs] == columns, 'column log ordering mismatch')
            require(all(c.mul(c.g, int(entry['log'])) == p
                        and 0 <= int(entry['log']) < c.r for entry, p in zip(logs, representatives)),
                    'incorrect column log certificate')
        else:
            require(queries == cfg['max_trials'], 'uncertified collection stopped early')
    return dict(columns=len(columns), accepted_rows=len(expected), duplicate_relations=duplicates,
                rank=len(pivots), solve_attempts=attempts, certified_logs=terminal,
                matrix_sha256=sha256(snapshot), rank_trajectory=trajectory)


def verify_stages(report, fixture, job):
    require(job.get('exclusive_phases') is True and type(report.get('generic_admission_schema')) is int
            and report['generic_admission_schema'] == 1,
            'missing scientific stage evidence')
    require(report.get('generic_runtime_policy') == 'default-environment-one-rayon-v1',
            'missing runtime override policy')
    require(report['fixture'] == fixture and job.get('public_targets') == fixture['targets']
            and len(fixture['targets']) == 1, 'stage fixture mismatch')
    require(job['mode'] in {'ic', 'inventory'} and report['mode'] == 'ic', 'stage mode mismatch')
    cfg = effective_config(job, report)
    base = verify_base(report, fixture, job)
    inventory = job['mode'] == 'inventory'
    law = observed_dispatch = certificate = None
    if not inventory:
        law = verify_query_law(report, fixture, job)
        observed_dispatch = dispatch(report, cfg)
        if report['status'] == 'complete':
            certificate = verify(report, fixture, summands=cfg['summands'])
    matrix = matrix_audit(report, fixture, cfg, inventory=inventory)
    return dict(schema_version=1, status='PASS', base=base, query_law=law,
                dispatch=observed_dispatch, matrix=matrix, certificate=certificate,
                scope='base, query, dispatch, matrix and correctness; build/accounting required separately',
                promotion_eligible=False)
