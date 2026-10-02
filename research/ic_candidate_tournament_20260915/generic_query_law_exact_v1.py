"""Versioned query-law admission with explicit exact three-sum negatives.

The pinned RNG implementation is reused; every query, budget and chronology
check from the legacy law remains. Closed historical evaluators are unchanged.
"""
import hashlib
from pathlib import Path
from generic_query_law import (SOLVERS, STREAM_VERSION, canonical_digest,
                               collection_coefficients, descent_coefficients, uint)
from generic_queries_exact_v1 import verify_queries
from oracle import require


def verify_query_law(report, fixture, job):
    """Check all recorded queries against the externally frozen worker job.

    The ordinary group/ledger checker runs first. A passing receipt does not
    establish a rank-based stop condition or certify an implementation/source.
    """
    require(job['mode'] == 'ic', 'query-law adapter requires IC')
    degree = uint(job['degree'], 32, 'degree')
    require(5 <= degree <= 31 and degree % 2 == 1, 'query-law adapter degree unsupported')
    a = uint(job['curve_a'], 8, 'curve coefficient')
    require(a <= 1 and degree == fixture['degree'] and a == fixture['curve_a'],
            'query-law curve mismatch')
    require(job.get('public_targets') == fixture['targets'] and len(fixture['targets']) == 1,
            'query-law requires the same one supplied point')
    seed = uint(job['algorithm_seed'], 64, 'algorithm seed')
    config = job['config']
    solver = config.get('solver', 'pair_table')
    require(solver in SOLVERS, 'query-law solver unsupported')
    summands = uint(config.get('summands', 3), 64, 'summands')
    require(2 <= summands <= 4 and report['summands'] == summands, 'query-law summands mismatch')
    batch_size = uint(config.get('batch_trials', 64), 64, 'batch size')
    limit = uint(config.get('max_trials', 4096), 64, 'trial limit')
    require(1 <= batch_size <= min(limit, 4096) and limit <= 65536, 'query-law limits unsupported')
    window = config.get('collection_window')
    if window is not None:
        uint(window, 64, 'collection window')
    order = int(fixture['subgroup_order'])
    group_receipt = verify_queries(report, fixture, summands)
    # Every supported bounded field has the packed backend. A pair-table worker
    # cannot start without its pair table. The actual geometric point count
    # controls the window predicate (not orbit columns or the usable census B).
    walked_collection = (solver == 'pair_table' and summands == 3 and window is not None
                         and 0 < window < len(report['factor_base']))
    expected = collection_coefficients(seed, order, walked_collection)
    trial = 0
    for batch in report['collection_reports']:
        require(trial < limit and batch['trials'] == min(batch_size, limit - trial),
                'query-law batch partition mismatch')
        for attempt in batch['attempts']:
            require((attempt['a'], attempt['b']) == next(expected),
                    f'query-law collection coefficient mismatch at {trial}')
            trial += 1
    require(0 < trial <= limit, 'query-law collection budget mismatch')
    walked_descent = solver == 'pair_table'
    solutions = report.get('solutions', [])
    require(report['status'] in {'complete', 'incomplete'}, 'query-law invalid terminal status')
    require(len(solutions) <= 1, 'query-law unexpected target count')
    if report['status'] == 'complete':
        require(len(solutions) == 1, 'query-law completed target missing')
    elif not solutions:
        require(trial == limit, 'query-law failed preparation stopped early')
    for solution in solutions:
        require(solution['index'] == 0 and 0 < solution['trials'] <= limit,
                'query-law descent budget mismatch')
        expected = descent_coefficients(seed, order, walked_descent)
        for t, attempt in enumerate(solution['attempts']):
            require((attempt['a'], attempt['b']) == next(expected),
                    f'query-law descent coefficient mismatch at {t}')
        if solution['recovered'] is None:
            require(solution['trials'] == limit, 'query-law failed descent stopped early')
    return dict(schema_version=1, stream_version=STREAM_VERSION,
                collection_law='windowed-64' if walked_collection else 'sampled',
                descent_law='walked-64' if walked_descent else 'sampled',
                job_sha256=canonical_digest(job), query_sha256=group_receipt['query_sha256'],
                checker_sha256={name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in ('generic_query_law_exact_v1.py', 'generic_queries_exact_v1.py',
                                             'generic_exact_negative_v1.py', 'generic_query_law.py', 'oracle.py')},
                collection_queries=trial, descent_queries=group_receipt['descent_queries'],
                negative_proofs=group_receipt['negative_proofs'],
                scope='query law, bounded exact negatives and group replay only', promotion_eligible=False)
