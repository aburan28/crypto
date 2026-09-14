"""Complete matched eight-target batch with one charged setup per variant."""
import json
import random
import time
import driver
import batch_trial as prior

HERE, ROOT, base = driver.HERE, driver.ROOT, driver.base


def solve(items, variant, budget):
    if variant != 'cached-spans':
        return prior.solve(items, variant, budget)
    start = time.perf_counter()
    deadline = start + budget
    f, c, v = base.old.context(items[0])
    setup = time.perf_counter() - start
    queries = []
    for item in items:
        qstart = time.perf_counter()
        solutions, first, complete, search = {}, None, False, None
        attempted = qstart < deadline
        if attempted:
            target = tuple(f.fromCoords(x) for x in item['target'])
            try:
                search = driver.rank_cache.Search(f, c, v, target, deadline)
                for xs, points in search.solve('enumerate'):
                    solutions[xs] = points
                    if first is None:
                        first = time.perf_counter() - qstart
                complete = True
            except TimeoutError:
                pass
        queries.append({'target': item['target'], 'stratum': item['stratum'],
                        'status': 'complete' if complete else 'timeout' if attempted else 'unattempted',
                        'all_phase_seconds': time.perf_counter() - qstart,
                        'first_verified_seconds': first,
                        'solutions': [list(x) for x in sorted(solutions)],
                        'certificates': [{'xs': list(x), 'points': solutions[x]} for x in sorted(solutions)],
                        'rejection_certificates': search.certificates if search else [],
                        'stats': search.stats if search else {}, 'phase_seconds': search.times if search else {}})
    elapsed = time.perf_counter() - start
    query_seconds = sum(q['all_phase_seconds'] for q in queries)
    return {'variant': variant, 'mode': 'batch-enumerate', 'budget_seconds': budget,
            'status': 'complete' if all(q['status'] == 'complete' for q in queries) else 'incomplete',
            'within_budget': elapsed <= budget, 'all_phase_seconds': elapsed,
            'setup_seconds': setup, 'query_seconds': query_seconds,
            'orchestration_seconds': elapsed - setup - query_seconds,
            'queries': queries, 'counters': f.fullReport(), 'peak_memory_bytes': None,
            'common_operation_speedup': None, 'full_dlp_S': None, 'rho_ratio': None}


def main():
    contract = json.loads((HERE / 'contract.json').read_text())
    target_path = ROOT / contract['batch_targets']
    assert base.digest(target_path) == contract['batch_targets_sha256']
    items = json.loads(target_path.read_text())
    f, c, v = base.old.context(items[0], False)
    oracle = base.old.Oracle(f, c, v)
    checkers = [base.Checker(item) for item in items]
    variants = contract['batch_variants'][:]
    random.Random(contract['seed'] + 2).shuffle(variants)
    with (HERE / 'batch_raw.jsonl').open('x') as out:
        def record(row):
            out.write(json.dumps(row) + '\n')
            out.flush()
        provenance = driver.provenance()
        provenance['sha256'][str(target_path.relative_to(ROOT))] = base.digest(target_path)
        record(provenance)
        record({'kind': 'targets', 'items': items})
        for variant in variants:
            row = solve(items, variant, contract['batch_budget_seconds'])
            t = time.perf_counter()
            for item, query, checker in zip(items, row['queries'], checkers):
                query['expected_count'] = base.check(query, item, oracle, checker)
            record({'kind': 'batch', **row, 'harness_audit_seconds': time.perf_counter() - t})
            print(variant, row['status'], round(row['all_phase_seconds'], 3), row['counters']['field_api']['totals'], flush=True)


if __name__ == '__main__':
    main()
