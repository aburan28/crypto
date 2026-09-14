"""Conditional eight-target follow-up with explicitly charged shared setup."""
import argparse
import importlib.metadata
import json
import platform
import random
import subprocess
import time
import experiment as exp

HERE, ROOT, old = exp.HERE, exp.ROOT, exp.old


def freeze():
    screening = json.loads((HERE / 'summary.json').read_text())
    if not screening['batch_gate']['passed']:
        raise ArithmeticError('predeclared cold screen failed; do not extend candidate')
    contract = json.loads((HERE / 'batch_contract.json').read_text())
    desc = exp.descriptor(contract['n'], contract['d'])
    f, c, v = old.context(desc, False)
    oracle = old.Oracle(f, c, v)
    rng = random.Random(contract['seed'])
    items = []
    prior_targets = {tuple(item['target']) for item in json.loads((HERE / 'targets.json').read_text())}
    for stratum in ('uniform', 'known_decomposable'):
        for _ in range(contract['targets_per_stratum']):
            p = old.previous.uniform(f, c, rng) if stratum == 'uniform' else oracle.supported(rng)
            item = {**desc, 'target': [f.toCoords(x) for x in p], 'stratum': stratum, 'cohort': 'fresh-batch-holdout'}
            if tuple(item['target']) in prior_targets:
                raise ArithmeticError('batch target unexpectedly repeats cold input; record before redesign')
            items.append(item)
    with (HERE / 'batch_targets.json').open('x') as out:
        json.dump(items, out, indent=2)
        out.write('\n')
    return {'targets': len(items), 'targets_sha256': exp.digest(HERE / 'batch_targets.json')}


def solve(items, variant, budget):
    start = time.perf_counter()
    deadline = start + budget
    f, c, v = old.context(items[0])
    table = cache = None
    queries = []
    setup_complete = False
    try:
        if variant == 'pair-invariants-s3':
            table = old.solvers.PairTable(f, c, v, deadline)
        elif variant == 'hybrid-filtered':
            cache = old.solvers.ImageCache(f, v, True, deadline)
        setup_complete = True
    except TimeoutError:
        pass
    setup_seconds = time.perf_counter() - start
    for item in items:
        qstart = time.perf_counter()
        if not setup_complete or qstart >= deadline:
            queries.append({'target': item['target'], 'stratum': item['stratum'], 'status': 'unattempted',
                            'solutions': [], 'certificates': [], 'rejection_certificates': [],
                            'all_phase_seconds': 0., 'first_verified_seconds': None})
            continue
        target = tuple(f.fromCoords(x) for x in item['target'])
        solutions, first, complete, search = {}, None, False, None
        try:
            if variant.startswith('bilinear-'):
                if not target[0]:
                    raise ValueError('batch contract uses nonzero target chart')
                search = exp.block.Search(f, c, v, target, deadline, variant == 'bilinear-blocks')
                for xs, points in search.solve('enumerate'):
                    solutions[xs] = points
                    if first is None:
                        first = time.perf_counter() - qstart
            elif variant == 'pair-invariants-s3':
                f.phase = 'queries_and_verification'
                for xs, points in table.solve(target, 'enumerate', deadline):
                    solutions[xs] = points
                    if first is None:
                        first = time.perf_counter() - qstart
            else:
                f.phase = 'query_setup'
                search = old.solvers.Search(f, c, v, target, deadline, True, cache)
                seen = set()
                iterator = search.candidates()
                while True:
                    f.phase = 'candidate_search'
                    try:
                        a, b, invB, z = next(iterator)
                    except StopIteration:
                        break
                    if (a, b) in seen:
                        continue
                    seen.add((a, b))
                    f.phase = 'extraction_verification'
                    result = search.recover(a, b, invB, z)
                    if result:
                        solutions[result[0]] = result[1]
                        if first is None:
                            first = time.perf_counter() - qstart
            complete = True
        except TimeoutError:
            pass
        queries.append({'target': item['target'], 'stratum': item['stratum'],
                        'status': 'complete' if complete else 'timeout',
                        'all_phase_seconds': time.perf_counter() - qstart, 'first_verified_seconds': first,
                        'solutions': [list(x) for x in sorted(solutions)],
                        'certificates': [{'xs': list(x), 'points': solutions[x]} for x in sorted(solutions)],
                        'rejection_certificates': search.certificates if isinstance(search, exp.block.Search) else [],
                        'stats': search.stats if isinstance(search, exp.block.Search) else None})
    elapsed = time.perf_counter() - start
    return {'variant': variant, 'mode': 'batch-enumerate', 'budget_seconds': budget,
            'status': 'complete' if all(q['status'] == 'complete' for q in queries) else 'incomplete',
            'within_budget': elapsed <= budget, 'all_phase_seconds': elapsed, 'setup_seconds': setup_seconds,
            'query_seconds': sum(q['all_phase_seconds'] for q in queries),
            'orchestration_seconds': elapsed - setup_seconds - sum(q['all_phase_seconds'] for q in queries),
            'queries': queries, 'counters': f.fullReport(), 'table': table.metadata() if table else None,
            'common_operation_speedup': None, 'full_dlp_S': None, 'rho_ratio': None}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    if args.freeze:
        print(json.dumps(freeze()), flush=True)
        return
    contract = json.loads((HERE / 'batch_contract.json').read_text())
    items = json.loads((HERE / 'batch_targets.json').read_text())
    if not json.loads((HERE / 'summary.json').read_text())['batch_gate']['passed']:
        raise ArithmeticError('no follow-up authorization from frozen screen')
    paths = [HERE / p for p in ('batch.py', 'batch_contract.json', 'batch_targets.json', 'block.py', 'experiment.py')]
    provenance = {'kind': 'provenance',
                  'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  'sha256': {str(p.relative_to(ROOT)): exp.digest(p) for p in paths},
                  'cold_raw_sha256': exp.digest(HERE / 'raw.jsonl'),
                  'python': platform.python_version(), 'pycryptosat': importlib.metadata.version('pycryptosat'),
                  'hardware': platform.uname()._asdict(), 'worker_threads': 1}
    f, c, v = old.context(items[0], False)
    oracle = old.Oracle(f, c, v)
    checkers = [exp.Checker(item) for item in items]
    variants = contract['variants'][:]
    random.Random(contract['seed'] + 1).shuffle(variants)
    with (HERE / 'batch_raw.jsonl').open('x') as out:
        def record(row):
            out.write(json.dumps(row) + '\n')
            out.flush()
        record(provenance)
        record({'kind': 'targets', 'items': items})
        for variant in variants:
            row = solve(items, variant, contract['batch_budget_seconds'])
            t = time.perf_counter()
            for item, query, checker in zip(items, row['queries'], checkers):
                query['expected_count'] = exp.check(query, item, oracle, checker)
            record({'kind': 'batch', **row, 'harness_audit_seconds': time.perf_counter() - t})
            print(variant, row['status'], round(row['all_phase_seconds'], 3),
                  sum(len(q['solutions']) for q in row['queries']), flush=True)


if __name__ == '__main__':
    main()
