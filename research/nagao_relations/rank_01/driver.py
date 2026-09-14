"""Frozen comparison and exhaustive truth checks for cached block spans."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REUSE = HERE.parent / 'reuse_01'
sys.path.insert(0, str(REUSE))
import trial
base = trial.base
import rank_cache


def cell(item, variant, mode, budget):
    if variant != 'cached-spans':
        return trial.cell(item, variant, mode, budget)
    if not item['target'][0]:
        row = trial.cell(item, 'reuse-batch-inverse', mode, budget)
        row['variant'] = variant
        return row
    start = time.perf_counter()
    f, c, v = base.old.context(item)
    target = tuple(f.fromCoords(x) for x in item['target'])
    phases = {'context_setup': time.perf_counter() - start}
    search, first, complete, solutions = None, None, False, {}
    try:
        search = rank_cache.Search(f, c, v, target, start + budget)
        for xs, points in search.solve(mode):
            solutions[xs] = points
            if first is None:
                first = time.perf_counter() - start
        complete = mode == 'enumerate' or not solutions
    except TimeoutError:
        pass
    elapsed = time.perf_counter() - start
    if search:
        phases.update(search.times)
    assert sum(phases.values()) <= elapsed + 1e-6
    phases['control_and_uninstrumented'] = max(0., elapsed - sum(phases.values()))
    return {'variant': variant, 'mode': mode,
            'status': 'first' if mode == 'first' and solutions else 'complete' if complete else 'timeout',
            'within_budget': elapsed <= budget, 'budget_seconds': budget, 'all_phase_seconds': elapsed,
            'first_verified_seconds': first, 'phase_seconds': phases,
            'solutions': [list(x) for x in sorted(solutions)],
            'certificates': [{'xs': list(x), 'points': solutions[x]} for x in sorted(solutions)],
            'rejection_certificates': search.certificates if search else [], 'stats': search.stats if search else {},
            'verified_unique_relations': len(solutions), 'counters': f.fullReport(), 'peak_memory_bytes': None,
            'common_operation_speedup': None, 'full_dlp_S': None, 'rho_ratio': None}


def blocks(row):
    return [{k: v for k, v in c.items() if k != 'separator'} for c in row['rejection_certificates']]


def validate():
    rng = random.Random(2026091482)
    fresh = base.old.space.independent(rng.sample(range(1, 64), 63))[:4]
    counts = {'target_space_pairs': 0, 'solver_trials': 0, 'relation_certificates': 0,
              'rejection_certificates': 0, 'identical_block_partitions': 0, 'failures': 0}
    for kind, basis in [('prefix', None), ('f4-stable', None), ('fresh', fresh)]:
        desc = base.descriptor(6, 4, kind, basis)
        f, c, v = base.old.context(desc, False)
        oracle = base.old.Oracle(f, c, v)
        truth = {}
        xs = sorted(oracle.base)
        for i, x in enumerate(xs):
            for j in range(i+1, len(xs)):
                for z in xs[j+1:]:
                    for mask in range(8):
                        total = None
                        for bit, xx in enumerate((x, xs[j], z)):
                            total = c.add(total, oracle.base[xx][mask >> bit & 1])
                        if total is not None and f.toCoords(total[0]) not in (x, xs[j], z):
                            truth.setdefault(total, set()).add((x, xs[j], z))
        for x in range(64):
            p = c.pointFromX(f.fromCoords(x))
            if p is None:
                continue
            for target in sorted({p, c.neg(p)}):
                assert oracle.expected(target) == truth.get(target, set())
                item = {**desc, 'target': [f.toCoords(x) for x in target]}
                checker = base.Checker(item)
                rows = [cell(item, variant, 'enumerate', 30) for variant in ('reuse-batch-inverse', 'cached-spans')]
                for row in rows:
                    assert row['status'] == 'complete'
                    base.check(row, item, oracle, checker)
                    counts['solver_trials'] += 1
                    counts['relation_certificates'] += len(row['certificates'])
                    counts['rejection_certificates'] += len(row['rejection_certificates'])
                assert rows[0]['solutions'] == rows[1]['solutions']
                assert blocks(rows[0]) == blocks(rows[1])
                for key in ('visited_admissible_leaves', 'pruned_admissible_branches', 'potential_branches'):
                    if rows[0].get('stats'):
                        assert rows[0]['stats'][key] == rows[1]['stats'][key]
                counts['target_space_pairs'] += 1
                counts['identical_block_partitions'] += 1
    return counts


def freeze():
    contract = json.loads((HERE / 'contract.json').read_text())
    assert base.digest(REUSE / 'targets.json') == contract['predecessor_targets_sha256']
    items = []
    for old in json.loads((REUSE / 'targets.json').read_text()):
        items.append({**{k: v for k, v in old.items() if k != 'instance_sha256'}, 'cohort': 'frozen-' + old['cohort']})
    rng = random.Random(contract['seed'])
    for n, d in ((18, 8), (18, 9), (30, 8), (30, 9)):
        desc = base.descriptor(n, d)
        f, c, v = base.old.context(desc, False)
        oracle = base.old.Oracle(f, c, v)
        for stratum in ('uniform', 'known_decomposable'):
            p = base.old.previous.uniform(f, c, rng) if stratum == 'uniform' else oracle.supported(rng)
            item = {**desc, 'target': [f.toCoords(x) for x in p], 'stratum': stratum, 'cohort': 'fresh-holdout'}
            assert not any(old['n'] == n and old['d'] == d and old['target'] == item['target'] for old in items)
            items.append(item)
    for item in items:
        item['instance_sha256'] = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()
    with (HERE / 'targets.json').open('x') as out:
        out.write(json.dumps(items, indent=2) + '\n')
    return {'targets': len(items), 'sha256': base.digest(HERE / 'targets.json')}


def provenance():
    paths = list(HERE.glob('*.py')) + [HERE / 'contract.json', HERE / 'targets.json']
    paths += [REUSE / name for name in ('trial.py', 'reuse.py', 'batch_trial.py')]
    paths += [trial.BLOCKS / name for name in ('block.py', 'experiment.py', 'batch.py')]
    paths += list(base.PRIOR.glob('*.py')) + list(base.old.OLD.glob('*.py')) + list(base.old.CODE.glob('*.py'))
    return {'kind': 'provenance', 'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'sha256': {str(p.relative_to(ROOT)): base.digest(p) for p in paths},
            'python': platform.python_version(), 'hardware': platform.uname()._asdict(), 'worker_threads': 1}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--freeze', action='store_true')
    p.add_argument('--validate-only', action='store_true')
    args = p.parse_args()
    if args.freeze:
        print(json.dumps(freeze()), flush=True)
        return
    if args.validate_only:
        print(json.dumps(validate()), flush=True)
        return
    contract = json.loads((HERE / 'contract.json').read_text())
    items = json.loads((HERE / 'targets.json').read_text())
    assert len(items)*len(contract['variants'])*len(contract['modes']) == contract['expected_trials']
    with (HERE / 'raw.jsonl').open('x') as out:
        def record(row):
            out.write(json.dumps(row) + '\n')
            out.flush()
        record(provenance())
        record({'kind': 'validation', **validate()})
        print('validation passed', flush=True)
        oracles = {}
        for item in items:
            key = item['n'], tuple(item['basis']), tuple(item['coefficients'])
            if key not in oracles:
                f, c, v = base.old.context(item, False)
                oracles[key] = base.old.Oracle(f, c, v)
            oracle, checker = oracles[key], base.Checker(item)
            expected = oracle.expected(tuple(oracle.f.fromCoords(x) for x in item['target']))
            record({'kind': 'instance', **item, 'expected': [list(x) for x in sorted(expected)]})
            for mode in contract['modes']:
                variants = contract['variants'][:]
                random.Random(item['instance_sha256'] + mode).shuffle(variants)
                for variant in variants:
                    row = cell(item, variant, mode, contract['cold_budget_seconds'])
                    t = time.perf_counter()
                    count = base.check(row, item, oracle, checker)
                    record({'kind': 'trial', **item, **row, 'expected_count': count, 'harness_audit_seconds': time.perf_counter() - t})
                    print(item['n'], item['d'], item['cohort'], mode, variant, row['status'], len(row['solutions']), round(row['all_phase_seconds'], 3), flush=True)


if __name__ == '__main__':
    main()
