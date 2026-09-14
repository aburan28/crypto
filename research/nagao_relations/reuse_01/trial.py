"""Matched coefficient-reuse trials with frozen order and rejection rules."""
import argparse
import gzip
import hashlib
import importlib.util
import importlib.metadata
import json
from pathlib import Path
import platform
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BLOCKS = HERE.parent / 'blocks_01'
sys.path.insert(0, str(BLOCKS))
spec = importlib.util.spec_from_file_location('blocks_experiment', BLOCKS / 'experiment.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
import reuse


def cell(instance, variant, mode, budget):
    if not variant.startswith('reuse-'):
        return base.cell(instance, variant, mode, budget)
    if not instance['target'][0]:
        row = base.cell(instance, 'bilinear-blocks', mode, budget)
        row['variant'] = variant
        return row
    start = time.perf_counter()
    f, c, v = base.old.context(instance)
    target = tuple(f.fromCoords(x) for x in instance['target'])
    phases = {'context_setup': time.perf_counter() - start}
    search, first, complete, solutions = None, None, False, {}
    try:
        search = reuse.Search(f, c, v, target, start + budget, variant == 'reuse-batch-inverse')
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
    if sum(phases.values()) > elapsed + 1e-6:
        raise ArithmeticError('reuse phase overlap')
    phases['control_and_uninstrumented'] = max(0., elapsed - sum(phases.values()))
    return {'variant': variant, 'mode': mode,
            'status': 'first' if mode == 'first' and solutions else 'complete' if complete else 'timeout',
            'within_budget': elapsed <= budget, 'budget_seconds': budget, 'all_phase_seconds': elapsed,
            'first_verified_seconds': first, 'phase_seconds': phases,
            'solutions': [list(x) for x in sorted(solutions)],
            'certificates': [{'xs': list(x), 'points': solutions[x]} for x in sorted(solutions)],
            'rejection_certificates': search.certificates if search else [],
            'stats': search.stats if search else {}, 'verified_unique_relations': len(solutions),
            'counters': f.fullReport(), 'peak_memory_bytes': None,
            'common_operation_speedup': None, 'full_dlp_S': None, 'rho_ratio': None}


def freeze():
    contract = json.loads((HERE / 'contract.json').read_text())
    data = gzip.decompress((ROOT / contract['predecessor']).read_bytes())
    if hashlib.sha256(data).hexdigest() != contract['predecessor_uncompressed_sha256']:
        raise ArithmeticError('predecessor hash differs')
    items = []
    for row in map(json.loads, data.splitlines()):
        if row['kind'] == 'instance':
            item = {k: row[k] for k in ('n', 'd', 'base_kind', 'basis', 'coefficients', 'target', 'stratum')}
            item['cohort'] = 'frozen-' + row['cohort']
            items.append(item)
    if len(items) != 16:
        raise ArithmeticError('incomplete frozen input coverage')
    rng = random.Random(contract['seed'])
    for n, d in ((18, 8), (18, 9), (30, 8), (30, 9)):
        desc = base.descriptor(n, d)
        f, c, v = base.old.context(desc, False)
        oracle = base.old.Oracle(f, c, v)
        for stratum in ('uniform', 'known_decomposable'):
            p = base.old.previous.uniform(f, c, rng) if stratum == 'uniform' else oracle.supported(rng)
            item = {**desc, 'target': [f.toCoords(x) for x in p], 'stratum': stratum, 'cohort': 'fresh-holdout'}
            if any(item['n'] == old['n'] and item['d'] == old['d'] and item['target'] == old['target'] for old in items):
                raise ArithmeticError('fresh target repeats predecessor; record before redesign')
            items.append(item)
    for item in items:
        item['instance_sha256'] = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()
    with (HERE / 'targets.json').open('x') as out:
        json.dump(items, out, indent=2)
        out.write('\n')
    return {'frozen': 16, 'fresh': 8, 'targets_sha256': base.digest(HERE / 'targets.json')}


def validate():
    prior = base.validate()
    counts = {'solver_trials': 0, 'relation_certificates': 0, 'rejection_certificates': 0, 'interpolated_circuits': 0}
    rng = random.Random(2026091474)
    fresh = base.old.space.independent(rng.sample(range(1, 64), 63))[:4]
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
        checked_tables = 0
        for x in range(64):
            point = c.pointFromX(f.fromCoords(x))
            if point is None:
                continue
            for target in sorted({point, c.neg(point)}):
                if oracle.expected(target) != truth.get(target, set()):
                    raise ArithmeticError('reuse validation oracle differs from exhaustive triples')
                item = {**desc, 'target': [f.toCoords(x) for x in target]}
                checker = base.Checker(item)
                for variant in ('reuse-circuit', 'reuse-batch-inverse'):
                    row = cell(item, variant, 'enumerate', 30)
                    if row['status'] != 'complete':
                        raise ArithmeticError('tiny reuse search incomplete')
                    base.check(row, item, oracle, checker)
                    counts['solver_trials'] += 1
                    counts['relation_certificates'] += len(row['certificates'])
                    counts['rejection_certificates'] += len(row['rejection_certificates'])
                if target[0] and checked_tables < 2:
                    counts['interpolated_circuits'] += check_table(item, None)
                    checked_tables += 1
    for n, d in ((18, 8), (18, 9), (30, 8), (30, 9)):
        desc = base.descriptor(n, d)
        f, c, v = base.old.context(desc, False)
        p = base.old.previous.uniform(f, c, rng)
        counts['interpolated_circuits'] += check_table({**desc, 'target': [f.toCoords(x) for x in p]}, 32)
    return {'predecessor': prior, **counts, 'failures': 0}


def check_table(item, limit):
    f, c, v = base.old.context(item)
    target = tuple(f.fromCoords(x) for x in item['target'])
    search = reuse.Search(f, c, v, target, time.perf_counter() + 30, True)
    search.prepare()
    entries = list(search.cache.items())
    if limit is not None:
        random.Random(2026091475).shuffle(entries)
        entries = entries[:limit]
    for b, vector in entries:
        h = f.add(f.add(f.sqr(b), b), target[0])
        direct = base.block.Circuit(f, c, v, target[0], b, h, search.powers)
        if vector != reuse.flatten(direct):
            raise ArithmeticError('quadratic interpolation differs from direct circuit')
        if b and f.mul(b, search.inverses[b]) != f.one():
            raise ArithmeticError('batch inverse differs')
    return len(entries)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    if args.freeze:
        print(json.dumps(freeze()), flush=True)
        return
    if args.validate_only:
        print(json.dumps(validate()), flush=True)
        return
    contract = json.loads((HERE / 'contract.json').read_text())
    items = json.loads((HERE / 'targets.json').read_text())
    if len(items) * len(contract['variants']) * len(contract['modes']) != contract['expected_trials']:
        raise ArithmeticError('trial count mismatch')
    paths = list(HERE.glob('*.py')) + [HERE / 'contract.json', HERE / 'targets.json']
    paths += [BLOCKS / name for name in ('block.py', 'experiment.py')]
    paths += list(base.PRIOR.glob('*.py')) + list(base.old.OLD.glob('*.py')) + list(base.old.CODE.glob('*.py'))
    provenance = {'kind': 'provenance',
                  'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  'sha256': {str(p.relative_to(ROOT)): base.digest(p) for p in paths},
                  'python': platform.python_version(), 'pycryptosat': importlib.metadata.version('pycryptosat'),
                  'hardware': platform.uname()._asdict(), 'worker_threads': 1}
    with (HERE / 'raw.jsonl').open('x') as out:
        def record(row):
            out.write(json.dumps(row) + '\n')
            out.flush()
        record(provenance)
        record({'kind': 'validation', **validate()})
        print('validation passed', flush=True)
        for item in items:
            f, c, v = base.old.context(item, False)
            oracle = base.old.Oracle(f, c, v)
            checker = base.Checker(item)
            expected = oracle.expected(tuple(f.fromCoords(x) for x in item['target']))
            record({'kind': 'instance', **item, 'expected': [list(x) for x in sorted(expected)]})
            for mode in contract['modes']:
                variants = contract['variants'][:]
                random.Random(item['instance_sha256'] + mode).shuffle(variants)
                for variant in variants:
                    row = cell(item, variant, mode, contract['cold_budget_seconds'])
                    t = time.perf_counter()
                    count = base.check(row, item, oracle, checker)
                    record({'kind': 'trial', **item, **row, 'expected_count': count,
                            'harness_audit_seconds': time.perf_counter() - t})
                    print(item['n'], item['d'], item['cohort'], item['stratum'], mode, variant,
                          row['status'], len(row['solutions']), round(row['all_phase_seconds'], 3), flush=True)


if __name__ == '__main__':
    main()
