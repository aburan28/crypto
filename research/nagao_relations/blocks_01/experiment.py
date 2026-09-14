"""Frozen block-rejection experiment; audits are outside solver budgets."""
import argparse
import hashlib
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
PRIOR = HERE.parent / 'structured_01'
sys.path.insert(0, str(PRIOR))
import run as old
import block


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cell(instance, variant, mode, budget):
    if not variant.startswith('bilinear-'):
        return old.cell(instance, variant, mode, budget)
    if instance['target'][0] == 0:
        row = old.cell(instance, 'hybrid-reference', mode, budget)
        row.update(variant=variant, zero_chart_fallback=True, rejection_certificates=[])
        return row
    start = time.perf_counter()
    f, c, v = old.context(instance)
    target = tuple(f.fromCoords(x) for x in instance['target'])
    phases = {'context_setup': time.perf_counter() - start}
    search, first, complete = None, None, False
    solutions = {}
    try:
        search = block.Search(f, c, v, target, start + budget, variant == 'bilinear-blocks')
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
        raise ArithmeticError('overlapping phase timers')
    phases['control_and_uninstrumented'] = max(0., elapsed - sum(phases.values()))
    return {'variant': variant, 'mode': mode,
            'status': 'first' if mode == 'first' and solutions else 'complete' if complete else 'timeout',
            'within_budget': elapsed <= budget, 'budget_seconds': budget,
            'all_phase_seconds': elapsed, 'first_verified_seconds': first, 'phase_seconds': phases,
            'solutions': [list(x) for x in sorted(solutions)],
            'certificates': [{'xs': list(x), 'points': solutions[x]} for x in sorted(solutions)],
            'verified_unique_relations': len(solutions), 'zero_chart_fallback': False,
            'rejection_certificates': search.certificates if search else [],
            'stats': search.stats if search else {}, 'counters': f.fullReport(),
            'common_operations': None, 'full_dlp_S': None, 'rho_ratio': None, 'speedup': None}


class Checker:
    """Derive coefficients by finite differences, never Circuit's formula."""
    def __init__(self, instance):
        self.f, self.c, self.v = old.context(instance, False)
        self.r = self.f.fromCoords(instance['target'][0])
        self.cache = {}

    def coefficients(self, b):
        if b not in self.cache:
            f, c, v = self.f, self.c, self.v
            def direct(z, w):
                return block.direct(f, c, self.r, b, z, w)
            constant = direct(0, 0)
            zz = [direct(z, 0) ^ constant for z in v.basis]
            ww = [direct(0, w) ^ constant for w in v.basis]
            mixed = [[direct(z, w) ^ zz[i] ^ ww[j] ^ constant
                      for j, w in enumerate(v.basis)] for i, z in enumerate(v.basis)]
            self.cache[b] = constant, zz, ww, mixed
        return self.cache[b]

    def check(self, cert):
        f, v = self.f, self.v
        b = f.fromCoords(cert['b'])
        h = f.add(f.add(f.sqr(b), b), self.r)
        free, prefix, mask = cert['free_bits'], cert['prefix'], cert['separator']
        if not b or f.toCoords(h) != cert['h2'] or h not in v.indices:
            raise ArithmeticError('invalid coefficient in rejection')
        if not 0 < free <= v.d or not 0 <= prefix < 2**v.d or prefix & (2**free - 1):
            raise ArithmeticError('invalid block prefix')
        constant, zz, ww, mixed = self.coefficients(b)
        columns = ww[:]
        for i in range(free, v.d):
            if prefix >> i & 1:
                constant ^= zz[i]
                columns = [a ^ b for a, b in zip(columns, mixed[i])]
        columns += zz[:free]
        columns += [x for row in mixed[:free] for x in row]
        if (constant & mask).bit_count() % 2 != 1 or any((x & mask).bit_count() % 2 for x in columns):
            raise ArithmeticError('false block rejection certificate')


def check(row, instance, oracle, checker):
    expected_count = old.check(row, instance, oracle)
    if {tuple(c['xs']) for c in row['certificates']} != {tuple(x) for x in row['solutions']}:
        raise ArithmeticError('missing relation certificate')
    for cert in row.get('rejection_certificates', []):
        checker.check(cert)
    if row.get('stats'):
        stats = row['stats']
        if stats['rejected_blocks'] != len(row['rejection_certificates']):
            raise ArithmeticError('missing rejection certificate')
        if row['status'] == 'complete' and stats['visited_admissible_leaves'] + stats['pruned_admissible_branches'] != stats['potential_branches']:
            raise ArithmeticError('incomplete branch partition')
    return expected_count


def descriptor(n, d, kind='prefix', basis=None):
    f = old.field.Onb(n)
    alpha = old.algebra.subfieldGenerator(f, 2)
    if basis is None:
        basis, _ = old.space.makeBasis(f, d, kind, alpha)
    return {'n': n, 'd': d, 'base_kind': kind, 'basis': basis,
            'coefficients': [f.toCoords(alpha)] * 2}


def freeze():
    contract = json.loads((HERE / 'contract.json').read_text())
    source = ROOT / contract['predecessor']
    if digest(source) != contract['predecessor_sha256']:
        raise ArithmeticError('predecessor hash mismatch')
    items = []
    for line in source.read_text().splitlines():
        row = json.loads(line)
        if row['kind'] == 'instance' and row['base_kind'] == 'prefix' and row['d'] == 8:
            item = {k: row[k] for k in ('n', 'd', 'base_kind', 'basis', 'coefficients', 'target', 'stratum')}
            item['cohort'] = 'frozen-' + row['cohort']
            items.append(item)
    if len(items) != 8:
        raise ArithmeticError('wrong frozen input count')
    rng = random.Random(contract['seed'])
    for n, d in ((18, 8), (18, 9), (30, 8), (30, 9)):
        desc = descriptor(n, d)
        f, c, v = old.context(desc, False)
        oracle = old.Oracle(f, c, v)
        for stratum in ('uniform', 'known_decomposable'):
            p = old.previous.uniform(f, c, rng) if stratum == 'uniform' else oracle.supported(rng)
            items.append({**desc, 'target': [f.toCoords(x) for x in p],
                          'cohort': 'fresh-holdout', 'stratum': stratum})
    for item in items:
        item['instance_sha256'] = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()
    with (HERE / 'targets.json').open('x') as out:
        json.dump(items, out, indent=2)
        out.write('\n')
    return {'frozen': 8, 'fresh': 8, 'targets_sha256': digest(HERE / 'targets.json')}


def validate():
    start = time.perf_counter()
    counts = {'exhaustive_target_space_pairs': 0, 'new_solver_trials': 0,
              'rejection_certificates': 0, 'circuit_identity_evaluations': 0, 'relation_certificates': 0}
    rng = random.Random(2026091448)
    fresh = old.space.independent(rng.sample(range(1, 64), 63))[:4]
    for kind, basis in [('prefix', None), ('f4-stable', None), ('fresh', fresh)]:
        desc = descriptor(6, 4, kind, basis)
        f, c, v = old.context(desc, False)
        oracle = old.Oracle(f, c, v)
        truth = {}
        xs = sorted(oracle.base)
        for i, x in enumerate(xs):
            for j in range(i + 1, len(xs)):
                for z in xs[j + 1:]:
                    for mask in range(8):
                        total = None
                        for bit, xx in enumerate((x, xs[j], z)):
                            total = c.add(total, oracle.base[xx][mask >> bit & 1])
                        if total is not None and f.toCoords(total[0]) not in (x, xs[j], z):
                            truth.setdefault(total, set()).add((x, xs[j], z))
        circuit_targets = 0
        for x in range(64):
            p = c.pointFromX(f.fromCoords(x))
            if p is None:
                continue
            for target in sorted({p, c.neg(p)}):
                if oracle.expected(target) != truth.get(target, set()):
                    raise ArithmeticError('pair oracle versus exhaustive triple mismatch')
                instance = {**desc, 'target': [f.toCoords(a) for a in target]}
                checker = Checker(instance)
                for variant in ('bilinear-reference', 'bilinear-blocks'):
                    row = cell(instance, variant, 'enumerate', 30)
                    if row['status'] != 'complete':
                        raise ArithmeticError('tiny validation incomplete')
                    check(row, instance, oracle, checker)
                    counts['new_solver_trials'] += 1
                    counts['rejection_certificates'] += len(row['rejection_certificates'])
                    counts['relation_certificates'] += len(row['certificates'])
                counts['exhaustive_target_space_pairs'] += 1
                if target[0] and circuit_targets < 2:
                    powers = [(w, f.sqr(w), f.frob(w, 2)) for w in v.basis]
                    for b in [f.fromCoords(i) for i in (1, 2, 3, 7)]:
                        h = f.add(f.add(f.sqr(b), b), target[0])
                        circuit = block.Circuit(f, c, v, target[0], b, h, powers)
                        for zi, z in enumerate(v.values):
                            for wi, w in enumerate(v.values):
                                if circuit.evaluate(zi, wi) != block.direct(f, c, target[0], b, z, w):
                                    raise ArithmeticError('small bilinear identity mismatch')
                                counts['circuit_identity_evaluations'] += 1
                    circuit_targets += 1
    for n in (18, 30):
        desc = descriptor(n, 9)
        f, c, v = old.context(desc, False)
        powers = [(w, f.sqr(w), f.frob(w, 2)) for w in v.basis]
        for _ in range(16):
            r, b = [f.fromCoords(rng.randrange(1, 2**n)) for _ in range(2)]
            h = f.add(f.add(f.sqr(b), b), r)
            circuit = block.Circuit(f, c, v, r, b, h, powers)
            for _ in range(32):
                zi, wi = [rng.randrange(2**v.d) for _ in range(2)]
                if circuit.evaluate(zi, wi) != block.direct(f, c, r, b, v.values[zi], v.values[wi]):
                    raise ArithmeticError('large bilinear identity mismatch')
                counts['circuit_identity_evaluations'] += 1
    return {**counts, 'failures': 0, 'seconds': time.perf_counter() - start}


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
        raise ArithmeticError('unexpected trial count')
    paths = list(HERE.glob('*.py')) + [HERE / 'contract.json', HERE / 'targets.json']
    paths += list(PRIOR.glob('*.py')) + list(old.OLD.glob('*.py')) + list(old.CODE.glob('*.py'))
    provenance = {'kind': 'provenance',
                  'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  'sha256': {str(p.relative_to(ROOT)): digest(p) for p in paths},
                  'python': platform.python_version(), 'pycryptosat': importlib.metadata.version('pycryptosat'),
                  'hardware': platform.uname()._asdict(), 'worker_threads': 1}
    with (HERE / 'raw.jsonl').open('x') as out:
        def record(row):
            out.write(json.dumps(row) + '\n')
            out.flush()
        record(provenance)
        record({'kind': 'validation', **validate()})
        print('validation passed', flush=True)
        for instance in items:
            f, c, v = old.context(instance, False)
            t = time.perf_counter()
            oracle = old.Oracle(f, c, v)
            checker = Checker(instance)
            target = tuple(f.fromCoords(x) for x in instance['target'])
            expected = oracle.expected(target)
            record({'kind': 'instance', **instance, 'expected': [list(x) for x in sorted(expected)],
                    'oracle_setup_seconds': time.perf_counter() - t})
            for mode in contract['modes']:
                variants = contract['variants'][:]
                random.Random(instance['instance_sha256'] + mode).shuffle(variants)
                for variant in variants:
                    row = cell(instance, variant, mode, contract['cold_budget_seconds'])
                    t = time.perf_counter()
                    count = check(row, instance, oracle, checker)
                    record({'kind': 'trial', **instance, **row, 'expected_count': count,
                            'harness_audit_seconds': time.perf_counter() - t})
                    print(instance['n'], instance['d'], instance['cohort'], instance['stratum'],
                          mode, variant, row['status'], len(row['solutions']),
                          round(row['all_phase_seconds'], 3), flush=True)


if __name__ == '__main__':
    main()
