#!/usr/bin/env python3
"""Independently census rational-support spans on every frozen base."""
import gzip
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('field', HERE.parent/'weil_factor_composition_20260914/structure.py')
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)


def frob(x, k, poly):
    for _ in range(k):
        x = f.square(x, poly)
    return x


def kernel(n, k, poly):
    pivots, out = [], []
    for i in range(n):
        img, pre = frob(1 << i, k, poly) ^ (1 << i), 1 << i
        for pi, pp in pivots:
            if img & (1 << (pi.bit_length()-1)):
                img, pre = img ^ pi, pre ^ pp
        if img:
            pivots.append((img, pre))
            pivots.sort(reverse=True)
        else:
            out.append(pre)
    return out


def combine(basis, mask):
    out = 0
    for i, b in enumerate(basis):
        if mask >> i & 1:
            out ^= b
    return out


def main():
    setups = {}
    for line in gzip.open(HERE/'results/iteration-2/processes.jsonl.gz', 'rt'):
        r = json.loads(line)
        if r['kind'] == 'dlp' and r['revision'] == 'candidate':
            s = next(json.loads(x) for x in r['stdout'].splitlines() if json.loads(x)['phase'] == 'setup')
            setups[(r['case'],r['seed'])] = s
    out = []
    for (ci, seed), s in sorted(setups.items()):
        c = s['configuration']
        n, k = c['n'], c['k']
        poly = (1 << n) | sum(1 << i for i in s['field_polynomial'])
        field_basis = kernel(n, k, poly)
        a, b = combine(field_basis, c['a']), combine(field_basis, c['b'])

        def lifts(x):
            if x == 0:
                return True
            rhs = x ^ a ^ f.mul(b, f.square(f.pow_field(x, (1 << n)-2, poly), poly), poly)
            tr, acc = rhs, rhs
            for _ in range(1, n):
                acc = f.square(acc, poly)
                tr ^= acc
            assert tr in (0, 1)
            return tr == 0

        basis = s['basis']
        seen, xs, dimensions, populations = set(), set(), [], []
        while f.canonical(basis) not in seen:
            seen.add(f.canonical(basis))
            support = []
            for mask in range(1 << len(basis)):
                x = combine(basis, mask)
                xs.add(x)
                if lifts(x):
                    support.append(mask)
            dimensions.append(len(f.canonical(support)))
            populations.append(len(support))
            basis = [frob(x, k, poly) for x in basis]
        points = sum(1 if x == 0 else 2 for x in xs if lifts(x))
        assert points == s['points'], (ci, seed, points, s['points'])
        out.append(dict(case=ci, seed=seed, nominal_dimension=c['ell'], support_dimensions=dimensions,
                        rational_abscissae_per_component=populations, factor_base_points=points))
    (HERE/'support-census.json').write_text(json.dumps(out, indent=2)+'\n')
    for r in out:
        if r['seed'] == 101:
            print(r['case'], r['nominal_dimension'], r['support_dimensions'], r['rational_abscissae_per_component'])


if __name__ == '__main__':
    main()
