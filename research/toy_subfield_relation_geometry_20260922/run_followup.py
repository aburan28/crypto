#!/usr/bin/env python3
"""A predeclared second toy population; keep the original experiment intact."""
import argparse
import hashlib
import json
from pathlib import Path
import random

from run import HERE, ToyCurve, ToyField, basis, inspect_base, random_basis


def largest_odd_prime_factor(value):
    while value % 2 == 0:
        value //= 2
    largest = 1
    divisor = 3
    while divisor * divisor <= value:
        while value % divisor == 0:
            largest = divisor
            value //= divisor
        divisor += 2
    return max(largest, value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('refusing to overwrite existing evidence')
    contract = json.loads((HERE / 'contract_followup.json').read_text())
    rows, selections, screen = [], [], []
    for n, k in contract['field_pairs_n_k']:
        f = ToyField(n)
        rng = random.Random(contract['seed'] + 100 * n + k)
        rejected = []
        for attempt in range(512):
            b = rng.randrange(1, f.size)
            if f.degree(b) != n:
                rejected.append({'b': b, 'reason': 'coefficient_in_proper_subfield'})
                continue
            curve = ToyCurve(f, b)
            prime = largest_odd_prime_factor(curve.order)
            if prime == 1 or curve.order > 8 * prime:
                rejected.append({'b': b, 'reason': 'odd_prime_factor_too_small', 'curve_order': curve.order, 'largest_odd_prime': prime})
                continue
            break
        else:
            raise RuntimeError(f'no admitted curve in frozen budget at n={n}')
        selections.append({'n': n, 'k': k, 'b': b, 'curve_order': curve.order,
                           'largest_odd_prime': prime, 'cofactor': curve.order // prime,
                           'rejected_candidates': rejected})
        u = basis(x for x in range(f.size) if f.frobenius(x, k) == x)
        gamma = next(x for x in range(2, f.size) if f.frobenius(x, k) != x)
        tower = basis(u + [f.mul(gamma, x) for x in u])
        for multiplier in contract['dimension_multipliers']:
            structured = u if multiplier == 1 else tower
            candidates = {
                'generic_linear': random_basis(f, multiplier * k, rng),
                'subfield_linear': structured,
                'scaled_subfield_linear': basis(f.mul(gamma, x) for x in structured),
            }
            for family, vectors in candidates.items():
                row = inspect_base(curve, vectors, k)
                row.update(n=n, k=k, extension_degree=n // k, curve_case='odd_prime_holdout',
                           curve_a=0, curve_b=b, coefficient_degree=f.degree(b), modulus=f.modulus,
                           curve_order=curve.order, largest_odd_prime=prime,
                           cofactor=curve.order // prime, family=family)
                rows.append(row)
        cell = [r for r in rows if r['n'] == n and r['x_dimension'] == k]
        reference = next(r for r in cell if r['family'] == 'generic_linear')
        passing = [r['family'] for r in cell if r['family'] != 'generic_linear'
                   and r['product_span_dimension'] < reference['product_span_dimension']
                   and r['distinct_nonidentity_targets'] >= reference['distinct_nonidentity_targets']]
        screen.append({'n': n, 'k': k, 'families_passing': passing})
    result = {'scope': contract['scope'], 'selections': selections, 'rows': rows,
              'holdout_screen': screen, 'screen_pass_cells': sum(bool(x['families_passing']) for x in screen),
              'screen_passed': sum(bool(x['families_passing']) for x in screen) >= 3,
              'source_hashes': {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                                for name in ('run.py', 'run_followup.py', 'contract.json', 'contract_followup.json')}}
    with args.output.open('x') as out:
        json.dump(result, out, indent=2, sort_keys=True)
        out.write('\n')
    print(json.dumps({k: result[k] for k in ('screen_pass_cells', 'screen_passed')}))


if __name__ == '__main__':
    main()
