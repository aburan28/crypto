#!/usr/bin/env python3
"""Exact product-span diagnostics. No curve/DLP cost extrapolation."""
import argparse
import hashlib
import json
import random
from pathlib import Path


def rem(a, f):
    n = f.bit_length()
    while a.bit_length() >= n:
        a ^= f << (a.bit_length() - n)
    return a


def mul(a, b, f):
    c = 0
    while b:
        bit = b & -b
        c ^= a << (bit.bit_length() - 1)
        b ^= bit
    return rem(c, f)


def square(a, f):
    c = 0
    while a:
        bit = a & -a
        c |= 1 << (2 * (bit.bit_length() - 1))
        a ^= bit
    return rem(c, f)


def gcd(a, b):
    while b:
        a, b = b, rem(a, b)
    return a


def irreducible(f):
    n = f.bit_length() - 1
    x = 2
    for i in range(1, n + 1):
        x = square(x, f)
        if i <= n // 2 and gcd(x ^ 2, f) != 1:
            return False
    return x == 2


def canonical(values):
    p = {}
    for v in values:
        while v:
            b = v.bit_length()
            if b not in p:
                p[b] = v
                break
            v ^= p[b]
    for b in sorted(p):
        for c in p:
            if c > b and p[c] >> (b - 1) & 1:
                p[c] ^= p[b]
    return tuple(p[b] for b in sorted(p))


def product_rank(a, b, f):
    n = f.bit_length() - 1
    pivots = {}
    for x in a:
        for y in b:
            z = mul(x, y, f)
            while z:
                k = z.bit_length()
                if k not in pivots:
                    pivots[k] = z
                    break
                z ^= pivots[k]
            if len(pivots) == n:
                return n
    return len(pivots)


def pow_field(a, e, f):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a, f)
        a = square(a, f)
        e >>= 1
    return r


def basis_for(f, ell, family, seed):
    n = f.bit_length() - 1
    if family == 'power_span':
        return [1 << i for i in range(ell)]
    if family == 'frobenius_span':
        out, x = [], 2
        for _ in range(n):
            if len(canonical(out + [x])) > len(out):
                out.append(x)
                if len(out) == ell:
                    return out
            x = square(x, f)
        return None
    if family == 'scaled_subfield':
        if ell == 1 or n % ell:
            return None
        out = []
        for v in range(1, 1 << n):
            x = pow_field(v, ((1 << n) - 1) // ((1 << ell) - 1), f)
            if len(canonical(out + [x])) > len(out):
                out.append(x)
                if len(out) == ell:
                    return [mul(2, x, f) for x in out]
    rng = random.Random(seed ^ 0x564543544F5253)
    out = []
    while len(out) < ell:
        x = rng.randrange(1, 1 << n)
        if len(canonical(out + [x])) > len(out):
            out.append(x)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('output', type=Path)
    args = ap.parse_args()
    contract = json.loads(Path(__file__).with_name('contract.json').read_text())
    rows = []
    fields = {}
    for n in contract['structural']['toy_degrees'] + [131]:
        if n == 131:
            f = (1 << 131) | (1 << 13) | 7
        else:
            f = next((1 << n) | lo for lo in range(1, 1 << n, 2)
                     if irreducible((1 << n) | lo))
        assert irreducible(f)
        fields[n] = {'polynomial_hex': hex(f), 'irreducible': True}
        dims = contract['structural']['large_dimensions'] if n == 131 else [max([d for d in range(2, 8) if n % d == 0] or [3])]
        for ell in dims:
            for family in contract['structural']['families']:
                for seed in contract['structural']['random_seeds']:
                    row = dict(n=n, ell=ell, family=family, seed=seed,
                               signed_orbit_count=None, projected_columns=None,
                               exact_pair_coverage=None, S=None, rho_ratio=None, floor_ratio=None)
                    basis = basis_for(f, ell, family, seed)
                    if basis is None:
                        row['unavailable'] = 'no subfield of requested dimension' if family == 'scaled_subfield' else 'seed orbit has insufficient rank'
                        rows.append(row)
                        continue
                    start = canonical(basis)
                    charts = []
                    current = basis
                    while True:
                        charts.append(current)
                        current = [square(x, f) for x in current]
                        if canonical(current) == start:
                            break
                        assert len(charts) < n
                    ranks = [product_rank(basis, b, f) for b in charts]
                    assert len(charts) == 1 or all(ranks[d] == ranks[-d] for d in range(1, len(charts)))
                    if family == 'scaled_subfield':
                        assert set(ranks) == {ell}
                    if n == 131:
                        assert min(ranks) >= min(n, 2 * ell - 1)
                    c = len(charts)
                    row.update(basis_hex=[hex(x) for x in basis], component_count=c,
                               relative_product_ranks=ranks, component_pairs=c * (c + 1) // 2,
                               pair_weighted_product_rank=sum((c - d) * r for d, r in enumerate(ranks)) / (c * (c + 1) / 2))
                    rows.append(row)
        print(f'completed field degree {n}', flush=True)
    out = dict(fields=fields, rows=rows, scope='field algebra only; curve census is in stage raw runs',
               contract_sha256=hashlib.sha256(Path(__file__).with_name('contract.json').read_bytes()).hexdigest())
    with args.output.open('x') as handle:
        json.dump(out, handle, indent=2)
        handle.write('\n')


if __name__ == '__main__':
    main()
