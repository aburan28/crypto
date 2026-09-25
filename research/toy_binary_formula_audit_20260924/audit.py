"""Bounded exhaustive algebra audit over GF(8), GF(16), GF(32).

No discrete-log solver; no large-field mode. Python standard library only.
Run: python3 audit.py --output results.json
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path


def multiply(a, b, modulus, n):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & (1 << n):
            a ^= modulus
    return out


def span(basis):
    values = [0]
    for b in basis:
        values += [v ^ b for v in values]
    return values


def anf(table):
    coefficients = list(table)
    for bit in range((len(table) - 1).bit_length()):
        for mask in range(len(table)):
            if mask & (1 << bit):
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return coefficients


def run():
    records = []
    checks = dict(field_inverse=0, anf_roundtrip=0,
                  constants_only=0, frobenius_covariance=0,
                  permutation_symmetry=0)
    for n, modulus in [(3, 0b1011), (4, 0b10011), (5, 0b100101)]:
        q = 1 << n
        mul = [[multiply(a, b, modulus, n) for b in range(q)] for a in range(q)]
        sq = [mul[a][a] for a in range(q)]
        for a in range(1, q):
            assert mul[a].count(1) == 1
            checks['field_inverse'] += 1

        def f(x, y, z, c):
            return sq[mul[x][y] ^ mul[x][z] ^ mul[y][z]] ^ mul[mul[x][y]][z] ^ c

        # All coordinate subspaces of dimensions 2 and 3 in this fixed basis.
        for d in (2, 3):
            for indices in itertools.combinations(range(n), d):
                basis = [1 << i for i in indices]
                values = span(basis)
                m = 1 << d
                triples = [(values[k & (m-1)], values[(k >> d) & (m-1)],
                            values[k >> (2*d)]) for k in range(1 << (3*d))]
                base_table = [f(x, y, z, 0) for x, y, z in triples]
                base_anf = anf(base_table)
                assert anf(base_anf) == base_table
                checks['anf_roundtrip'] += 1
                fingerprint = hashlib.sha256(bytes(base_anf[1:])).hexdigest()
                degrees = [max((mask.bit_count() for mask, c in enumerate(base_anf)
                                if (c >> bit) & 1), default=-1) for bit in range(n)]
                for c in range(1, q):
                    table = [v ^ c for v in base_table]
                    coefficients = anf(table)
                    assert coefficients[1:] == base_anf[1:]
                    assert coefficients[0] == base_anf[0] ^ c
                    checks['constants_only'] += 1
                    transported_zeros = 0
                    for (x, y, z), v in zip(triples, table):
                        conjugate = f(sq[x], sq[y], sq[z], sq[c])
                        assert conjugate == sq[v]
                        transported_zeros += conjugate == 0
                        checks['frobenius_covariance'] += 1
                        assert f(y, x, z, c) == v == f(y, z, x, c)
                        checks['permutation_symmetry'] += 1
                    zeros = table.count(0)
                    assert transported_zeros == zeros
                    records.append(dict(n=n, modulus=modulus, dimension=d,
                                        basis=basis, a6=c, assignments=len(table),
                                        polynomial_zeros=zeros,
                                        coordinate_degrees=degrees,
                                        nonconstant_fingerprint=fingerprint))
    summary = []
    for n in (3, 4, 5):
        for d in (2, 3):
            rows = [r for r in records if r['n'] == n and r['dimension'] == d]
            groups = {}
            for r in rows:
                groups.setdefault(tuple(r['basis']), []).append(r['polynomial_zeros'])
            summary.append(dict(n=n, dimension=d, systems=len(rows),
                                subspaces=len(groups),
                                min_zeros=min(r['polynomial_zeros'] for r in rows),
                                max_zeros=max(r['polynomial_zeros'] for r in rows),
                                subspaces_with_varying_zero_counts=sum(len(set(v)) > 1 for v in groups.values()),
                                max_boolean_degree=max(max(r['coordinate_degrees']) for r in rows)))
    return dict(scope='Toy coefficient audit; no isogeny-neighbor construction or DLP solver',
                caveats=['Polynomial zeros are not verified rational-point relations.',
                         'Boolean algebraic degree is not Groebner solving degree or degree of regularity.',
                         'Coefficient sweep does not hold point count or isogeny class fixed.',
                         'Frobenius control transports the subspace along with the coefficient.',
                         'Only coordinate subspaces in the specified polynomial basis are sampled.'],
                checks=checks, summary=summary, records=records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results.json'))
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'summary')}, indent=2))
