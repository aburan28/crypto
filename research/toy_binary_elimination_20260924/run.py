"""Exhaustive six-variable Boolean ideal audit; fixed tiny fields only."""
import argparse
import itertools
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       'toy_binary_formula_audit_20260924'))
from audit import anf, multiply, span


def terms(poly):
    while poly:
        bit = poly & -poly
        yield bit.bit_length() - 1
        poly ^= bit


def matrix_profile(coefficients, zeros):
    """Literal monomial multiples, not an F4/F5 implementation."""
    generators = []
    for i in range(5):
        poly = sum(1 << m for m, c in enumerate(coefficients) if c & (1 << i))
        if poly:
            generators.append((poly, max(m.bit_count() for m in terms(poly))))
    evaluations = [sum(1 << m for m in range(64) if m & z == m) for z in zeros]
    buckets = [[] for _ in range(10)]
    for poly, degree in generators:
        for multiplier in range(64):
            product = 0
            for m in terms(poly):
                product ^= 1 << (m | multiplier)
            buckets[degree + multiplier.bit_count()].append(product)
    pivots = {}
    xors = rows = 0
    saturation = None
    profile = []
    contradiction = None
    high_columns = sum(1 << m for m in range(64) if m.bit_count() > 1)
    for degree, bucket in enumerate(buckets):
        for row in bucket:
            assert all((row & e).bit_count() % 2 == 0 for e in evaluations)
            rows += 1
            while row:
                lead = row.bit_length() - 1
                if lead not in pivots:
                    pivots[lead] = row
                    break
                row ^= pivots[lead]
                xors += 1
        rank = len(pivots)
        assert rank <= 64 - len(zeros)
        if rank == 64 - len(zeros) and saturation is None:
            saturation = degree
        projected = {}
        for row in pivots.values():
            v = row & high_columns
            while v:
                lead = v.bit_length() - 1
                if lead not in projected:
                    projected[lead] = v
                    break
                v ^= projected[lead]
        remainder = 1
        while remainder and remainder.bit_length()-1 in pivots:
            remainder ^= pivots[remainder.bit_length()-1]
        if not remainder and contradiction is None:
            contradiction = degree
        profile.append(dict(degree=degree, rank=rank, rows=rows, row_xors=xors,
                            linear_consequence_dimension=rank-len(projected)))
    assert len(pivots) == 64 - len(zeros)
    assert (contradiction is not None) == (not zeros)
    final_linear = profile[-1]['linear_consequence_dimension']
    linear_complete = next(p['degree'] for p in profile
                           if p['linear_consequence_dimension'] == final_linear)
    return dict(saturation_degree=saturation, stages=profile,
                contradiction_degree=contradiction,
                linear_completion_degree=linear_complete)


def run():
    curves, records = [], []
    pair_checks = 0
    for n, modulus in [(3, 11), (4, 19), (5, 37)]:
        q = 1 << n
        mul = [[multiply(a, b, modulus, n) for b in range(q)] for a in range(q)]
        sq = [mul[x][x] for x in range(q)]
        inv = [0] + [mul[x].index(1) for x in range(1, q)]

        def f(x, y, z, c):
            return sq[mul[x][y] ^ mul[x][z] ^ mul[y][z]] ^ mul[mul[x][y]][z] ^ c

        def add(p, r):
            if p is None:
                return r
            if r is None:
                return p
            x, y = p
            u, v = r
            if x == u:
                if y != v or x == 0:
                    return None
                slope = x ^ mul[y][inv[x]]
            else:
                slope = mul[y ^ v][inv[x ^ u]]
            w = sq[slope] ^ slope ^ x ^ u
            return (w, mul[slope][x ^ w] ^ w ^ y)

        def inspect(c, basis):
            values = span(basis)
            triples = [(values[k & 3], values[(k >> 2) & 3], values[k >> 4])
                       for k in range(64)]
            table = [f(*t, c) for t in triples]
            zeros = [k for k, v in enumerate(table) if not v]
            lifted, nonlifted = [], []
            for k in zeros:
                xs = triples[k]
                valid = any(add(add(p, r), s) is None
                            for p, r, s in itertools.product(*(fibers[x] for x in xs)))
                if valid:
                    lifted.append(k)
                else:
                    nonlifted.append(k)
                if all(fibers[x] for x in xs):
                    assert valid
            coefficients = anf(table)
            assert anf(coefficients) == table
            return dict(polynomial_zeros=zeros, rational_lifted_triples=lifted,
                        nonlifted_triples=nonlifted,
                        elimination=matrix_profile(coefficients, zeros))

        for c in range(1, q):
            points = [(x, y) for x in range(q) for y in range(q)
                      if sq[y] ^ mul[x][y] == mul[sq[x]][x] ^ c]
            points_set = set(points) | {None}
            fibers = [[p for p in points if p[0] == x] for x in range(q)]
            for p in points_set:
                assert add(p, None) == p == add(None, p)
                if p is not None:
                    assert add(p, (p[0], p[0] ^ p[1])) is None
                for r in points_set:
                    s = add(p, r)
                    assert s in points_set and s == add(r, p)
                    if p is not None and r is not None and s is not None:
                        assert f(p[0], r[0], s[0], c) == 0
                    pair_checks += 1
            orbit = {c}
            v = sq[c]
            while v != c:
                orbit.add(v)
                v = sq[v]
            curves.append(dict(n=n, modulus=modulus, a2=0, a6=c,
                               point_count=len(points_set), trace=q+1-len(points_set),
                               points=points, frobenius_orbit=sorted(orbit)))
            for label, basis in [('fixed', [1, 2]), ('squared', [sq[1], sq[2]])]:
                records.append(dict(n=n, a6=c, basis=basis, basis_kind=label,
                                    point_count=len(points_set), orbit_rep=min(orbit),
                                    **inspect(c, basis)))
        lookup = {(r['a6'], r['basis_kind']): r for r in records if r['n'] == n}
        for c in range(1, q):
            a, b = lookup[c, 'fixed'], lookup[sq[c], 'squared']
            assert a['point_count'] == b['point_count']
            assert a['polynomial_zeros'] == b['polynomial_zeros']
            assert a['rational_lifted_triples'] == b['rational_lifted_triples']

    groups = []
    for n in (3, 4, 5):
        rows = [r for r in records if r['n'] == n and r['basis_kind'] == 'fixed']
        for count in sorted({r['point_count'] for r in rows}):
            group = [r for r in rows if r['point_count'] == count]
            groups.append(dict(n=n, point_count=count, coefficients=[r['a6'] for r in group],
                               frobenius_orbits=len({r['orbit_rep'] for r in group}),
                               zero_counts=sorted({len(r['polynomial_zeros']) for r in group}),
                               lift_counts=sorted({len(r['rational_lifted_triples']) for r in group}),
                               saturation_degrees=sorted({r['elimination']['saturation_degree'] for r in group})))
    matched = []
    for a, b in itertools.combinations([r for r in records if r['basis_kind']=='fixed'], 2):
        if (a['n'], a['point_count'], len(a['polynomial_zeros'])) != (
                b['n'], b['point_count'], len(b['polynomial_zeros'])):
            continue
        if a['orbit_rep'] == b['orbit_rep']:
            continue
        matched.append(dict(n=a['n'], point_count=a['point_count'], a6_pair=[a['a6'], b['a6']],
                            polynomial_zeros=len(a['polynomial_zeros']),
                            saturation_pair=[a['elimination']['saturation_degree'], b['elimination']['saturation_degree']],
                            final_xor_pair=[a['elimination']['stages'][-1]['row_xors'], b['elimination']['stages'][-1]['row_xors']]))
    return dict(protocol='EXP-TOY-BINARY-ELIMINATION-001', curves=curves, records=records,
                summary=dict(curves=len(curves), systems=len(records), point_pair_checks=pair_checks,
                             frobenius_pairs=len(curves), groups=groups,
                             nonconjugate_equal_zero_count_pairs=matched),
                end_to_end_operations=None, S=None, rho_ratio=None, speedup=None,
                degree_of_regularity=None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results.json'))
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n')
    print(json.dumps(result['summary'], indent=2))
