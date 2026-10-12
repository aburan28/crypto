#!/usr/bin/env python3
"""Exact-target structural preflight; no challenge logarithm is computed."""
from __future__ import annotations

import hashlib
import argparse
import json
import random
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
DEPENDENCIES = Path(__file__).resolve().parents[1] / 'ecc2k130_relations'
sys.path.insert(0, str(DEPENDENCIES))
from fastfield import FastGF2m, IRR131
from relations import Koblitz, CHALLENGE_ELL, CHALLENGE_PX, CHALLENGE_PY, CHALLENGE_QX, CHALLENGE_QY

OUT = Path(__file__).resolve().parent / 'twist_torsion_results.json'
ELL = 263


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class CountedField(FastGF2m):
    def __init__(self):
        self.counts = dict(mul=0, sqr=0, inv=0)
        super().__init__(131, IRR131)

    def mul(self, a, b):
        self.counts['mul'] += 1
        return super().mul(a, b)

    def sqr(self, a):
        self.counts['sqr'] += 1
        return super().sqr(a)

    def inv(self, a):
        self.counts['inv'] += 1
        assert a
        u, v, g, h = a, self.irr, 1, 0
        while u != 1:
            shift = u.bit_length() - v.bit_length()
            if shift < 0:
                u, v, g, h = v, u, h, g
                shift = -shift
            u ^= v << shift
            g ^= h << shift
        return self.reduce(g)


class CountedCurve(Koblitz):
    def __init__(self, field, a=0, b=1):
        super().__init__(field, a, b)
        self.add_calls = 0

    def add(self, p, q):
        self.add_calls += 1
        return super().add(p, q)


def point_json(p):
    return None if p is None else [hex(p[0]), hex(p[1])]


def recurrence():
    a, b = 1, 0
    for _ in range(131):
        a, b = -2*b, a-b
    t = 2*a-b
    return a, b, t, (1 << 131)+1-t, (1 << 131)+1+t


def run(seed):
    start = time.monotonic()
    rng = random.Random(seed)
    a, b, trace, order, twist_order = recurrence()
    assert order == 4*CHALLENGE_ELL
    assert b == ELL*146505763881528721
    assert a % ELL == ELL-1 and b % ELL == 0
    h = twist_order//ELL**2
    assert twist_order % ELL**2 == 0 and h % ELL != 0
    f = CountedField()
    twist = CountedCurve(f, a=1)
    original = CountedCurve(f, a=0)
    inversion_checked = 0
    for _ in range(12):
        value = rng.randrange(1, 1 << 131)
        assert f.inv(value) == FastGF2m.inv(f, value)
        assert f.mul(value, f.inv(value)) == 1
        inversion_checked += 1
    public_p = (CHALLENGE_PX, CHALLENGE_PY)
    public_q = (CHALLENGE_QX, CHALLENGE_QY)
    assert original.on_curve(public_p) and original.on_curve(public_q)
    assert not twist.on_curve(public_p)
    assert original.mul(public_p, CHALLENGE_ELL) is None
    attempts = 0

    def project():
        nonlocal attempts
        while True:
            assert time.monotonic()-start < 180, 'contract time limit'
            attempts += 1
            pts = twist.points_over(rng.randrange(1, 1 << 131))
            if not pts:
                continue
            point = pts[rng.randrange(2)]
            assert twist.on_curve(point)
            torsion = twist.mul(point, h)
            if torsion is not None:
                assert twist.on_curve(torsion)
                assert twist.mul(torsion, ELL) is None
                return torsion

    u = project()
    multiples_u = [None]
    for _ in range(1, ELL):
        multiples_u.append(twist.add(multiples_u[-1], u))
    assert len(set(multiples_u)) == ELL
    lookup_u = {point: i for i, point in enumerate(multiples_u)}
    v = project()
    while v in lookup_u:
        v = project()
    assert v not in lookup_u

    def coordinates(point):
        current = point
        minus_v = twist.neg(v)
        for bv in range(ELL):
            if current in lookup_u:
                return lookup_u[current], bv
            current = twist.add(current, minus_v)
        raise AssertionError('Frobenius image outside rational torsion span')

    mu = coordinates(twist.frobenius(u))
    mv = coordinates(twist.frobenius(v))
    # Matrix columns are images of u,v; twist trace is +1 and norm is 2.
    assert (mu[0]+mv[1]) % ELL == 1
    assert (mu[0]*mv[1]-mv[0]*mu[1]) % ELL == 2
    for p in (u, v):
        tau_p = twist.frobenius(p)
        assert twist.add(twist.frobenius(tau_p), twist.mul(p, 2)) == tau_p

    lines = [(1, slope) for slope in range(ELL)] + [(0, 1)]

    def normalize(x, y):
        x, y = x % ELL, y % ELL
        assert x or y
        return (1, y*pow(x, -1, ELL) % ELL) if x else (0, 1)

    def next_line(line):
        x, y = line
        return normalize(mu[0]*x+mv[0]*y, mu[1]*x+mv[1]*y)

    orbits, seen = [], set()
    for line in lines:
        if line in seen:
            continue
        orbit, current = [], line
        while current not in orbit:
            assert current not in seen
            orbit.append(current)
            seen.add(current)
            current = next_line(current)
        assert current == line
        orbits.append(orbit)
    assert sorted(map(len, orbits)) == [1, 1, 131, 131]

    roots, codomain_b, generators = {}, {}, {}
    multiples_v = [None]
    for _ in range(1, ELL):
        multiples_v.append(twist.add(multiples_v[-1], v))
    for line in lines:
        assert time.monotonic()-start < 180, 'contract time limit'
        gen = v if line == (0, 1) else twist.add(u, multiples_v[line[1]])
        generators[line] = gen
        point, xs, total = gen, set(), 0
        for _ in range((ELL-1)//2):
            assert point is not None and twist.on_curve(point)
            assert point[0] not in xs
            xs.add(point[0])
            total ^= point[0]
            point = twist.add(point, gen)
        assert point == twist.neg(twist.mul(gen, (ELL-1)//2))
        roots[line] = frozenset(xs)
        codomain_b[line] = 1 ^ total ^ f.sqr(total)
        assert codomain_b[line] != 0
    assert len(set(roots.values())) == ELL+1
    assert len(set().union(*roots.values())) == (ELL**2-1)//2
    for line in lines:
        following = next_line(line)
        assert frozenset(f.sqr(x) for x in roots[line]) == roots[following]
        assert f.sqr(codomain_b[line]) == codomain_b[following]
    loops = [line for line in lines if codomain_b[line] == 1]
    assert len(loops) == 2
    assert all(len(orbit) != 1 or orbit[0] in loops for orbit in orbits)
    assert len(set(codomain_b.values())) == ELL
    assert len(set(codomain_b.values())-{1}) == ELL-1

    def x_map(x, kernel_xs):
        kernel_xs = sorted(kernel_xs)
        denominators = [x ^ xq for xq in kernel_xs]
        assert all(denominators), 'evaluation in kernel'
        inverses = f.batch_inv(denominators)
        accum = 0
        for xq, inv in zip(kernel_xs, inverses):
            accum ^= f.mul(xq, inv)
        return x ^ accum ^ f.sqr(accum)

    map_checks = []
    for orbit in orbits:
        line = orbit[0]
        xs = roots[line]
        destination = CountedCurve(f, a=0, b=codomain_b[line])
        image_px = x_map(public_p[0], xs)
        lifts = destination.points_over(image_px)
        assert len(lifts) == 2
        image_p = lifts[0]
        assert destination.mul(image_p, CHALLENGE_ELL) is None
        image_qx = x_map(public_q[0], xs)
        image_q_lifts = destination.points_over(image_qx)
        assert len(image_q_lifts) == 2
        assert destination.mul(image_q_lifts[0], CHALLENGE_ELL) is None
        for scalar in (1, 2, 3, 17, 263):
            source = original.mul(public_p, scalar)
            expected = destination.mul(image_p, scalar)
            assert expected is not None
            assert x_map(source[0], xs) == expected[0]
        # Original-curve 263-torsion is not rational, although x is rational.
        assert not original.has_point(next(iter(xs)))
        map_checks.append({
            'line': list(line), 'orbit_length': len(orbit),
            'kernel_generator_on_twist': point_json(generators[line]),
            'kernel_abscissae': [hex(x) for x in sorted(xs)],
            'codomain_a2': 0, 'codomain_b': hex(codomain_b[line]),
            'codomain_j': hex(f.inv(codomain_b[line])),
            'image_P_x': hex(image_px), 'image_Q_x': hex(image_qx),
            'planted_scalar_checks': [1, 2, 3, 17, 263],
            'images_in_prime_order_subgroup': True,
            'sign_scope': 'x map only; full y-map and oriented relations are not constructed',
            'destination_group_add_calls': destination.add_calls,
        })
    return {
        'seed': seed, 'elapsed_seconds': time.monotonic()-start,
        'integer_parameters': {'q': str(1 << 131), 'trace': str(trace),
             'order': str(order), 'twist_order': str(twist_order),
             'twist_263_valuation': 2, 'torsion_projection_H': str(h),
             'tau131_A': str(a), 'tau131_B': str(b),
             'q_frobenius_on_original_263_torsion': -1},
        'field_modulus_hex': hex(IRR131), 'projection_attempts': attempts,
        'basis_on_twist': [point_json(u), point_json(v)],
        'frobenius_matrix_columns_mod263': [list(mu), list(mv)],
        'frobenius_kernel_line_orbit_lengths': sorted(map(len, orbits)),
        'kernel_lines': len(lines), 'distinct_kernel_x_coordinates': (ELL**2-1)//2,
        'horizontal_loop_count': len(loops), 'distinct_descending_j_count': ELL-1,
        'distinct_j_count_including_original': ELL,
        'quotient_formula': 'b_prime = b + t + t^2; t=sum of half-kernel x coordinates',
        'frobenius_covariance_verified_all_lines': True,
        'independent_inversion_checks': inversion_checked,
        'representative_x_maps': map_checks,
        'operation_counts': {'field': f.counts, 'twist_group_add_calls': twist.add_calls,
                             'original_group_add_calls': original.add_calls},
        'counter_scope': 'diagnostic native counters; inversions include variable-time polynomial Euclid; no calibrated total-work unit',
        'PDP_cost': None, 'full_ECDLP_cost': None, 'end_to_end_speedup': None,
        'status': 'PASS_EXACT_TARGET_STRUCTURAL_PREFLIGHT',
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=OUT)
    args = parser.parse_args()
    assert not args.out.exists(), 'preserve prior evidence; use a fresh --out'
    payload = {
        'schema': 'ecc2k130_twist_torsion_preflight_v1',
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'source_sha256': sha(__file__),
        'contract_sha256': sha(Path(__file__).with_name('twist_torsion_contract.md')),
        'dependency_sha256': {name: sha(DEPENDENCIES/name) for name in ('fastfield.py', 'relations.py')},
        'scope': 'Exact public-target torsion and explicit x-coordinate isogeny checks, not a DLP solution or PDP improvement.',
        'runs': [],
    }
    for seed in (20260924, 20260925):
        result = run(seed)
        payload['runs'].append(result)
        print(json.dumps({k: result[k] for k in ('seed', 'status', 'elapsed_seconds', 'frobenius_kernel_line_orbit_lengths', 'distinct_descending_j_count')}), flush=True)
    args.out.write_text(json.dumps(payload, indent=2)+'\n')
    print(str(args.out))


if __name__ == '__main__':
    main()
