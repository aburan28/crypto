# Pair-enumeration decomposition, the cheapest oracle this repository has
# measured for three-summand index calculus (RESEARCH_ECC2K130_RR_SOLVER_PANEL.md
# §7).  Used as the CPU reference for the G7e CUDA instrument and as a
# small-field end-to-end check that does not need a SAT solver.
#
# No type hints, camelCase identifiers, no itertools (project convention).

import argparse
import json
import math
import random
import time

import indexcalc_e2e as engine


RHO_LOG2 = 60.809
R_131 = 680564733841876926932320129493409985129
ORDER_131 = 4 * R_131
# Itoh–Tsujii inv131 is eight field products; affine addition uses two more.
PRODUCTS_PER_AFFINE_ADD = 10
# One pair search does P_i + P_j and then R - S.
PRODUCTS_PER_PAIR = 2 * PRODUCTS_PER_AFFINE_ADD


def log2(value):
    if value <= 0:
        return float('-inf')
    return math.log(value, 2)


def combinations(n, k):
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def expectedYield(baseSize, summands, groupOrder):
    return combinations(baseSize, summands) / float(groupOrder)


def streamingProducts(baseSize, relations, groupOrder, summands=3):
    """Field products to collect `relations` by streaming pair enumeration.

    Each trial scans every unordered pair.  A random target in a group of
    `groupOrder` that contains every sum of `summands` base points hits with
    probability C(|F|, m) / groupOrder, so the expected trial count is the
    reciprocal.  This is the product-law cost of the oracle, not a GPU
    constant."""
    hits = expectedYield(baseSize, summands, groupOrder)
    if hits <= 0:
        return float('inf')
    pairs = combinations(baseSize, 2)
    return relations * (1.0 / hits) * pairs * PRODUCTS_PER_PAIR


def sScore(operations, subgroupOrder):
    return operations / math.sqrt(subgroupOrder)


def pairDecompose(curve, points, target, summands=3):
    """Return one unordered tuple of `summands` points adding to target, or None.

    `points` is a list.  Distinctness is by object identity of the chosen
    list entries; the caller must not store a point and its negative as the
    same Python object.
    """
    if target is None or summands != 3:
        return None
    n = len(points)
    for i in range(n):
        pi = points[i]
        for j in range(i + 1, n):
            pj = points[j]
            if pi == curve.neg(pj):
                continue
            partial = curve.add(pi, pj)
            if partial is None:
                continue
            t = curve.add(target, curve.neg(partial))
            if t is None:
                continue
            for k in range(j + 1, n):
                if t == points[k]:
                    if t == curve.neg(pi) or t == curve.neg(pj):
                        continue
                    return (pi, pj, t)
    return None


def pairDecomposeLookup(curve, lookup, target, summands=3):
    """Meet-in-the-middle form: remainder membership in `lookup`."""
    if target is None or summands != 3:
        return None
    points = list(lookup)
    n = len(points)
    for i in range(n):
        pi = points[i]
        for j in range(i + 1, n):
            pj = points[j]
            if pi == curve.neg(pj):
                continue
            partial = curve.add(pi, pj)
            if partial is None:
                continue
            t = curve.add(target, curve.neg(partial))
            if t is None or t not in lookup:
                continue
            if t == pi or t == pj or t == curve.neg(pi) or t == curve.neg(pj):
                continue
            return (pi, pj, t)
    return None


def relationRow(lookup, reps, ell, triple):
    row = [0] * len(reps)
    for point in triple:
        column, coeff = lookup[point]
        row[column] = (row[column] + coeff) % ell
    return row


def recoverLog(m, weight, seed, attempts=256):
    """Planted-log recovery on a tiny K_0 by pair enumeration.

    Degree 131 is refused: the pair table does not fit a correctness test
    that is meant to finish in unit-test time.
    """
    if m > 11:
        raise ValueError('pair-enumeration DLP recovery is a toy-field check')
    meter = engine.Ledger()
    context = engine.setup(m, weight, 3, meter)
    onb, curve, ell, generator, eigen, reps, lookup, _, _ = context
    rng = random.Random(seed)
    expected = rng.randrange(1, ell)
    target = curve.mul(generator, expected)
    matrix = engine.RelationMatrix(len(reps), ell, meter)
    records = []
    t0 = time.perf_counter_ns()
    # Factor-base logarithms from known multiples of G, then one mixed descent
    # onto Q.  The secret is not used as a right-hand side.
    for i in range(attempts):
        scalar = rng.randrange(1, ell)
        probe = curve.mul(generator, scalar)
        triple = pairDecomposeLookup(curve, lookup, probe)
        detail = {'index': i, 'status': 'miss' if triple is None else 'hit'}
        if triple is not None:
            row = relationRow(lookup, reps, ell, triple)
            detail['independent'] = matrix.push(row, scalar)
            got = None
            for p in triple:
                got = curve.add(got, p)
            assert got == probe
        records.append(detail)
        if matrix.solution() is not None:
            break
    elapsed = time.perf_counter_ns() - t0
    logs = matrix.solution()
    report = {
        'degree': m, 'weight': weight, 'seed': seed, 'attempts': attempts,
        'subgroup_order': str(ell), 'factor_base_points': len(lookup),
        'orbit_columns': len(reps), 'records': records,
        'elapsed_ns': elapsed, 'status': 'insufficient_relations',
        'scalar_recovered': None, 'verified': False,
        'oracle_class': 'enumerate_pairs',
    }
    if logs is None:
        return report
    for log, point in zip(logs, reps):
        assert curve.mul(generator, log) == point
    descentRng = random.Random(seed ^ 0xDE5C)
    for _ in range(64):
        shift = descentRng.randrange(1, ell)
        residual = curve.add(target, curve.mul(generator, shift))
        triple = pairDecomposeLookup(curve, lookup, residual)
        if triple is None:
            continue
        row = relationRow(lookup, reps, ell, triple)
        recovered = (sum(c * log for c, log in zip(row, logs)) - shift) % ell
        if curve.mul(generator, recovered) != target:
            continue
        report.update({
            'status': 'complete',
            'scalar_recovered': str(recovered),
            'verified': recovered == expected,
            'expected_scalar': str(expected),
        })
        assert recovered == expected
        return report
    report['status'] = 'descent_budget'
    return report


def projectAttack(baseSize, relations, subgroupOrder, curveOrder, rhoLog2=RHO_LOG2):
    products = streamingProducts(baseSize, relations, subgroupOrder)
    s = sScore(products, subgroupOrder)
    sRho = sScore(2 ** rhoLog2, subgroupOrder)
    return {
        'base_size': baseSize,
        'relations_charged': relations,
        'group_order': str(curveOrder),
        'subgroup_order': str(subgroupOrder),
        'expected_yield_per_target': expectedYield(baseSize, 3, subgroupOrder),
        'streaming_field_products': products,
        'log2_streaming_field_products': log2(products) if products < float('inf') else None,
        'S': s if products < float('inf') else None,
        'S_rho': sRho,
        'ratio_to_rho': (products / (2 ** rhoLog2)) if products < float('inf') else None,
        'log2_ratio_to_rho': (log2(products) - rhoLog2) if products < float('inf') else None,
        'products_per_affine_add': PRODUCTS_PER_AFFINE_ADD,
        'products_per_pair': PRODUCTS_PER_PAIR,
        'class': 'engineering',
        'note': 'GPU pair-add throughput changes wall-clock, not S; the product law is the floor',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--m', type=int, default=9)
    parser.add_argument('--weight', type=int, default=2)
    parser.add_argument('--seed', type=int, default=20260918)
    parser.add_argument('--attempts', type=int, default=256)
    parser.add_argument('--json', action='store_true')
    args = parser.parse_args()
    report = recoverLog(args.m, args.weight, args.seed, args.attempts)
    if args.json:
        print(json.dumps(report, sort_keys=True))
    else:
        print('m=%d weight=%d status=%s recovered=%s verified=%s points=%d orbits=%d'
              % (report['degree'], args.weight, report['status'],
                 report['scalar_recovered'], report['verified'],
                 report['factor_base_points'], report['orbit_columns']))
    return 0 if report['status'] == 'complete' and report['verified'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
