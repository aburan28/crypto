#!/usr/bin/env python3
"""Reconstruct pair-table occupancy and compare sealed coverage probes."""
import argparse
import json
import random
from collections import defaultdict
from pathlib import Path


def affine_coefficients(points, p):
    first = points[0]
    second = next(point for point in points if point[0] != first[0])
    x1, y1 = first
    x2, y2 = second
    a = ((y2 * y2 - x2**3) - (y1 * y1 - x1**3)) * pow(x2 - x1, -1, p) % p
    b = (y1 * y1 - x1**3 - a * x1) % p
    assert all((y * y - x**3 - a * x - b) % p == 0 for x, y in points)
    return a, b


def add_x(left, right, p, a):
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if left == right:
        slope = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        slope = (y2 - y1) * pow(x2 - x1, -1, p) % p
    return (slope * slope - x1 - x2) % p


def table_size(fb):
    icv1 = fb['curve']['icv1']
    p = int(icv1.split(':')[1].split('-')[1])
    r = fb['curve']['r']
    points = [(int(point['x'], 16), int(point['y'], 16)) for point in fb['points']]
    a, b = affine_coefficients(points, p)
    by_col = defaultdict(list)
    for point, record in zip(points, fb['points']):
        by_col[record['col']].append(point)
    keys = set()
    for i, group in sorted(by_col.items()):
        positive = group[0]
        for j, other in sorted(by_col.items()):
            if j < i:
                continue
            for point in other:
                x = add_x(positive, point, p, a)
                if x is not None:
                    keys.add(x)
    assert len(keys) <= (r - 1) // 2
    return {'curve_p': p, 'curve_a': a, 'curve_b': b, 'subgroup_r': r,
            'unique_pair_keys': len(keys), 'uniform_hit_probability': 2 * len(keys) / (r - 1),
            'factor_base_sha256': fb['factor_base']['fb_sha256']}


def probe_counts(row):
    relations = next(phase for phase in row['phases'] if phase['name'] == 'relations')
    return relations['native']['lookups'], row['counters']['relations_found']


def cluster_interval(rows, probability, seed):
    by_target = defaultdict(lambda: [0, 0])
    for row in rows:
        key = row['workload']['target_index']
        probes, hits = probe_counts(row)
        by_target[key][0] += probes
        by_target[key][1] += hits
    clusters = list(by_target.values())
    rng = random.Random(seed)
    ratios = []
    for _ in range(10_000):
        draw = rng.choices(clusters, k=len(clusters))
        probes = sum(pair[0] for pair in draw)
        hits = sum(pair[1] for pair in draw)
        ratios.append(hits / (probes * probability))
    ratios.sort()
    return [ratios[249], ratios[9749]]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('session', type=Path)
    parser.add_argument('factor_bases', nargs='+', type=Path)
    args = parser.parse_args()
    bases = {}
    for path in args.factor_bases:
        fb = json.loads(path.read_text())
        stride = int(fb['factor_base']['params'].get('stride', 1))
        bases[stride] = table_size(fb)
    rows = [json.loads(line) for line in (args.session / 'records.jsonl').read_text().splitlines()]
    by_arm = defaultdict(list)
    for row in rows:
        if not row['warmup'] and row['arm'].startswith('b48-stride-'):
            by_arm[row['arm']].append(row)
    summary = {}
    for arm, arm_rows in sorted(by_arm.items()):
        stride = int(arm.rsplit('-', 1)[-1])
        base = bases[stride]
        assert all(row['factor_base']['fb_sha256'] == base['factor_base_sha256'] for row in arm_rows)
        probes = sum(probe_counts(row)[0] for row in arm_rows)
        hits = sum(probe_counts(row)[1] for row in arm_rows)
        probability = base['uniform_hit_probability']
        summary[arm] = {**base, 'measured_runs': len(arm_rows), 'targets': len({row['workload']['target_index'] for row in arm_rows}),
                        'lookups': probes, 'relations': hits,
                        'observed_hit_probability': hits / probes,
                        'observed_to_uniform_ratio': hits / (probes * probability),
                        'target_cluster_ratio_interval_95': cluster_interval(arm_rows, probability, 20261023 + stride),
                        'mean_cold_gae': sum(row['cost']['total_gae'] for row in arm_rows) / len(arm_rows),
                        'peak_rss_kib': max(row['time']['max_rss_kib'] for row in arm_rows),
                        'mean_independent_rank': sum(row['counters']['matrix_rank'] for row in arm_rows) / len(arm_rows)}
    print(json.dumps({'session': str(args.session), 'arms': summary,
                      'inference_limit': 'Walk probes are correlated and stop at target pin; the uniform probability is exact for independent uniform subgroup points, while the target-cluster intervals are descriptive.'}, indent=2))


if __name__ == '__main__':
    main()
