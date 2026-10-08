#!/usr/bin/env python3
"""Independently replay compact-orbit rank traces in Python group arithmetic."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[4]
REPLAY = ROOT / 'research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924'
sys.path.insert(0, str(REPLAY))
from independent_replay import Curve, Field  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point(value):
    return None if value is None else tuple(value)


def read_jsonl(path: Path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


class Echelon:
    def __init__(self, columns: int, modulus: int):
        self.columns = columns
        self.modulus = modulus
        self.pivots = {}

    def insert(self, original: list[int]) -> bool:
        row = list(original)
        r = self.modulus
        for column in range(self.columns):
            if row[column] == 0:
                continue
            pivot = self.pivots.get(column)
            if pivot is not None:
                scale = row[column]
                for index in range(column, self.columns + 1):
                    row[index] = (row[index] - scale * pivot[index]) % r
            else:
                inverse = pow(row[column], -1, r)
                for index in range(column, self.columns + 1):
                    row[index] = row[index] * inverse % r
                self.pivots[column] = row
                return True
        return False

    def solve(self) -> list[int]:
        assert len(self.pivots) == self.columns
        logs = [0] * self.columns
        r = self.modulus
        for column in range(self.columns - 1, -1, -1):
            row = self.pivots[column]
            logs[column] = (row[-1] - sum(row[j] * logs[j]
                                           for j in range(column + 1, self.columns))) % r
        return logs


def verify(trace_path: Path, base_path: Path, summary_path: Path) -> dict:
    trace = read_jsonl(trace_path)
    assert len(trace) >= 3
    header, solution = trace[0], trace[-1]
    attempts = trace[1:-1]
    base = read_jsonl(base_path)
    assert len(base) == 1
    base = base[0]
    summaries = read_jsonl(summary_path)
    assert len(summaries) == 1
    summary = summaries[0]
    assert header['kind'] == 'compact_orbit_rank_header'
    assert header['schema_version'] == 1
    assert solution['kind'] == 'compact_orbit_rank_solution'
    assert base['kind'] == 'point_defined_factor_base'
    assert summary['kind'] == 'compact_orbit_dlp_summary'
    n, a, r, k = (int(header[name]) for name in
                  ('n', 'a', 'subgroup_order', 'orbit_columns'))
    for record in (base, summary):
        assert record['n'] == n and record['a'] == a
        assert record['base_hash'] == header['base_hash']
        assert record['orbit_columns'] == k
    assert base['subgroup_order'] == r
    assert header['factor_base_points'] == base['factor_base_points']
    field = Field(n, base['field_modulus_low_terms'])
    curve = Curve(field, a)
    generator = point(header['generator'])
    assert generator is not None and curve.on_curve(generator)
    assert curve.mul(r, generator) is None
    representatives = [point(raw) for raw in base['factor_base_representatives']]
    points = [point(raw) for raw in base['factor_base_point_coordinates']]
    labels = base['factor_base_point_labels']
    assert len(representatives) == k
    assert len(points) == len(labels) == base['factor_base_points']
    for rep in representatives:
        assert rep is not None and curve.on_curve(rep) and curve.mul(r, rep) is None
    verified_points = set()
    echelon = Echelon(k, r)
    relations = failures = rows_without_gain = 0
    for index, record in enumerate(attempts):
        assert record['kind'] == 'compact_orbit_rank_attempt'
        assert record['attempt_index'] == index
        assert record['rank_before'] == len(echelon.pivots)
        scalar = record['scalar']
        column = record['pivotless_column']
        assert 1 <= scalar < r and 0 <= column < k
        assert column not in echelon.pivots
        target = point(record['target'])
        expected = curve.add(curve.mul(scalar, generator), curve.neg(representatives[column]))
        assert target == expected and curve.on_curve(target)
        if not record['found']:
            failures += 1
            assert record['rank_after'] == len(echelon.pivots)
            continue
        relations += 1
        indices, codes = record['point_indices'], record['x_codes']
        assert len(indices) == len(codes) == 4
        chosen = []
        expected_row = [0] * (k + 1)
        for position, point_index in enumerate(indices):
            assert 0 <= point_index < len(points)
            base_point = points[point_index]
            assert base_point is not None and codes[position] == base_point[0]
            base_column, coefficient = labels[point_index]
            assert 0 <= base_column < k and 0 <= coefficient < r
            if point_index not in verified_points:
                assert curve.on_curve(base_point)
                assert curve.mul(r, base_point) is None
                assert curve.mul(coefficient, representatives[base_column]) == base_point
                verified_points.add(point_index)
            expected_row[base_column] = (expected_row[base_column] + coefficient) % r
            chosen.append(base_point)
        expected_row[column] = (expected_row[column] + 1) % r
        expected_row[-1] = scalar
        assert record['row'] == expected_row
        left = curve.add(chosen[0], chosen[1])
        right = curve.add(chosen[2], chosen[3])
        assert curve.add(left, right) == target
        assert left is not None and right is not None
        assert {left[0], right[0]} == set(record['pinned_intermediates'])
        gained = echelon.insert(expected_row)
        assert record['gained'] == gained
        assert record['rank_after'] == len(echelon.pivots)
        if not gained:
            rows_without_gain += 1
    assert solution['attempts'] == len(attempts) == summary['rank_attempts']
    assert solution['relations'] == relations == summary['rank_relations']
    assert solution['failures'] == failures == summary['rank_failures']
    assert solution['rows_without_gain'] == rows_without_gain == summary['rank_rows_without_gain']
    assert solution['rank'] == len(echelon.pivots) == summary['rank'] == k
    logs = echelon.solve()
    assert solution['logs'] == logs
    for log, rep in zip(logs, representatives):
        assert curve.mul(log, generator) == rep
    return {
        'status': 'PASS', 'schema': 'compact-orbit-rank-independent-replay-v1',
        'n': n, 'a': a, 'subgroup_order': r, 'orbit_columns': k,
        'attempts': len(attempts), 'relations': relations,
        'failed_searches': failures, 'rows_without_gain': rows_without_gain,
        'rank': len(echelon.pivots), 'representative_logs_verified': k,
        'used_base_points_verified': len(verified_points),
        'miss_nonexistence_proved': False,
        'trace_sha256': sha(trace_path), 'base_sha256': sha(base_path),
        'summary_sha256': sha(summary_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--trace', type=Path, required=True)
    parser.add_argument('--base', type=Path, required=True)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), 'refusing to overwrite replay receipt'
    try:
        result = verify(args.trace, args.base, args.summary)
    except BaseException as exc:
        result = {'status': 'FAIL', 'error_type': type(exc).__name__,
                  'error': str(exc), 'traceback': traceback.format_exc()}
        args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
        raise
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
