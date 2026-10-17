#!/usr/bin/env python3
"""Apply the preregistered fine-grid arm rule to independently replayed data."""

import argparse
import hashlib
import json
import math
import statistics
from pathlib import Path

K_VALUES = (2752, 2816, 2880, 2944, 3072)
TARGET_COUNT = 128
BATCH_TARGETS_SHA256 = '72486b7d708023838db1104a4a5e75a04b915d6f6d4013057c934fb31000aa87'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def select(root):
    replay_path = root / 'batch_full_replay.json'
    try:
        replay = json.loads(replay_path.read_text())
    except (OSError, json.JSONDecodeError):
        replay = {}
    complete = (replay.get('schema') == 'n83.intermediate-k-batch-replay/v1'
                and replay.get('status') == 'PASS_FULL_INDEPENDENT_REPLAY'
                and replay.get('targets_sha256') == BATCH_TARGETS_SHA256
                and len(replay.get('rows', [])) == len(K_VALUES)
                and [row.get('k') for row in replay['rows']] == list(K_VALUES))
    records = []
    if complete:
        for k, proof in zip(K_VALUES, replay['rows']):
            path = root / f'batch_k{k}' / 'batch_result.json'
            try:
                row = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError):
                complete = False
                break
            times = row.get('per_target_online_ms', [])
            if (proof.get('status') != 'PASS_FULL_INDEPENDENT_REPLAY'
                    or proof.get('target_count') != TARGET_COUNT
                    or proof.get('batch_result_sha256') != digest(path)
                    or row.get('k') != k
                    or row.get('status') !=
                    'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY'
                    or row.get('target_count') != TARGET_COUNT
                    or len(times) != TARGET_COUNT
                    or any(not isinstance(t, (int, float)) or
                           not math.isfinite(t) or t <= 0 for t in times)):
                complete = False
                break
            geometric_mean = math.exp(statistics.mean(math.log(t) for t in times))
            if not math.isclose(geometric_mean,
                                row['online_geometric_mean_ms'],
                                rel_tol=1e-12):
                complete = False
                break
            records.append({'k': k, 'candidate_id': proof['candidate_id'],
                            'geometric_mean_online_ms': geometric_mean,
                            'arithmetic_mean_online_ms': statistics.mean(times),
                            'batch_result_sha256': digest(path),
                            'raw_summary_sha256': proof['summary_sha256'],
                            'raw_base_dump_sha256': proof['base_dump_sha256'],
                            'raw_target_results_sha256':
                                proof['target_results_sha256']})
    selected = min(records, key=lambda row: (
        row['geometric_mean_online_ms'], row['k']))['k'] if complete else None
    return {'schema': 'n83.root11-finegrid-local-selection/v1',
            'status': 'SELECTED_BY_COMPLETE_REPLAY' if complete else
                      'BOUNDED_UNKNOWN_INCOMPLETE_REPLAY',
            'rule': 'minimum_geometric_mean_128_online_ms_then_lower_k',
            'selected_k': selected,
            'ordered_k_values': list(K_VALUES),
            'target_count': TARGET_COUNT,
            'batch_replay_sha256': digest(replay_path) if replay_path.exists()
                                   else None,
            'batch_targets_sha256': replay.get('targets_sha256'),
            'rows': records if complete else []}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('batch_root', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = select(args.batch_root)
    blob = json.dumps(result, sort_keys=True, indent=2) + '\n'
    if args.output.exists() and args.output.read_text() != blob:
        raise ValueError('selection already exists with different bytes')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(blob)
    print(json.dumps({'status': result['status'],
                      'selected_k': result['selected_k']}))


if __name__ == '__main__':
    main()
