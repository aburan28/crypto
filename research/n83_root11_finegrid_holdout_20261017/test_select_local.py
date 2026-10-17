#!/usr/bin/env python3
"""Exercise the preregistered choice, tie break, and incomplete-replay gate."""

import json
import tempfile
from pathlib import Path

import select_local


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + '\n')


def fixture(root, means):
    rows = []
    for k, mean in zip(select_local.K_VALUES, means):
        path = root / f'batch_k{k}' / 'batch_result.json'
        write(path, {
            'k': k, 'status': 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY',
            'target_count': select_local.TARGET_COUNT,
            'online_geometric_mean_ms': mean,
            'per_target_online_ms': [mean] * select_local.TARGET_COUNT})
        rows.append({'k': k, 'candidate_id': f'exact-{k}',
                     'status': 'PASS_FULL_INDEPENDENT_REPLAY',
                     'target_count': select_local.TARGET_COUNT,
                     'batch_result_sha256': select_local.digest(path),
                     'summary_sha256': f'summary-{k}',
                     'base_dump_sha256': f'base-{k}',
                     'target_results_sha256': f'targets-{k}'})
    write(root / 'batch_full_replay.json', {
        'schema': 'n83.intermediate-k-batch-replay/v1',
        'status': 'PASS_FULL_INDEPENDENT_REPLAY',
        'targets_sha256': select_local.BATCH_TARGETS_SHA256, 'rows': rows})


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    fixture(root, [50.0, 40.0, 30.0, 20.0, 10.0])
    result = select_local.select(root)
    assert result['status'] == 'SELECTED_BY_COMPLETE_REPLAY'
    assert result['selected_k'] == 3072
    assert len(result['rows']) == 5

    fixture(root, [50.0, 40.0, 10.0, 10.0, 30.0])
    assert select_local.select(root)['selected_k'] == 2880

    replay = json.loads((root / 'batch_full_replay.json').read_text())
    replay['rows'][0]['status'] = 'VERIFICATION_FAILED'
    write(root / 'batch_full_replay.json', replay)
    result = select_local.select(root)
    assert result['status'] == 'BOUNDED_UNKNOWN_INCOMPLETE_REPLAY'
    assert result['selected_k'] is None
    assert result['rows'] == []

    (root / 'batch_full_replay.json').unlink()
    assert select_local.select(root)['status'] == 'BOUNDED_UNKNOWN_INCOMPLETE_REPLAY'

print(json.dumps({'status': 'PASS', 'cases': 4}))
