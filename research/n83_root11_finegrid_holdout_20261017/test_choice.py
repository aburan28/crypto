#!/usr/bin/env python3
"""Check distinct-arm resolution and fixed-candidate identity gates."""

import json
import tempfile
from pathlib import Path

import arms
import select_local


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + '\n')


with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    results = root / 'results'
    original_results, original_manifest = arms.RESULTS, arms.manifest_path
    arms.RESULTS = results
    arms.manifest_path = lambda k, candidate_id: root / f'manifest_{k}.json'
    ids = {k: f'candidate-{k}' for k in arms.SCREEN_K}
    ids.update({k: candidate_id for k, candidate_id in arms.FIXED})
    try:
        def fixture(means):
            proofs = []
            for k, mean in zip(arms.SCREEN_K, means):
                write(root / f'manifest_{k}.json', {
                    'candidate_id': ids[k], 'factor_base': {
                        'actual_usable_points': 166 * k,
                        'effective_columns': k}})
                raw = {'k': k, 'target_count': 128,
                       'status': 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY',
                       'online_geometric_mean_ms': mean,
                       'per_target_online_ms': [mean] * 128}
                path = root / f'batch_k{k}' / 'batch_result.json'
                write(path, raw)
                write(results / f'batch_k{k}_result.json', raw)
                proofs.append({'k': k, 'candidate_id': ids[k],
                               'status': 'PASS_FULL_INDEPENDENT_REPLAY',
                               'target_count': 128,
                               'batch_result_sha256': arms.digest(path),
                               'summary_sha256': f'summary-{k}',
                               'base_dump_sha256': f'base-{k}',
                               'target_results_sha256': f'targets-{k}'})
            write(results / 'batch_full_replay.json', {
                'schema': 'n83.intermediate-k-batch-replay/v1',
                'status': 'PASS_FULL_INDEPENDENT_REPLAY',
                'targets_sha256': arms.BATCH_TARGETS_SHA256,
                'rows': proofs})
            # select_local expects the raw batch layout beside its replay.
            for k in arms.SCREEN_K:
                write(results / f'batch_k{k}' / 'batch_result.json',
                      json.loads((results / f'batch_k{k}_result.json').read_text()))
            choice = select_local.select(results)
            write(results / 'local_selection.json', choice)
            return choice

        fixture([40.0, 30.0, 10.0, 20.0, 25.0])
        write(root / 'manifest_2543.json', {
            'candidate_id': ids[2543], 'factor_base': {
                'actual_usable_points': 166 * 2543,
                'effective_columns': 2543}})
        chosen, choice = arms.choose()
        assert len(chosen) == 4 and choice['selected_fourth_k'] == 2880

        fixture([40.0, 10.0, 20.0, 30.0, 25.0])
        chosen, choice = arms.choose()
        assert len(chosen) == 3 and choice['selected_fourth_k'] == 2816

        ids[2816] = 'wrong-fixed-identity'
        fixture([40.0, 10.0, 20.0, 30.0, 25.0])
        try:
            arms.choose()
            raise AssertionError('mismatched fixed candidate was admitted')
        except ValueError:
            pass
    finally:
        arms.RESULTS, arms.manifest_path = original_results, original_manifest

print(json.dumps({'status': 'PASS', 'cases': 3}))
