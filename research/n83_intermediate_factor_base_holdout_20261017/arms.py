#!/usr/bin/env python3
"""Resolve the sealed intermediate-size holdout arms from published evidence."""

import hashlib
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
BATCH_TARGETS_SHA256 = 'ac14c0326532b81041d357e921d99f86ef7ff4f6cf5e475b9c3a2bcc8997ac7c'
IDS = {
    2304: 'IC1N83Ckb1a1fb382464PDP4rootRCguidedLAgaussTDdirectISO0hf7c27ecf3b60',
    2543: 'IC1N83Ckb1a1fb422138PDP4rootRCguidedLAgaussTDdirectISO0h7a345a837c55',
    2688: 'IC1N83Ckb1a1fb446208PDP4rootRCguidedLAgaussTDdirectISO0hc92369996197',
    2816: 'IC1N83Ckb1a1fb467456PDP4rootRCguidedLAgaussTDdirectISO0hf693b1b33c8b',
    3072: 'IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a',
}
FIXED_K = (2304, 2543, 3072)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest_path(k):
    directory = (HERE if k in (2304, 2688, 2816) else
                 ROOT / 'research/n83_smallk_challenge_20261017')
    return directory / f'candidate_{IDS[k]}.json'


def choose():
    choice_path = RESULTS / 'batch_selection.json'
    choice = json.loads(choice_path.read_text())
    replay_path = RESULTS / 'batch_full_replay.json'
    replay = json.loads(replay_path.read_text())
    if (choice['schema'] != 'n83.intermediate-batch-arm-selection/v1' or
            choice['status'] != 'SELECTED_BY_COMPLETE_REPLAY' or
            choice['rule'] !=
            'minimum_geometric_mean_128_online_ms_then_lower_k' or
            choice['ordered_k_values'] != list(IDS) or
            choice['target_count'] != 128 or
            choice['batch_targets_sha256'] != BATCH_TARGETS_SHA256 or
            choice['batch_replay_sha256'] != digest(replay_path) or
            replay['status'] != 'PASS_FULL_INDEPENDENT_REPLAY' or
            replay['targets_sha256'] != BATCH_TARGETS_SHA256 or
            [row['k'] for row in replay['rows']] != list(IDS)):
        raise ValueError('published batch selection or replay differs')
    if len(choice['rows']) != len(IDS):
        raise ValueError('five independently replayed batch arms required')
    for k, row, proof in zip(IDS, choice['rows'], replay['rows']):
        path = RESULTS / f'batch_k{k}_result.json'
        data = json.loads(path.read_text())
        times = data['per_target_online_ms']
        geometric_mean = math.exp(statistics.mean(math.log(t) for t in times))
        manifest = json.loads(manifest_path(k).read_text())
        if (row['k'] != k or row['candidate_id'] != IDS[k] or
                proof['candidate_id'] != IDS[k] or
                proof['status'] != 'PASS_FULL_INDEPENDENT_REPLAY' or
                proof['target_count'] != 128 or
                row['batch_result_sha256'] != digest(path) or
                proof['batch_result_sha256'] != digest(path) or
                row['raw_summary_sha256'] != proof['summary_sha256'] or
                row['raw_base_dump_sha256'] != proof['base_dump_sha256'] or
                row['raw_target_results_sha256'] !=
                proof['target_results_sha256'] or
                len(times) != 128 or any(t <= 0 for t in times) or
                not math.isclose(row['geometric_mean_online_ms'],
                                 geometric_mean, rel_tol=1e-12) or
                not math.isclose(row['arithmetic_mean_online_ms'],
                                 statistics.mean(times), rel_tol=1e-12) or
                manifest['candidate_id'] != IDS[k] or
                manifest['factor_base']['actual_usable_points'] != 166 * k or
                manifest['factor_base']['effective_columns'] != k):
            raise ValueError(f'K{k} published batch or manifest differs')
    winner = min(choice['rows'], key=lambda row: (
        row['geometric_mean_online_ms'], row['k']))['k']
    if winner != choice['selected_k']:
        raise ValueError('published winner differs from frozen rule')
    chosen_k = tuple(dict.fromkeys((*FIXED_K, winner)))
    if len(chosen_k) != 4:
        raise ValueError('complete screen must add one new fourth arm')
    arms = tuple((f'k{k}root11', k, 'root11_dynamic', IDS[k])
                 for k in chosen_k)
    return arms, {
        'schema': 'n83.intermediate-holdout-choice/v1',
        'source_sha256': SOURCE_SHA256,
        'batch_selection_sha256': digest(choice_path),
        'batch_replay_sha256': digest(replay_path),
        'ordered_candidate_ids': [arm[3] for arm in arms],
        'selected_fourth_k': winner,
    }
