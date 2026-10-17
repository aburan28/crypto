#!/usr/bin/env python3
"""Resolve the prospective 24-point arms only after full local replay."""

import hashlib
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
BATCH_TARGETS_SHA256 = '72486b7d708023838db1104a4a5e75a04b915d6f6d4013057c934fb31000aa87'
SCREEN_K = (2752, 2816, 2880, 2944, 3072)
FIXED = (
    (2543, 'IC1N83Ckb1a1fb422138PDP4rootRCguidedLAgaussTDdirectISO0h7a345a837c55'),
    (2816, 'IC1N83Ckb1a1fb467456PDP4rootRCguidedLAgaussTDdirectISO0hf693b1b33c8b'),
    (3072, 'IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a'),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def manifest_path(k, candidate_id):
    directory = (ROOT / 'research/n83_smallk_challenge_20261017' if
                 k in (2543, 3072) else
                 ROOT / 'research/n83_intermediate_factor_base_holdout_20261017'
                 if k == 2816 else HERE)
    return directory / f'candidate_{candidate_id}.json'


def choose():
    choice_path = RESULTS / 'local_selection.json'
    replay_path = RESULTS / 'batch_full_replay.json'
    choice = json.loads(choice_path.read_text())
    replay = json.loads(replay_path.read_text())
    if (choice['schema'] != 'n83.root11-finegrid-local-selection/v1' or
            choice['status'] != 'SELECTED_BY_COMPLETE_REPLAY' or
            choice['rule'] != 'minimum_geometric_mean_128_online_ms_then_lower_k' or
            choice['ordered_k_values'] != list(SCREEN_K) or
            choice['target_count'] != 128 or
            choice['batch_targets_sha256'] != BATCH_TARGETS_SHA256 or
            choice['batch_replay_sha256'] != digest(replay_path) or
            replay['status'] != 'PASS_FULL_INDEPENDENT_REPLAY' or
            replay['targets_sha256'] != BATCH_TARGETS_SHA256 or
            [row['k'] for row in replay['rows']] != list(SCREEN_K) or
            len(choice['rows']) != len(SCREEN_K)):
        raise ValueError('fine-grid replay or selection incomplete')
    by_k = {}
    for k, row, proof in zip(SCREEN_K, choice['rows'], replay['rows']):
        path = RESULTS / f'batch_k{k}_result.json'
        data = json.loads(path.read_text())
        times = data['per_target_online_ms']
        manifest = json.loads(manifest_path(k, row['candidate_id']).read_text())
        geometric_mean = math.exp(statistics.mean(math.log(t) for t in times))
        if (row['k'] != k or proof['k'] != k or
                row['candidate_id'] != proof['candidate_id'] or
                row['batch_result_sha256'] != digest(path) or
                proof['batch_result_sha256'] != digest(path) or
                row['raw_summary_sha256'] != proof['summary_sha256'] or
                row['raw_base_dump_sha256'] != proof['base_dump_sha256'] or
                row['raw_target_results_sha256'] !=
                proof['target_results_sha256'] or
                len(times) != 128 or any(t <= 0 for t in times) or
                not math.isclose(row['geometric_mean_online_ms'],
                                 geometric_mean, rel_tol=1e-12) or
                manifest['candidate_id'] != row['candidate_id'] or
                manifest['factor_base']['actual_usable_points'] != 166 * k or
                manifest['factor_base']['effective_columns'] != k):
            raise ValueError(f'K{k} exact batch candidate differs')
        by_k[k] = row['candidate_id']
    selected = min(choice['rows'], key=lambda row: (
        row['geometric_mean_online_ms'], row['k']))['k']
    if selected != choice['selected_k']:
        raise ValueError('fourth arm differs from the frozen rule')
    ordered = tuple(dict.fromkeys((*FIXED, (selected, by_k[selected]))))
    if len(ordered) not in (3, 4):
        raise ValueError('three or four distinct exact arms required')
    arms = tuple((f'k{k}root11', k, 'root11_dynamic', candidate_id)
                 for k, candidate_id in ordered)
    return arms, {
        'schema': 'n83.root11-finegrid-holdout-choice/v1',
        'source_sha256': SOURCE_SHA256,
        'local_selection_sha256': digest(choice_path),
        'batch_replay_sha256': digest(replay_path),
        'selected_fourth_k': selected,
        'ordered_candidate_ids': [arm[3] for arm in arms],
    }
