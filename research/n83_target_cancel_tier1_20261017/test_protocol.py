#!/usr/bin/env python3
"""Exercise successful and censored 24-point cancellation selection paths."""

import json

import run
import selector


def sample():
    ids = [arm[3] for arm in run.ARMS]
    points = []
    for index in range(24):
        points.append({
            'point_index': index,
            'arms': [
                {'candidate_id': ids[0], 'status': 'PASS_INDEPENDENT_REPLAY',
                 'resource_status': 'PASS', 'one_target_online_ms': 100.0,
                 'cold_process_ms': 100_000.0},
                {'candidate_id': ids[1], 'status': 'PASS_INDEPENDENT_REPLAY',
                 'resource_status': 'PASS', 'one_target_online_ms': 50.0,
                 'cold_process_ms': 90_000.0},
            ],
            'rho_by_candidate': {
                candidate: {'status': 'PASS_INDEPENDENT_REPLAY',
                            'resource_status': 'PASS', 'online_ms': 1000.0}
                for candidate in ids},
            'calibration_by_candidate': {
                candidate: {'status': 'PASS_10_PERCENT'} for candidate in ids},
        })
    return {
        'schema': 'n83.target-cancel-tier1-holdout-validated/v1',
        'source_sha256': run.SOURCE_SHA256,
        'targets_sha256': run.TARGETS_SHA256,
        'evidence_status': 'PASS_INDEPENDENT_REPLAY',
        'different_cloud_instances': True,
        'candidate_ids': ids,
        'choice': run.choose()[1],
        'host_strata': [
            {'point_indices': list(range(6 * i, 6 * i + 6)),
             'tier_status': 'PASS_TIER1',
             'aa_status': 'PASS_10_PERCENT_PER_PAIR'}
            for i in range(4)],
        'points': points,
    }


def main():
    data = sample()
    selected = selector.select(data)
    cancellation = run.ARMS[1][3]
    assert selected['online']['selected_candidate_id'] == cancellation
    assert selected['cold']['selected_candidate_id'] == cancellation
    assert selected['online']['comparisons'][0]['two_sided_sign_flip_p'] == 2 / 2**24
    assert selected['cold']['comparisons'][0]['two_sided_sign_flip_p'] == 2 / 2**24
    data['different_cloud_instances'] = False
    rejected = selector.select(data)
    assert rejected['online']['selected_candidate_id'] is None
    assert rejected['admission_status'] == 'SOURCE_TARGET_OR_REPLAY_GATE'
    print(json.dumps({'status': 'PASS_PROTOCOL_SELECTION_AND_PROVIDER_GATE',
                      'candidate_id': cancellation,
                      'exact_p': selected['online']['comparisons'][0][
                          'two_sided_sign_flip_p']}, sort_keys=True))


if __name__ == '__main__':
    main()
