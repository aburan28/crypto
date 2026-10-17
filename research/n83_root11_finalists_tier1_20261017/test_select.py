#!/usr/bin/env python3
"""Exercise the frozen sign-test decision and Tier 1 admission gates."""

import json

import arms
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    'n83_finalists_selector', Path(__file__).with_name('select.py'))
selector = importlib.util.module_from_spec(spec)
spec.loader.exec_module(selector)
chosen, choice = arms.choose()
ids = [arm[3] for arm in chosen]


def data(faster_second):
    points = []
    for index in range(32):
        first_ms, second_ms = ((20.0, 10.0) if index in faster_second
                               else (10.0, 20.0))
        points.append({
            'point_index': index,
            'arms': [
                {'candidate_id': ids[0], 'status': 'PASS_INDEPENDENT_REPLAY',
                 'resource_status': 'PASS', 'one_target_online_ms': first_ms,
                 'cold_process_ms': 1000.0},
                {'candidate_id': ids[1], 'status': 'PASS_INDEPENDENT_REPLAY',
                 'resource_status': 'PASS', 'one_target_online_ms': second_ms,
                 'cold_process_ms': 1100.0}],
            'rho_by_candidate': {
                candidate: {'status': 'PASS_INDEPENDENT_REPLAY',
                            'resource_status': 'PASS', 'online_ms': 1000.0}
                for candidate in ids},
            'calibration_by_candidate': {
                candidate: {'status': 'PASS_10_PERCENT'} for candidate in ids},
        })
    return {
        'schema': 'n83.root11-finalists-validated/v1',
        'source_sha256': selector.SOURCE_SHA256,
        'targets_sha256': selector.TARGETS_SHA256,
        'evidence_status': 'PASS_INDEPENDENT_REPLAY',
        'different_cloud_instances': True,
        'same_candidate_bases_across_strata': True,
        'candidate_ids': ids, 'choice': choice,
        'host_strata': [
            {'point_indices': list(range(i * 8, i * 8 + 8)),
             'tier_status': 'PASS_TIER1',
             'aa_status': 'PASS_10_PERCENT_PER_PAIR'}
            for i in range(4)],
        'points': points,
    }


all_second = selector.select(data(set(range(32))))
assert all_second['online']['selected_candidate_id'] == ids[1]
assert all_second['cold']['selected_candidate_id'] == ids[0]
assert all_second['online']['two_sided_exact_sign_p'] == 2 / 2**32
assert selector.select(data(set(range(23))))['online']['selected_candidate_id'] == ids[1]
assert selector.select(data(set(range(22))))['online']['status'] == 'BOUNDED_UNKNOWN'
assert selector.select(data(set()))['online']['selected_candidate_id'] == ids[0]
noisy = data(set(range(32)))
noisy['host_strata'][0]['aa_status'] = 'NOISY_OR_INCOMPLETE'
assert selector.select(noisy)['admission_status'] == 'HOST_STRATA_GATE'
print(json.dumps({'status': 'PASS', 'cases': 6}))
