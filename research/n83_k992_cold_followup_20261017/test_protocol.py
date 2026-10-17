#!/usr/bin/env python3
"""Check the sealed K992 fixture, three-arm schedule and selection gates."""

import importlib.util
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module('k992_test_run', HERE / 'run.py')
selector = module('k992_test_select', HERE / 'select.py')


def synthetic_validated():
    ids = [arm[3] for arm in runner.ARMS]
    points = []
    for index in range(10):
        points.append({
            'point_index': index,
            'arms': [{'candidate_id': candidate_id,
                      'status': 'PASS_INDEPENDENT_REPLAY',
                      'resource_status': 'PASS',
                      'one_target_online_ms': online,
                      'cold_process_ms': cold}
                     for candidate_id, online, cold in
                     zip(ids, (50.0, 100.0, 120.0),
                         (90_000.0, 100_000.0, 110_000.0))],
            'rho_by_candidate': {candidate_id: {
                'status': 'PASS_INDEPENDENT_REPLAY',
                'resource_status': 'PASS', 'online_ms': 500.0}
                for candidate_id in ids},
            'calibration_by_candidate': {
                candidate_id: {'status': 'PASS_10_PERCENT'}
                for candidate_id in ids},
        })
    return {
        'schema': 'n83.k992-cold-followup-validated/v1',
        'source_sha256': runner.SOURCE_SHA256,
        'targets_sha256': runner.TARGETS_SHA256,
        'evidence_status': 'PASS_INDEPENDENT_REPLAY',
        'candidate_ids': ids,
        'host_strata': [{'point_indices': list(range(half * 5, half * 5 + 5)),
                         'tier_status': 'PASS_TIER1',
                         'aa_status': 'PASS_10_PERCENT_PER_PAIR'}
                        for half in (0, 1)],
        'points': points,
    }


def main():
    runner.fixture()
    for arm in runner.ARMS:
        manifest = json.loads((HERE / f'candidate_{arm[3]}.json').read_text())
        assert manifest['candidate_id'] == arm[3]
        assert manifest['factor_base']['actual_usable_points'] == 166 * arm[1]
        assert manifest['implementation']['source_sha256']['ic'] == runner.SOURCE_SHA256
    for half in (0, 1):
        events = runner.schedule(half)
        assert len(events) == 60
        assert Counter(event['kind'] for event in events) == {
            'ic': 15, 'rho': 15, 'calibration': 30}
        for candidate in (arm[0] for arm in runner.ARMS):
            positions = []
            for index in range(half * 5, half * 5 + 5):
                pair = [event for event in events
                        if event['point_index'] == index and
                        event['candidate'] == candidate]
                fixed = next(i for i, arm in enumerate(runner.ARMS)
                             if arm[0] == candidate)
                expected = (['ic', 'rho'] if (index + fixed) % 2 == 0
                            else ['rho', 'ic'])
                assert [event['kind'] for event in pair] == [
                    'calibration', *expected, 'calibration']
                order = [event['candidate'] for event in events
                         if event['point_index'] == index and event['kind'] == 'ic']
                positions.append(order.index(candidate))
            assert set(positions) == {0, 1, 2}
    data = synthetic_validated()
    selected = selector.select(data)
    assert selected['admission_status'] == 'PASS'
    assert selected['online']['selected_candidate_id'] == runner.K992[3]
    assert selected['cold']['selected_candidate_id'] == runner.K992[3]
    assert len(selected['online']['comparisons']) == 3
    assert selected['online']['bonferroni_alpha_per_pair'] == 0.05 / 3
    data['points'][0]['arms'][0]['status'] = 'VERIFICATION_FAILED'
    rejected = selector.select(data)
    assert rejected['admission_status'] == 'IC_REPLAY_RESOURCE_OR_TIMING_GATE'
    assert rejected['online']['selected_candidate_id'] is None
    print(json.dumps({'status': 'PASS', 'events_per_host': 60,
                      'candidate_count': 3, 'paired_comparisons': 3},
                     sort_keys=True))


if __name__ == '__main__':
    main()
