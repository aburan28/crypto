#!/usr/bin/env python3
"""Check sealed five-arm schedule, candidate binding and selection gates."""
import importlib.util
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module('smallk_test_run', HERE / 'run.py')
selector = module('smallk_test_select', HERE / 'select.py')


def synthetic_validated():
    ids = [arm[3] for arm in (*runner.FIXED_ARMS, next(
        arm for arm in runner.finalist.ARMS if arm[3] == runner.FALLBACK_ID))]
    times = [100.0, 20.0, 90.0, 95.0, 110.0]
    points = []
    for index in range(10):
        points.append({
            'point_index': index,
            'arms': [{'candidate_id': candidate_id,
                      'status': 'PASS_INDEPENDENT_REPLAY',
                      'resource_status': 'PASS',
                      'one_target_online_ms': value,
                      'cold_process_ms': value * 1000}
                     for candidate_id, value in zip(ids, times)],
            'rho_by_candidate': {candidate_id: {
                'status': 'PASS_INDEPENDENT_REPLAY',
                'resource_status': 'PASS', 'online_ms': 500.0}
                for candidate_id in ids},
            'calibration_by_candidate': {candidate_id: {'status': 'PASS_10_PERCENT'}
                                         for candidate_id in ids},
        })
    return {
        'schema': 'n83.smallk-challenge-validated/v1',
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
    for arm, study, filename in (
            (runner.K1248, 'n83_smallk_current_source_20261017',
             f'candidate_{runner.K1248[3]}.json'),
            (runner.K1280, 'n83_smallk_current_source_20261017',
             f'candidate_{runner.K1280[3]}.json'),
            (runner.K2543, 'n83_k2543_current_source_20261017',
             f'candidate_{runner.K2543[3]}.json'),
            (runner.K3072, 'n83_k3072_current_source_20261017',
             f'candidate_{runner.K3072[3]}.json')):
        manifest = json.loads((ROOT / 'research' / study / filename).read_text())
        assert manifest['candidate_id'] == arm[3]
        assert manifest['factor_base']['actual_usable_points'] == 166 * arm[1]
    arms = (*runner.FIXED_ARMS, next(
        arm for arm in runner.finalist.ARMS if arm[3] == runner.FALLBACK_ID))
    for half in (0, 1):
        events = runner.schedule(half, arms)
        assert len(events) == 100
        assert Counter(event['kind'] for event in events) == {
            'ic': 25, 'rho': 25, 'calibration': 50}
        for candidate in (arm[0] for arm in arms):
            positions = []
            for index in range(half * 5, half * 5 + 5):
                arm_events = [event for event in events
                              if event['point_index'] == index and
                              event['candidate'] == candidate]
                fixed = next(i for i, arm in enumerate(arms) if arm[0] == candidate)
                expected_pair = (['ic', 'rho'] if (index + fixed) % 2 == 0
                                 else ['rho', 'ic'])
                assert [event['kind'] for event in arm_events] == [
                    'calibration', *expected_pair, 'calibration']
                order = [event['candidate'] for event in events
                         if event['point_index'] == index and event['kind'] == 'ic']
                positions.append(order.index(candidate))
            assert sorted(positions) == list(range(5))
    data = synthetic_validated()
    selected = selector.select(data)
    assert selected['admission_status'] == 'PASS'
    assert selected['online']['selected_candidate_id'] == runner.K1280[3]
    assert selected['cold']['selected_candidate_id'] == runner.K1280[3]
    assert len(selected['online']['comparisons']) == 10
    assert selected['online']['bonferroni_alpha_per_pair'] == 0.005
    data['points'][0]['arms'][0]['status'] = 'VERIFICATION_FAILED'
    rejected = selector.select(data)
    assert rejected['admission_status'] == 'IC_REPLAY_RESOURCE_OR_TIMING_GATE'
    assert rejected['online']['selected_candidate_id'] is None
    print(json.dumps({'status': 'PASS', 'events_per_host': 100,
                      'candidate_count': 5, 'paired_comparisons': 10}, sort_keys=True))


if __name__ == '__main__':
    main()
