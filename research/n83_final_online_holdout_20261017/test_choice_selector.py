#!/usr/bin/env python3
"""Check the preregistered fallbacks and exact twelve-point selector."""

import json
import importlib.util
import tempfile
from pathlib import Path

import arms
import run

spec = importlib.util.spec_from_file_location(
    'n83_final_holdout_selector', Path(__file__).with_name('select.py'))
holdout_select = importlib.util.module_from_spec(spec)
spec.loader.exec_module(holdout_select)


def main():
    k_root = Path(__file__).resolve().parents[1] / 'n83_k992_cold_followup_20261017/results'
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        five_validated = base / 'five_validated.json'
        five_selection = base / 'five_selection.json'
        record = {
            'schema': 'n83.smallk-challenge-validated/v1',
            'source_sha256': arms.SOURCE_SHA256,
            'targets_sha256': arms.five_run.TARGETS_SHA256,
            'evidence_status': 'VERIFICATION_FAILED',
            'run_id': 1, 'head_sha': '0' * 40,
            'candidate_ids': [arm[3] for arm in arms.five_run.FIXED_ARMS] +
                             [arms.five_run.FALLBACK_ID],
            'choice': {'comparison_candidate_id': arms.five_run.FALLBACK_ID},
        }
        five_validated.write_text(json.dumps(record))
        five_selection.write_text(json.dumps(arms.five_select.select(record)))
        chosen, choice = arms.choose(
            five_validated, five_selection, k_root / 'validated_rows.json',
            k_root / 'selection.json')
        assert [item[3] for item in chosen] == [
            arms.K992_ID, arms.K1107_ID, arms.K1280_ID, arms.K1024_ID]
        assert choice['five_arm_fallback_used']
        assert choice['k992_fallback_used']
        repeated = arms.unique_names((arms.k992_run.K1107,
                                      arms.five_run.finalist.ARMS[2]))
        assert repeated[0][0] != repeated[1][0]
        assert repeated[0][3] != repeated[1][3]
        for half in (0, 1):
            events = run.schedule(half, chosen)
            assert len(events) == 96
            assert {event['point_index'] for event in events} == set(
                range(half * 6, half * 6 + 6))
        assert holdout_select.sign_flip_p([1.0] * 12) == 2 / 4096
        assert holdout_select.sign_flip_p([1.0, -1.0] * 6) == 1.0
        ids = [item[3] for item in chosen]
        admitted = {
            'schema': 'n83.final-online-holdout-validated/v1',
            'source_sha256': arms.SOURCE_SHA256,
            'targets_sha256': holdout_select.TARGETS_SHA256,
            'evidence_status': 'PASS_INDEPENDENT_REPLAY',
            'candidate_ids': ids, 'choice': choice,
            'host_strata': [
                {'point_indices': list(range(half * 6, half * 6 + 6)),
                 'tier_status': 'PASS_TIER1',
                 'aa_status': 'PASS_10_PERCENT_PER_PAIR'} for half in (0, 1)],
            'points': [{
                'point_index': index,
                'arms': [
                    {'candidate_id': candidate, 'status': 'PASS_INDEPENDENT_REPLAY',
                     'resource_status': 'PASS',
                     'one_target_online_ms': 10.0 * (rank + 1),
                     'cold_process_ms': 1000.0 * (rank + 1)}
                    for rank, candidate in enumerate(ids)],
                'rho_by_candidate': {
                    candidate: {'status': 'PASS_INDEPENDENT_REPLAY',
                                'resource_status': 'PASS', 'online_ms': 100.0}
                    for candidate in ids},
                'calibration_by_candidate': {
                    candidate: {'status': 'PASS_10_PERCENT'} for candidate in ids},
            } for index in range(12)],
        }
        selection = holdout_select.select(admitted)
        assert selection['admission_status'] == 'PASS'
        assert selection['online']['selected_candidate_id'] == ids[0]
        assert selection['cold']['selected_candidate_id'] == ids[0]
        assert len(selection['online']['comparisons']) == 6
        admitted['points'][0]['rho_by_candidate'][ids[0]]['status'] = 'VERIFICATION_FAILED'
        assert holdout_select.select(admitted)['admission_status'] == \
            'RHO_REPLAY_RESOURCE_OR_TIMING_GATE'
        broken = dict(record, source_sha256='0' * 64)
        five_validated.write_text(json.dumps(broken))
        try:
            arms.choose(five_validated, five_selection,
                        k_root / 'validated_rows.json', k_root / 'selection.json')
        except ValueError:
            pass
        else:
            raise AssertionError('source mutation was admitted')
    print(json.dumps({'status': 'PASS', 'fallback_arm_count': 4,
                      'minimum_exact_p': 2 / 4096}, sort_keys=True))


if __name__ == '__main__':
    main()
