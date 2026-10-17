#!/usr/bin/env python3
"""Replay the five-arm aggregate with provider-backed runner identity.

This successor leaves the dispatched aggregation code and its original
result intact. It only corrects a rejected shared workflow-label identity
after all remaining evidence gates and the provider's distinct concurrent
runner assignments have been checked.
"""

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


frozen = module('n83_smallk_frozen_aggregate', HERE / 'aggregate.py')
identity = module('n83_smallk_provider_runner_identity',
                  ROOT / 'research/n83_runner_identity_replay_20261017/jobs.py')


def main():
    parser = argparse.ArgumentParser()
    for half in (0, 1):
        parser.add_argument(f'--half{half}-results', type=Path, required=True)
        parser.add_argument(f'--half{half}-zip', type=Path, required=True)
        parser.add_argument(f'--half{half}-provider-sha256', required=True)
    parser.add_argument('--provider-jobs', type=Path, required=True)
    parser.add_argument('--original-aggregate', type=Path, required=True)
    parser.add_argument('--finalist-validated', type=Path, required=True)
    parser.add_argument('--finalist-selection', type=Path, required=True)
    parser.add_argument('--run-id', type=int, required=True)
    parser.add_argument('--run-attempt', type=int, default=1)
    parser.add_argument('--head-sha', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    halves = []
    for half in (0, 1):
        results = getattr(args, f'half{half}_results').resolve()
        item = frozen.bind_half(
            half, results, getattr(args, f'half{half}_zip').resolve(),
            getattr(args, f'half{half}_provider_sha256'), args.run_id,
            args.run_attempt, args.head_sha,
            args.finalist_validated.resolve(),
            args.finalist_selection.resolve(), True)
        item['results'] = str(results)
        halves.append(item)
    original, original_selection = frozen.aggregate(
        halves, args.run_id, args.run_attempt, args.head_sha)
    original_dir = args.original_aggregate.resolve()
    if (original != json.loads((original_dir / 'validated_rows.json').read_text()) or
            original_selection != json.loads((original_dir / 'selection.json').read_text())):
        raise ValueError('frozen original aggregate differs from retained result')
    if (original['evidence_status'] != 'VERIFICATION_FAILED' or
            original['different_cloud_instances'] is not False or
            original['same_source_and_binary'] is not True or
            original['same_candidate_bases_across_strata'] is not True or
            len(original['points']) != 10 or
            {point['point_index'] for point in original['points']} != set(range(10)) or
            halves[0]['host']['parent_instance_id'] !=
            f'gh-{args.run_id}-{args.run_attempt}-smallk-challenge-half' or
            halves[0]['host']['parent_instance_id'] !=
            halves[1]['host']['parent_instance_id'] or
            any(item['analysis']['status'] != 'HALF_ADMITTED' or
                item['analysis']['evidence_status'] != 'PASS_INDEPENDENT_REPLAY' or
                item['analysis']['tier_status'] != 'PASS_TIER1' or
                item['analysis']['aa_status'] != 'PASS_10_PERCENT_PER_PAIR' or
                item['local_replay_exit'] != 0 for item in halves)):
        raise ValueError('frozen aggregate has another failed evidence gate')
    proof = identity.verify(args.provider_jobs.resolve(), halves, args.run_id,
                            args.run_attempt, 'smallk-challenge-half')
    corrected = dict(original)
    corrected['different_cloud_instances'] = True
    corrected['evidence_status'] = 'PASS_INDEPENDENT_REPLAY'
    corrected['identity_replay'] = proof
    corrected['original_failed_aggregate_sha256'] = {
        name: frozen.sha256(original_dir / name)
        for name in ('validated_rows.json', 'selection.json')}
    selection = frozen.selector.select(corrected)
    if selection['admission_status'] != 'PASS':
        raise ValueError('corrected identity did not admit exact selector')
    args.output.mkdir(parents=True, exist_ok=True)
    frozen.write_once(args.output / 'validated_rows.json', corrected)
    frozen.write_once(args.output / 'selection.json', selection)
    print(json.dumps({'evidence_status': corrected['evidence_status'],
                      'runner_identity': proof['status'],
                      'online_selected': selection['online']['selected_candidate_id'],
                      'cold_selected': selection['cold']['selected_candidate_id']},
                     sort_keys=True))


if __name__ == '__main__':
    main()
