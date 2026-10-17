#!/usr/bin/env python3
"""Freeze the K4096 indexed-layout screen and prospective holdout choice."""

import argparse
import json
from pathlib import Path

import run

SOURCE_SHA256 = 'c73985a15892ea1674a3ffee8f3fdb00fdfd4728a41b3edda4645b054d573050'
HOLDOUT = run.ROOT / 'research/n83_indexed4_id30_holdout_20261017'


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'frozen receipt differs: {path}')
    path.write_text(blob)


def seal(binary):
    if run.digest(run.SOURCE) != SOURCE_SHA256:
        raise ValueError('four-byte ID30 source differs')
    files = run.frozen_paths(binary)
    plan = {
        'schema': 'n83.indexed4-id30-screen-plan/v1',
        'role': 'ordinary_local_candidate_diagnostic',
        'arm_names': list(run.ARMS),
        'k': run.K,
        'target_count': 128,
        'binary_path': str(binary),
        'cap_seconds_per_arm': run.CAP,
        'sha256': {name: run.digest(path) for name, path in files.items()},
    }
    write_once(run.HERE / 'plan.json', plan)
    paths = [run.HERE / 'plan.json', HOLDOUT / 'generate.py',
             HOLDOUT / 'manifest.json', HOLDOUT / 'targets.jsonl',
             HOLDOUT / 'answers.jsonl', HOLDOUT / 'PROTOCOL.md',
             HOLDOUT / 'choose.py']
    choice = {
        'schema': 'n83.indexed4-id30-holdout-selection-plan/v1',
        'arm_names': list(run.ARMS),
        'k': run.K,
        'choice_rule':
            'dispatch_only_if_both_complete_and_four_byte_geometric_online_lower',
        'screen_plan_sha256': run.digest(run.HERE / 'plan.json'),
        'sha256': {str(path.relative_to(run.ROOT)): run.digest(path)
                   for path in paths},
    }
    write_once(HOLDOUT / 'selection_plan.json', choice)
    return {'screen_plan_sha256': choice['screen_plan_sha256'],
            'selection_plan_sha256': run.digest(HOLDOUT / 'selection_plan.json')}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(seal(args.binary.resolve()), sort_keys=True))


if __name__ == '__main__':
    main()
