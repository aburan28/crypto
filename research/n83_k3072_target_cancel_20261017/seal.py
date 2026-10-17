#!/usr/bin/env python3
"""Freeze source, binary, candidates and inputs before the K3072 screen."""

import argparse
import json
from pathlib import Path

import run


def seal(binary):
    paths = run.frozen_files(binary)
    if (run.digest(paths['source']) !=
            '8109fb33374ce2aa66c6a5a0b6e144d3b0143168ac756f9e0bd70af60a911829' or
            run.digest(binary) !=
            'df8340eb76ac90d2520e0202159505426b9d7879978f10c597731c4bfe73045a'):
        raise ValueError('checked source or compiled binary differs')
    plan = {
        'schema': 'n83.k3072-target-cancel-plan/v1',
        'workload_role': 'ordinary_local_shared_setup_stage_diagnostic',
        'k': 3072,
        'target_count': 128,
        'policies': {name: list(value) for name, value in run.POLICIES.items()},
        'candidate_ids': json.loads((run.HERE / 'candidate_ids.json').read_text()),
        'binary_path': str(binary),
        'cap_seconds': 1200,
        'sha256': {name: run.digest(path) for name, path in paths.items()},
    }
    output = run.HERE / 'plan.json'
    blob = json.dumps(plan, sort_keys=True, indent=2) + '\n'
    if output.exists() and output.read_text() != blob:
        raise ValueError('frozen plan differs')
    output.write_text(blob)
    return plan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, required=True)
    args = parser.parse_args()
    plan = seal(args.binary.resolve())
    print(json.dumps({'plan_sha256': run.digest(run.HERE / 'plan.json'),
                      'candidate_ids': plan['candidate_ids']}, sort_keys=True))


if __name__ == '__main__':
    main()
