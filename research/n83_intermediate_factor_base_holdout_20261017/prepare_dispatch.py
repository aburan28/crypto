#!/usr/bin/env python3
"""Seal the published arm choice and executable inputs before Tier 1 work."""

import hashlib
import json
from pathlib import Path

import arms
import run

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FILES = (
    'generate.py', 'manifest.json', 'targets.jsonl', 'answers.jsonl',
    'seal.json', 'PROTOCOL.md', 'RESULT.md', 'arms.py', 'select_batch.py',
    'test_select_batch.py', 'run.py', 'verify.py', 'select.py',
    'test_select.py', 'provider_jobs.py', 'test_provider_jobs.py', 'aggregate.py',
    'prepare_dispatch.py', 'check_dispatch.py',
    'results/batch_selection.json', 'results/batch_full_replay.json',
    'results/verifier_correction.json',
    'results/batch_k2304_result.json', 'results/batch_k2543_result.json',
    'results/batch_k2688_result.json', 'results/batch_k2816_result.json',
    'results/batch_k3072_result.json',
    '../../.github/workflows/n83-intermediate-factor-base-holdout.yml',
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def record():
    chosen, choice = arms.choose()
    frozen_hashes, _ = run.fixture()
    if len(chosen) != 4 or choice['selected_fourth_k'] != 2816:
        raise ValueError('sealed candidate selection differs')
    manifests = {}
    for k, candidate_id in arms.IDS.items():
        path = arms.manifest_path(k)
        if json.loads(path.read_text())['candidate_id'] != candidate_id:
            raise ValueError('exact candidate manifest differs')
        manifests[candidate_id] = {'path': relative(path),
                                   'sha256': digest(path)}
    return {
        'schema': 'n83.intermediate-holdout-dispatch/v1',
        'source_sha256': arms.SOURCE_SHA256,
        'targets_sha256': run.TARGETS_SHA256,
        'answers_sha256': run.ANSWERS_SHA256,
        'choice': choice,
        'ordered_arms': [list(arm) for arm in chosen],
        'candidate_manifests': manifests,
        'frozen_source_hashes': frozen_hashes,
        'input_sha256': {name: digest((HERE / name).resolve()) for name in FILES},
        'pair_rho_seed_base': run.PAIR_RHO_SEED_BASE,
        'ic_cap_seconds': run.IC_CAP_SECONDS,
        'rho_cap_seconds': run.RHO_CAP_SECONDS,
        'producer_cap_seconds': run.GLOBAL_PRODUCER_CAP_SECONDS,
        'points_per_host': 5,
        'host_strata': 4,
    }


def main():
    item = record()
    path = HERE / 'dispatch.json'
    blob = json.dumps(item, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError('sealed dispatch differs')
    path.write_text(blob)
    print(json.dumps({'status': 'SEALED',
                      'candidate_ids': item['choice']['ordered_candidate_ids']}))


if __name__ == '__main__':
    main()
