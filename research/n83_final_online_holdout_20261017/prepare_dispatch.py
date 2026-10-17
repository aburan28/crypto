#!/usr/bin/env python3
"""Seal upstream evidence and all executable inputs before cloud dispatch."""

import argparse
import hashlib
import json
from pathlib import Path

import arms
import run

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FILES = ('generate.py', 'arms.py', 'run.py', 'verify.py', 'select.py',
         'aggregate.py', 'prepare_dispatch.py', 'check_dispatch.py',
         'test_choice_selector.py',
         '../n83_runner_identity_replay_20261017/jobs.py',
         '../../.github/workflows/n83-final-online-holdout.yml')
MANIFEST_DIRS = (
    'research/n83_k992_cold_followup_20261017',
    'research/n83_smallk_current_source_20261017',
    'research/n83_k2543_current_source_20261017',
    'research/n83_k3072_current_source_20261017',
    'research/n83_e2e_final_holdout_20261016',
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def within_repo(path):
    return str(path.resolve().relative_to(ROOT))


def manifest(candidate_id):
    for directory in MANIFEST_DIRS:
        path = ROOT / directory / f'candidate_{candidate_id}.json'
        if path.exists():
            if json.loads(path.read_text()).get('candidate_id') != candidate_id:
                raise ValueError('candidate manifest identity differs')
            return {'path': within_repo(path), 'sha256': digest(path)}
    raise ValueError(f'missing exact candidate manifest {candidate_id}')


def record(paths):
    chosen, choice = arms.choose(*paths)
    frozen_source_hashes, _ = run.fixture()
    for path in paths[::2]:
        panel = json.loads(path.read_text())
        strata = panel.get('host_strata')
        if (not isinstance(panel.get('run_id'), int) or panel['run_id'] <= 0 or
                not isinstance(panel.get('head_sha'), str) or
                len(panel['head_sha']) != 40 or
                not isinstance(strata, list) or len(strata) != 2 or
                any(not isinstance(item.get('raw_artifact_sha256'), str) or
                    len(item['raw_artifact_sha256']) != 64
                    for item in strata)):
            raise ValueError('upstream provider archive provenance is incomplete')
    item = {
        'schema': 'n83.final-online-holdout-dispatch/v1',
        'source_sha256': arms.SOURCE_SHA256,
        'targets_sha256': digest(HERE / 'targets.jsonl'),
        'answers_sha256': digest(HERE / 'answers.jsonl'),
        'fixture_manifest_sha256': digest(HERE / 'fixture_manifest.json'),
        'panel_paths': [within_repo(path) for path in paths],
        'panel_sha256': [digest(path) for path in paths],
        'choice': choice,
        'ordered_arms': [list(arm) for arm in chosen],
        'candidate_manifests': {arm[3]: manifest(arm[3]) for arm in chosen},
        'frozen_source_hashes': frozen_source_hashes,
        'input_sha256': {path: digest((HERE / path).resolve()) for path in FILES},
        'pair_rho_seed_base': 2026102900,
        'ic_cap_small_seconds': 480,
        'ic_cap_large_seconds': 900,
        'rho_cap_seconds': 120,
        'producer_cap_seconds': 190 * 60,
        'points_per_host': 6,
    }
    if (item['targets_sha256'] !=
            '4d167f8b2b85bcf6c599494daf3b3a61210a9dc8c0d65ad3ebcdf07cd6d5c4c7' or
            item['answers_sha256'] !=
            'f94a313ea6e7c414a4fc85c3ec55e384f6c6be14febbe3304d4b9dde1e54e484'):
        raise ValueError('sealed holdout fixture differs')
    return item


def main():
    parser = argparse.ArgumentParser()
    for name in ('five-validated', 'five-selection',
                 'k992-validated', 'k992-selection'):
        parser.add_argument(f'--{name}', type=Path, required=True)
    args = parser.parse_args()
    paths = tuple(getattr(args, name).resolve() for name in
                  ('five_validated', 'five_selection',
                   'k992_validated', 'k992_selection'))
    item = record(paths)
    path = HERE / 'dispatch.json'
    blob = json.dumps(item, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError('sealed dispatch differs')
    path.write_text(blob)
    print(json.dumps(item, sort_keys=True))


if __name__ == '__main__':
    main()
