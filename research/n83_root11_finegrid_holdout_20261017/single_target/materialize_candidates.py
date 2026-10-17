#!/usr/bin/env python3
"""Bind each independently replayed single-target screen to an exact IC1 base."""

import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path('/Volumes/SSD990/crypto/worktrees/n83-final-online-holdout-20261017')
HERE = Path(__file__).resolve().parent
K_VALUES = (2752, 2880, 2944)
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
SOURCE_COMMIT = '416489e09c3d8885b6c12f262935e6ba29373460'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


prior = module('finegrid_manifest_prior',
               ROOT / 'research/n83_e2e_final_holdout_20261016/verify.py')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'candidate manifest changed: {path}')
    path.write_text(blob)


def main():
    if (digest(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') != SOURCE_SHA256 or
            json.loads((HERE / 'plan.json').read_text())['screen_script_sha256'] !=
            digest(HERE / 'screen.py')):
        raise ValueError('single-target source or screen plan differs')
    prior.finegrid.panel.SOURCE_COMMIT = SOURCE_COMMIT
    field, curve = prior.finegrid.panel.curve_record()
    rows = []
    for k in K_VALUES:
        folder = HERE / f'k{k}'
        receipt = json.loads((folder / 'screen_result.json').read_text())
        summary = json.loads((folder / 'stdout.jsonl').read_text())
        if (receipt['k'] != k or
                receipt['producer_status'] != 'PROCESS_EXIT_0' or
                receipt['correctness'] != 'PASS_INDEPENDENT_SCALAR_BASE_RELATION' or
                summary['base_hash'] != receipt['base_hash'] or
                summary['rank'] != k or summary['factor_base_points'] != 166 * k):
            raise ValueError(f'K{k}: single-target replay is incomplete')
        manifest = prior.finegrid.dynamic.candidate(
            k, 'root11_dynamic', summary, folder / 'base_dump.jsonl',
            field, curve)
        if (manifest['factor_base']['source_base_hash_blake3'] !=
                summary['base_hash'] or
                manifest['factor_base']['actual_usable_points'] != 166 * k or
                manifest['factor_base']['effective_columns'] != k):
            raise ValueError(f'K{k}: candidate base identity differs')
        path = HERE / f'candidate_k{k}.json'
        write_once(path, manifest)
        rows.append({'k': k, 'candidate_id': manifest['candidate_id'],
                     'manifest_sha256': digest(path),
                     'single_target_receipt_sha256':
                         digest(folder / 'screen_result.json')})
    write_once(HERE / 'candidate_manifest_receipt.json', {
        'schema': 'n83.root11-finegrid-candidate-manifests/v1',
        'source_sha256': SOURCE_SHA256,
        'screen_plan_sha256': digest(HERE / 'plan.json'),
        'rows': rows})
    print(json.dumps({'status': 'PASS', 'candidate_ids':
                      [row['candidate_id'] for row in rows]}))


if __name__ == '__main__':
    main()
