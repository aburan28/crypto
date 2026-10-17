#!/usr/bin/env python3
"""Freeze ten N83 public points before the current-source comparison."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / 'research/n83_e2e_final_holdout_20261016/generate.py'
PARENT_SHA256 = '1ed595d1265574c9613657ea89515f44455223ef9520d2790917abbe7c121dce'
FINAL_TARGETS = ROOT / 'research/n83_e2e_final_holdout_20261016/targets.jsonl'
FINAL_TARGETS_SHA256 = '407cd248797100214ff039193bb04bf89b5c2708a8cd339ca98f7aca98174062'
FINAL_ANSWERS = ROOT / 'research/n83_e2e_final_holdout_20261016/answers.jsonl'
FINAL_ANSWERS_SHA256 = 'fe771284a3b0e5c0301309f678751f0630f0fa23844dd0ce1dbc9fa47f243076'
DOMAIN = b'n83-k3072-current-source-ten-20261017-v1'
COUNT = 10


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module():
    if digest(PARENT) != PARENT_SHA256:
        raise ValueError('parent fixture generator differs')
    spec = importlib.util.spec_from_file_location('n83_k3072_prior_fixture', PARENT)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def canonical_lines(rows):
    return ''.join(json.dumps(row, separators=(',', ':')) + '\n' for row in rows)


def expected_rows():
    parent = module()
    precision = parent.precision_module()
    group = precision.group
    if (digest(FINAL_TARGETS) != FINAL_TARGETS_SHA256 or
            digest(FINAL_ANSWERS) != FINAL_ANSWERS_SHA256):
        raise ValueError('finalist fixture differs')
    final_targets, final_answers = parent.expected_rows()
    if (canonical_lines(final_targets) != FINAL_TARGETS.read_text() or
            canonical_lines(final_answers) != FINAL_ANSWERS.read_text()):
        raise ValueError('finalist fixture replay differs')
    prior = precision.prior_points()
    precision_targets, _ = precision.expected_rows()
    prior.update(tuple(item) for item in precision_targets)
    prior.update(tuple(item) for item in final_targets)
    targets, answers = [], []
    for index in range(COUNT):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, 'little')).digest()
        scalar = 1 + int.from_bytes(seed, 'little') % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if (tuple(encoded) in prior or not group.on_curve(point) or
                group.multiply(point, group.R) is not None):
            raise ValueError(f'fresh point {index} repeats or fails subgroup replay')
        prior.add(tuple(encoded))
        targets.append(encoded)
        answers.append({'fixture_index': index, 'point': encoded,
                        'scalar': str(scalar), 'seed_sha256': seed.hex()})
    return targets, answers


def main(verify):
    targets, answers = expected_rows()
    target_blob = canonical_lines(targets)
    answer_blob = canonical_lines(answers)
    manifest = {
        'schema': 'n83.k3072-current-source-fixture/v1',
        'domain': DOMAIN.decode(), 'target_count': COUNT,
        'parent_generator_sha256': PARENT_SHA256,
        'final_targets_sha256': FINAL_TARGETS_SHA256,
        'final_answers_sha256': FINAL_ANSWERS_SHA256,
        'targets_sha256': hashlib.sha256(target_blob.encode()).hexdigest(),
        'answers_sha256': hashlib.sha256(answer_blob.encode()).hexdigest(),
    }
    files = {HERE / 'targets.jsonl': target_blob,
             HERE / 'answers.jsonl': answer_blob,
             HERE / 'fixture_manifest.json': json.dumps(manifest, sort_keys=True,
                                                        indent=2) + '\n'}
    for path, blob in files.items():
        if path.exists():
            if path.read_text() != blob:
                raise ValueError(f'sealed fixture differs: {path.name}')
        elif verify:
            raise ValueError(f'missing sealed fixture: {path.name}')
        else:
            path.write_text(blob)
    print(json.dumps({'mode': 'sealed_replay' if verify else 'freeze',
                      **manifest}, sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify-sealed', action='store_true')
    args = parser.parse_args()
    main(args.verify_sealed)
