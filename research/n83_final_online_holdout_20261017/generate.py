#!/usr/bin/env python3
"""Seal twelve N83 public points for the final online factor-base comparison."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / 'research/n83_k992_cold_followup_20261017/generate.py'
PARENT_SHA256 = '514329656510e88d3316b688eef4d15104c64a5be7edb42f56ace269461d1dec'
PARENT_TARGETS = PARENT.parent / 'targets.jsonl'
PARENT_TARGETS_SHA256 = '89c7208dd4fc7a5621cf727b3fac9a6714f04f66a312e516b418c31764e7547e'
PARENT_ANSWERS = PARENT.parent / 'answers.jsonl'
PARENT_ANSWERS_SHA256 = '363adc7f9bded3f4a25a80c7b27da5a2305afb30ee6741932cbe76f455c2865b'
DOMAIN = b'n83-final-online-holdout-twelve-20261017-v1'
COUNT = 12


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def canonical_lines(rows):
    return ''.join(json.dumps(row, separators=(',', ':')) + '\n' for row in rows)


def expected_rows():
    if (digest(PARENT) != PARENT_SHA256 or
            digest(PARENT_TARGETS) != PARENT_TARGETS_SHA256 or
            digest(PARENT_ANSWERS) != PARENT_ANSWERS_SHA256):
        raise ValueError('parent fixture source or bytes differ')
    parent = module('n83_final_online_parent', PARENT)
    parent_targets, parent_answers = parent.expected_rows()
    if (canonical_lines(parent_targets) != PARENT_TARGETS.read_text() or
            canonical_lines(parent_answers) != PARENT_ANSWERS.read_text()):
        raise ValueError('parent fixture replay differs')
    smallk = module('n83_final_online_smallk', parent.PARENT)
    final = smallk.module()
    precision = final.precision_module()
    group = precision.group
    prior = precision.prior_points()
    precision_targets, _ = precision.expected_rows()
    final_targets, _ = final.expected_rows()
    smallk_targets, _ = smallk.expected_rows()
    for set_of_points in (precision_targets, final_targets,
                          smallk_targets, parent_targets):
        prior.update(tuple(point) for point in set_of_points)
    targets, answers = [], []
    for index in range(COUNT):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, 'little')).digest()
        scalar = 1 + int.from_bytes(seed, 'little') % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if (tuple(encoded) in prior or not group.on_curve(point) or
                group.multiply(point, group.R) is not None):
            raise ValueError(f'holdout point {index} repeats or fails subgroup replay')
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
        'schema': 'n83.final-online-holdout-fixture/v1',
        'domain': DOMAIN.decode(), 'target_count': COUNT,
        'parent_generator_sha256': PARENT_SHA256,
        'parent_targets_sha256': PARENT_TARGETS_SHA256,
        'parent_answers_sha256': PARENT_ANSWERS_SHA256,
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
