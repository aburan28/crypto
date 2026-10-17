#!/usr/bin/env python3
"""Freeze ten N83 public points disjoint from earlier benchmark fixtures."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / 'research/n83_k3072_current_source_20261017/generate.py'
PARENT_SHA256 = 'b7211a3beaf6e349443c7708746043d780782a47f0959cb95560a85666a890d4'
PRIOR_TARGETS = PARENT.parent / 'targets.jsonl'
PRIOR_ANSWERS = PARENT.parent / 'answers.jsonl'
PRIOR_TARGETS_SHA256 = '9bf9bfca72d4e82a5780759d52ac000b812ca25e831af04238983c8b72d2d777'
PRIOR_ANSWERS_SHA256 = '3f7a75dc0134f883a77fa2b87eb6839b6df804c53950d8dc77e7f29b2500284c'
DOMAIN = b'n83-highk-tiebreak-ten-20261017-v1'
COUNT = 10


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
    if digest(PARENT) != PARENT_SHA256:
        raise ValueError('parent fixture generator differs')
    if (digest(PRIOR_TARGETS) != PRIOR_TARGETS_SHA256 or
            digest(PRIOR_ANSWERS) != PRIOR_ANSWERS_SHA256):
        raise ValueError('prior comparison fixture differs')
    parent = module('n83_highk_parent_generator', PARENT)
    prior_targets, prior_answers = parent.expected_rows()
    if (canonical_lines(prior_targets) != PRIOR_TARGETS.read_text() or
            canonical_lines(prior_answers) != PRIOR_ANSWERS.read_text()):
        raise ValueError('prior comparison fixture replay differs')
    precision = parent.module().precision_module()
    group = precision.group
    used = precision.prior_points()
    precision_targets, _ = precision.expected_rows()
    used.update(tuple(item) for item in precision_targets)
    finalist_targets, _ = parent.module().expected_rows()
    used.update(tuple(item) for item in finalist_targets)
    used.update(tuple(item) for item in prior_targets)
    targets, answers = [], []
    for index in range(COUNT):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, 'little')).digest()
        scalar = 1 + int.from_bytes(seed, 'little') % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if (tuple(encoded) in used or not group.on_curve(point) or
                group.multiply(point, group.R) is not None):
            raise ValueError(f'fresh point {index} repeats or fails subgroup replay')
        used.add(tuple(encoded))
        targets.append(encoded)
        answers.append({'fixture_index': index, 'point': encoded,
                        'scalar': str(scalar), 'seed_sha256': seed.hex()})
    return targets, answers


def main(verify):
    targets, answers = expected_rows()
    target_blob = canonical_lines(targets)
    answer_blob = canonical_lines(answers)
    manifest = {
        'schema': 'n83.highk-tiebreak-fixture/v1',
        'domain': DOMAIN.decode(), 'target_count': COUNT,
        'parent_generator_sha256': digest(PARENT),
        'prior_targets_sha256': PRIOR_TARGETS_SHA256,
        'prior_answers_sha256': PRIOR_ANSWERS_SHA256,
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
