#!/usr/bin/env python3
"""Seal 32 fresh N83 public points for the root11 finalist comparison."""

import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE
DOMAIN = b'n83-root11-finalists-32-20261017-v1'
COUNT = 32


spec = importlib.util.spec_from_file_location(
    'n83_holdout_group', ROOT / 'research/n83_indexed5_20261013/analyze.py')
group = importlib.util.module_from_spec(spec)
spec.loader.exec_module(group)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_once(path, blob):
    if path.exists() and path.read_bytes() != blob:
        raise ValueError(f'sealed fixture differs: {path}')
    path.write_bytes(blob)


def main():
    OUT.mkdir(exist_ok=True)
    targets, answers, points = [], [], []
    for index in range(COUNT):
        digest = hashlib.sha256(DOMAIN + index.to_bytes(4, 'little')).digest()
        scalar = 1 + int.from_bytes(digest, 'little') % (group.R - 1)
        point = group.multiply(group.G, scalar)
        if point is None or not group.on_curve(point):
            raise ValueError('point outside subgroup curve')
        points.append(point)
        targets.append(json.dumps([str(v) for v in point],
                                  separators=(',', ':')) + '\n')
        answers.append(json.dumps({'fixture_index': index,
                                   'scalar': str(scalar),
                                   'seed_sha256': digest.hex()},
                                  separators=(',', ':'), sort_keys=True) + '\n')
    if len(set(points)) != COUNT:
        raise ValueError('duplicate target point')
    target = OUT / 'targets.jsonl'
    answer = OUT / 'answers.jsonl'
    write_once(target, ''.join(targets).encode())
    write_once(answer, ''.join(answers).encode())
    manifest = {'schema': 'n83.root11-finalists-fixture/v1',
                'curve_id': 'EC1N83Ckb1a1haacb69ee9f37',
                'subgroup_order': group.R,
                'point_count': COUNT,
                'domain': DOMAIN.decode(),
                'scalar_rule':
                    '1+little_endian(SHA256(domain||u32_le(index))) mod (r-1)',
                'targets_sha256': sha256(target),
                'answers_sha256': sha256(answer),
                'generator_sha256': sha256(Path(__file__)),
                'status': 'SEALED_UNMEASURED'}
    write_once(OUT / 'manifest.json',
               (json.dumps(manifest, sort_keys=True, indent=2) + '\n').encode())
    print(json.dumps(manifest, sort_keys=True))


if __name__ == '__main__':
    main()
