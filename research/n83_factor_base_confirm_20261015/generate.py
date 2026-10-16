#!/usr/bin/env python3
"""Freeze 4000 disjoint N83 public targets before the confirmation run."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ANALYZER = ROOT / 'research/n83_indexed5_20261013/analyze.py'
ANALYZER_SHA256 = 'e34d06b3275d86dcb3a659b79adf87dcdc078cebb6f510c4b5c42403d83a1804'
DOMAIN = b'n83-factor-base-confirm-4000-20261015-v1'
TARGETS_SHA256 = '451df8524b2c356dd4a523260620c4dff406975a122762b98ddddb7987b75c2c'
ANSWERS_SHA256 = '8e18c9c97ba9f57e5eb43f1b39b9e9575b5cfab87840c25c428dac150c868101'

spec = importlib.util.spec_from_file_location('n83_100_group', ANALYZER)
group = importlib.util.module_from_spec(spec)
spec.loader.exec_module(group)


def frozen(path, rows, expected):
    content = ''.join(json.dumps(row, separators=(',', ':')) + '\n' for row in rows)
    if hashlib.sha256(content.encode()).hexdigest() != expected:
        raise ValueError(f'{path.name}: expected digest differs')
    if path.exists():
        if path.read_text() != content:
            raise ValueError(f'{path.name}: frozen file differs')
    else:
        path.write_text(content)


def main():
    if hashlib.sha256(ANALYZER.read_bytes()).hexdigest() != ANALYZER_SHA256:
        raise ValueError('independent N83 group implementation changed')
    prior = set()
    for path in (
        ROOT / 'research/n83_evaluation_20261012/optimization/fixtures/holdout_5_targets.jsonl',
        ROOT / 'research/n83_factor_base_100_20261014/excluded_final_targets.jsonl',
        ROOT / 'research/n83_factor_base_100_20261014/targets.jsonl',
        ROOT / 'research/n83_factor_base_100_20261014/reserve_targets.jsonl',
        ROOT / 'research/n83_indexed5_20261013/primary_target.jsonl',
        ROOT / 'research/n83_factor_base_refine_20261014/targets.jsonl',
        ROOT / 'research/n83_factor_base_refine_20261014/reserve_targets.jsonl',
    ):
        prior.update(tuple(json.loads(line)) for line in path.read_text().splitlines() if line)
    targets, answers = [], []
    for index in range(4000):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, 'little')).digest()
        scalar = 1 + int.from_bytes(seed, 'little') % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if (tuple(encoded) in prior or not group.on_curve(point)
                or group.multiply(point, group.R) is not None):
            raise ValueError(f'target {index}: duplicate or invalid subgroup point')
        prior.add(tuple(encoded))
        targets.append(encoded)
        answers.append({'fixture_index': index, 'point': encoded,
                        'scalar': str(scalar), 'seed_sha256': seed.hex()})
    frozen(HERE / 'targets.jsonl', targets, TARGETS_SHA256)
    frozen(HERE / 'answers.jsonl', answers, ANSWERS_SHA256)
    print(json.dumps({'count': len(targets),
                      'targets_sha256': hashlib.sha256((HERE / 'targets.jsonl').read_bytes()).hexdigest(),
                      'answers_sha256': hashlib.sha256((HERE / 'answers.jsonl').read_bytes()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
