#!/usr/bin/env python3
"""Freeze five new N83 subgroup points for Bloom-on Tier 1 pairs."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ANALYZER = ROOT / "research/n83_indexed5_20261013/analyze.py"
ANALYZER_SHA = "e34d06b3275d86dcb3a659b79adf87dcdc078cebb6f510c4b5c42403d83a1804"
DOMAIN = b"n83-root11-bloom-tier1-fresh-20261015-v1"
PRIOR = (
    "research/n83_factor_base_100_20261014/targets.jsonl",
    "research/n83_factor_base_100_20261014/reserve_targets.jsonl",
    "research/n83_factor_base_refine_20261014/targets.jsonl",
    "research/n83_factor_base_refine_20261014/reserve_targets.jsonl",
    "research/n83_factor_base_confirm_20261015/targets.jsonl",
    "research/n83_indexed5_20261013/primary_target.jsonl",
)

spec = importlib.util.spec_from_file_location("n83_bloom_fresh_group", ANALYZER)
group = importlib.util.module_from_spec(spec)
spec.loader.exec_module(group)


def freeze(path, rows):
    blob = "".join(json.dumps(row, separators=(",", ":")) + "\n" for row in rows)
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"frozen fixture changed: {path}")
    if not path.exists():
        path.write_text(blob)
    return hashlib.sha256(blob.encode()).hexdigest()


def main():
    if hashlib.sha256(ANALYZER.read_bytes()).hexdigest() != ANALYZER_SHA:
        raise ValueError("independent group implementation changed")
    prior = set()
    for relative in PRIOR:
        path = ROOT / relative
        prior.update(tuple(json.loads(line)) for line in path.read_text().splitlines() if line)
    points, answers = [], []
    for index in range(5):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if (tuple(encoded) in prior or not group.on_curve(point)
                or group.multiply(point, group.R) is not None):
            raise ValueError(f"new target {index} is repeated or outside the subgroup")
        prior.add(tuple(encoded))
        points.append(encoded)
        answers.append({"fixture_index": index, "point": encoded,
                        "scalar": str(scalar), "seed_sha256": seed.hex()})
    point_sha = freeze(HERE / "targets.jsonl", points)
    answer_sha = freeze(HERE / "answers.jsonl", answers)
    print(json.dumps({"count": 5, "targets_sha256": point_sha,
                      "answers_sha256": answer_sha}, sort_keys=True))


if __name__ == "__main__":
    main()
