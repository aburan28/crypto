#!/usr/bin/env python3
"""Freeze 4,000 new N83 points for the balanced K3300/K3350 comparison."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ANALYZER = ROOT / "research/n83_indexed5_20261013/analyze.py"
ANALYZER_SHA = "e34d06b3275d86dcb3a659b79adf87dcdc078cebb6f510c4b5c42403d83a1804"
DOMAIN = b"n83-k-abba-fresh-4000-20261015-v1"
EXCLUDED_BLOOM_SHA = "ef616c5254a9d883e962024b9f4bd7f794cfe5576c73a93b099ac5be6fe45668"
PRIOR_20000_SHA = "c80f5217350d3b89bbf252298295b66e81dea673c222ce36816c77c182e9cda7"
PRIOR_K3350_SHA = "fc8be4c2c946bdbcb78a6845c8e2be3d1a51728c2ed4d812f6546ea5b000055e"
PRIOR = (
    "research/n83_factor_base_100_20261014/targets.jsonl",
    "research/n83_factor_base_100_20261014/reserve_targets.jsonl",
    "research/n83_factor_base_refine_20261014/targets.jsonl",
    "research/n83_factor_base_refine_20261014/reserve_targets.jsonl",
    "research/n83_factor_base_confirm_20261015/targets.jsonl",
    "research/n83_indexed5_20261013/primary_target.jsonl",
    "research/n83_k3250_sparse_fresh_20261015/targets.jsonl",
    "research/n83_k3350_sparse_fresh_20261015/targets.jsonl",
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
    if not group.on_curve(group.G) or group.multiply(group.G, group.R) is not None:
        raise ValueError("frozen subgroup generator or order changed")
    prior = set()
    prior_panel = ROOT / "research/n83_k3250_sparse_fresh_20261015/targets.jsonl"
    if hashlib.sha256(prior_panel.read_bytes()).hexdigest() != PRIOR_20000_SHA:
        raise ValueError("prior K3250/K3300 target exclusion changed")
    k3350_panel = ROOT / "research/n83_k3350_sparse_fresh_20261015/targets.jsonl"
    if hashlib.sha256(k3350_panel.read_bytes()).hexdigest() != PRIOR_K3350_SHA:
        raise ValueError("prior K3300/K3350 target exclusion changed")
    for relative in PRIOR:
        path = ROOT / relative
        prior.update(tuple(json.loads(line)) for line in path.read_text().splitlines() if line)
    excluded = HERE / "excluded_bloom_targets.jsonl"
    if hashlib.sha256(excluded.read_bytes()).hexdigest() != EXCLUDED_BLOOM_SHA:
        raise ValueError("earlier Bloom Tier 1 reserve point exclusion changed")
    prior.update(tuple(json.loads(line)) for line in excluded.read_text().splitlines() if line)
    points, answers = [], []
    for index in range(4_000):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if tuple(encoded) in prior or not group.on_curve(point):
            raise ValueError(f"target {index} is repeated or outside the subgroup")
        prior.add(tuple(encoded))
        points.append(encoded)
        answers.append({"fixture_index": index, "point": encoded,
                        "scalar": str(scalar), "seed_sha256": seed.hex()})
    point_sha = freeze(HERE / "targets.jsonl", points)
    answer_sha = freeze(HERE / "answers.jsonl", answers)
    print(json.dumps({"count": len(points), "targets_sha256": point_sha,
                      "answers_sha256": answer_sha}, sort_keys=True))


if __name__ == "__main__":
    main()
