#!/usr/bin/env python3
"""Freeze 1,000 unseen N83 points for the root4 wide3 Bloom comparison."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ANALYZER = ROOT / "research/n83_indexed5_20261013/analyze.py"
ANALYZER_SHA = "e34d06b3275d86dcb3a659b79adf87dcdc078cebb6f510c4b5c42403d83a1804"
DOMAIN = b"n83-root4-wide3-fresh-1000-20261016-v1"
MEMORY_FRONTIER_DOMAIN = b"n83-root4-memory-frontier-fresh-1000-20261015-v1"
MEMORY_FRONTIER_SHA256 = "0626a709fc285817b63c1a6a46310cfe3577824ff7a02ba5f0c9d6f2a7b120c3"
EXCLUDED_BLOOM_SHA = "ef616c5254a9d883e962024b9f4bd7f794cfe5576c73a93b099ac5be6fe45668"
PRIOR_20000_SHA = "c80f5217350d3b89bbf252298295b66e81dea673c222ce36816c77c182e9cda7"
PRIOR_K3350_SHA = "fc8be4c2c946bdbcb78a6845c8e2be3d1a51728c2ed4d812f6546ea5b000055e"
PRIOR_ABBA_SHA = "35e8b43966d15065cf0186ace762f823a602aa1be7d7c0a9039e020550d01491"
PRIOR_HIGHK_SHA = "48a243619e270abc096e172fc669bb03494fd8c09d4e39c55935121408d44eea"
PRIOR_BEST_E2E_SHA = "bbb3a2fffdcb51db6105d00f12f9dc58450520a5bae646bd3dfb08d8dc72eab4"
PRIOR_LEGACY_SHA = "1603eea50d7ddebea7fe23d6c3ba7f480736f2682bcba819f5592daa8081c884"
PRIOR = (
    "research/n83_factor_base_100_20261014/targets.jsonl",
    "research/n83_factor_base_100_20261014/reserve_targets.jsonl",
    "research/n83_factor_base_refine_20261014/targets.jsonl",
    "research/n83_factor_base_refine_20261014/reserve_targets.jsonl",
    "research/n83_factor_base_confirm_20261015/targets.jsonl",
    "research/n83_indexed5_20261013/primary_target.jsonl",
    "research/n83_k3250_sparse_fresh_20261015/targets.jsonl",
    "research/n83_k3350_sparse_fresh_20261015/targets.jsonl",
    "research/n83_k_abba_20261015/targets.jsonl",
    "research/n83_root4_highk_20261015/targets.jsonl",
    "research/n83_best_e2e_20261015/targets.jsonl",
    "research/n83_legacy_challenger_20261015/targets.jsonl",
)

MORE_PRIOR = {
    "research/n83_root4_bloom_20261015/targets.jsonl": "71b93fd3a397ee2c1204dc14fc451894eee7537b99f50d75c98d65dee6b1e3d1",
    "research/n83_root11_highk_20261015/targets.jsonl": "fec0ad0c82b0608bbc7fdf4b9b29ba47e2c4761124e80ac7827d7e4260e78b42",
    "research/n83_root11_edgek_20261015/targets.jsonl": "c75cbb92ef1a23057272cf3481cc649b12495135c75affd1c2610205a42b7aef",
    "research/n83_final_tournament_20261015/targets.jsonl": "cfa943ae538fcc41d2af8bb8fd1c8fc0f29adc340b55644e90382da25b213424",
    "research/n83_edge_confirmation_20261015/targets.jsonl": "66c4683e252f1cd321b4b6c7c9d06744c82140222ddce160faa4e7a0e0984b7e",
    "research/n83_final_tier1_20261015/targets.jsonl": "34c530851eb135b7e348ac976caaae37c5883c52b52962ad496028c058b1dcda",
    "research/n83_unresolved_tier1_20261015/targets.jsonl": "dfea95ee4e423dd627adbd116828e9f478779c51476f0fa172cf4c160e645d12",
}

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
    abba_panel = ROOT / "research/n83_k_abba_20261015/targets.jsonl"
    if hashlib.sha256(abba_panel.read_bytes()).hexdigest() != PRIOR_ABBA_SHA:
        raise ValueError("prior balanced K target exclusion changed")
    highk_panel = ROOT / "research/n83_root4_highk_20261015/targets.jsonl"
    if hashlib.sha256(highk_panel.read_bytes()).hexdigest() != PRIOR_HIGHK_SHA:
        raise ValueError("prior high-K target exclusion changed")
    best_e2e_panel = ROOT / "research/n83_best_e2e_20261015/targets.jsonl"
    if hashlib.sha256(best_e2e_panel.read_bytes()).hexdigest() != PRIOR_BEST_E2E_SHA:
        raise ValueError("prior six-arm target exclusion changed")
    legacy_panel = ROOT / "research/n83_legacy_challenger_20261015/targets.jsonl"
    if hashlib.sha256(legacy_panel.read_bytes()).hexdigest() != PRIOR_LEGACY_SHA:
        raise ValueError("prior root22/root21 challenger exclusion changed")
    for relative in PRIOR:
        path = ROOT / relative
        prior.update(tuple(json.loads(line)) for line in path.read_text().splitlines() if line)
    excluded = ROOT / "research/n83_best_e2e_20261015/excluded_bloom_targets.jsonl"
    if hashlib.sha256(excluded.read_bytes()).hexdigest() != EXCLUDED_BLOOM_SHA:
        raise ValueError("earlier Bloom Tier 1 reserve point exclusion changed")
    prior.update(tuple(json.loads(line)) for line in excluded.read_text().splitlines() if line)
    for relative, expected_sha in MORE_PRIOR.items():
        path = ROOT / relative
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected_sha:
            raise ValueError(f"prior target exclusion changed: {relative}")
        prior.update(tuple(json.loads(line)) for line in path.read_text().splitlines() if line)
    # The memory-frontier source is still on a separate frozen branch. Recreate
    # its public points from its pinned input law to keep this panel disjoint.
    memory_points = []
    for index in range(1_000):
        seed = hashlib.sha256(MEMORY_FRONTIER_DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = (str(point[0]), str(point[1]))
        if encoded in prior:
            raise ValueError(f"memory-frontier target {index} duplicates earlier panel")
        prior.add(encoded)
        memory_points.append(list(encoded))
    memory_blob = "".join(json.dumps(point, separators=(",", ":")) + "\n"
                          for point in memory_points).encode()
    if hashlib.sha256(memory_blob).hexdigest() != MEMORY_FRONTIER_SHA256:
        raise ValueError("concurrent memory-frontier public points changed")
    points, answers = [], []
    for index in range(1_000):
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
