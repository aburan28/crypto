#!/usr/bin/env python3
"""Freeze public N83 points for the post-screen, same-host tournament."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ANALYZER = ROOT / "research/n83_indexed5_20261013/analyze.py"
ANALYZER_SHA256 = "e34d06b3275d86dcb3a659b79adf87dcdc078cebb6f510c4b5c42403d83a1804"
DOMAIN = b"n83-final-tournament-fresh-4000-20261015-v1"
COUNT = 4_000
PRIOR_SHA256 = {
    "n83_factor_base_100_20261014/targets.jsonl": "fa0044f7cad85185812ac4c79f4feabb300ddcdc1732b6a0890778d5b8e72af6",
    "n83_factor_base_100_20261014/reserve_targets.jsonl": "644c56a7bae2c92085704bb6f58ab70b38b278710de90152ece7631058c41dc9",
    "n83_factor_base_refine_20261014/targets.jsonl": "35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8",
    "n83_factor_base_refine_20261014/reserve_targets.jsonl": "2ccaac2baf367598e0d4e724b500ba66c5b4a711f0be8207565ce3acd6b1d0b4",
    "n83_factor_base_confirm_20261015/targets.jsonl": "451df8524b2c356dd4a523260620c4dff406975a122762b98ddddb7987b75c2c",
    "n83_indexed5_20261013/primary_target.jsonl": "21c8dfc7d488909e5c4a42c839fcc01bdefa9821c11ba001626a1d15697240e8",
    "n83_k3250_sparse_fresh_20261015/targets.jsonl": "c80f5217350d3b89bbf252298295b66e81dea673c222ce36816c77c182e9cda7",
    "n83_k3350_sparse_fresh_20261015/targets.jsonl": "fc8be4c2c946bdbcb78a6845c8e2be3d1a51728c2ed4d812f6546ea5b000055e",
    "n83_k_abba_20261015/targets.jsonl": "35e8b43966d15065cf0186ace762f823a602aa1be7d7c0a9039e020550d01491",
    "n83_root4_highk_20261015/targets.jsonl": "48a243619e270abc096e172fc669bb03494fd8c09d4e39c55935121408d44eea",
    "n83_best_e2e_20261015/targets.jsonl": "bbb3a2fffdcb51db6105d00f12f9dc58450520a5bae646bd3dfb08d8dc72eab4",
    "n83_legacy_challenger_20261015/targets.jsonl": "1603eea50d7ddebea7fe23d6c3ba7f480736f2682bcba819f5592daa8081c884",
    "n83_root4_bloom_20261015/targets.jsonl": "71b93fd3a397ee2c1204dc14fc451894eee7537b99f50d75c98d65dee6b1e3d1",
    "n83_root11_highk_20261015/targets.jsonl": "fec0ad0c82b0608bbc7fdf4b9b29ba47e2c4761124e80ac7827d7e4260e78b42",
    "n83_best_e2e_20261015/excluded_bloom_targets.jsonl": "ef616c5254a9d883e962024b9f4bd7f794cfe5576c73a93b099ac5be6fe45668",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def freeze(path, rows):
    blob = "".join(json.dumps(row, separators=(",", ":")) + "\n" for row in rows)
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"frozen fixture changed: {path}")
    if not path.exists():
        path.write_text(blob)
    return hashlib.sha256(blob.encode()).hexdigest()


def main():
    if digest(ANALYZER) != ANALYZER_SHA256:
        raise ValueError("independent group implementation changed")
    spec = importlib.util.spec_from_file_location("n83_final_group", ANALYZER)
    group = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(group)
    if not group.on_curve(group.G) or group.multiply(group.G, group.R) is not None:
        raise ValueError("frozen subgroup generator or order changed")
    prior = set()
    for relative, expected in PRIOR_SHA256.items():
        path = ROOT / "research" / relative
        if digest(path) != expected:
            raise ValueError(f"prior target exclusion changed: {relative}")
        for line in path.read_text().splitlines():
            row = json.loads(line)
            if isinstance(row, list) and len(row) == 2:
                prior.add(tuple(row))
            elif isinstance(row, dict):
                point = row.get("target") or row.get("point")
                if isinstance(point, list) and len(point) == 2:
                    prior.add(tuple(point))
    points, answers = [], []
    for index in range(COUNT):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        if point is None or not group.on_curve(point):
            raise ValueError(f"target {index} is outside the subgroup")
        encoded = [str(point[0]), str(point[1])]
        if tuple(encoded) in prior:
            raise ValueError(f"target {index} is repeated")
        prior.add(tuple(encoded))
        points.append(encoded)
        answers.append({"fixture_index": index, "point": encoded,
                        "scalar": str(scalar), "seed_sha256": seed.hex()})
    print(json.dumps({"count": COUNT,
                      "targets_sha256": freeze(HERE / "targets.jsonl", points),
                      "answers_sha256": freeze(HERE / "answers.jsonl", answers)},
                     sort_keys=True))


if __name__ == "__main__":
    main()
