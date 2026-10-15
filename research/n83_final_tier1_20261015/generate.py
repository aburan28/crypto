#!/usr/bin/env python3
"""Freeze five unseen same-point N83 IC/rho targets before finalist selection."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FINAL = ROOT / "research/n83_final_tournament_20261015/generate.py"
ANALYZER = ROOT / "research/n83_indexed5_20261013/analyze.py"
FINAL_TARGET_SHA256 = "cfa943ae538fcc41d2af8bb8fd1c8fc0f29adc340b55644e90382da25b213424"
ANALYZER_SHA256 = "e34d06b3275d86dcb3a659b79adf87dcdc078cebb6f510c4b5c42403d83a1804"
DOMAIN = b"n83-final-selected-tier1-fresh-5-20261015-v1"
COUNT = 5


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_once(path, rows):
    blob = "".join(json.dumps(row, separators=(",", ":")) + "\n" for row in rows)
    if path.exists() and path.read_text() != blob:
        raise ValueError(f"frozen fixture changed: {path}")
    if not path.exists():
        path.write_text(blob)
    return hashlib.sha256(blob.encode()).hexdigest()


def main():
    if digest(ANALYZER) != ANALYZER_SHA256:
        raise ValueError("independent N83 group implementation changed")
    fspec = importlib.util.spec_from_file_location("n83_final_targets", FINAL)
    final = importlib.util.module_from_spec(fspec)
    fspec.loader.exec_module(final)
    gspec = importlib.util.spec_from_file_location("n83_tier1_group", ANALYZER)
    group = importlib.util.module_from_spec(gspec)
    gspec.loader.exec_module(group)
    if not group.on_curve(group.G) or group.multiply(group.G, group.R) is not None:
        raise ValueError("frozen subgroup generator or order changed")
    prior = set()
    sources = dict(final.PRIOR_SHA256)
    sources["n83_final_tournament_20261015/targets.jsonl"] = FINAL_TARGET_SHA256
    for relative, expected in sources.items():
        path = ROOT / "research" / relative
        if digest(path) != expected:
            raise ValueError(f"prior point panel changed: {relative}")
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
            raise ValueError(f"fresh point {index} failed subgroup construction")
        encoded = [str(point[0]), str(point[1])]
        if tuple(encoded) in prior:
            raise ValueError(f"fresh point {index} repeats a prior public point")
        prior.add(tuple(encoded))
        points.append(encoded)
        answers.append({"fixture_index": index, "point": encoded,
                        "scalar": str(scalar), "seed_sha256": seed.hex()})
    print(json.dumps({"count": COUNT,
                      "targets_sha256": write_once(HERE / "targets.jsonl", points),
                      "answers_sha256": write_once(HERE / "answers.jsonl", answers)},
                     sort_keys=True))


if __name__ == "__main__":
    main()
