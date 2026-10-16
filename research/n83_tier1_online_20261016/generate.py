#!/usr/bin/env python3
"""Freeze five previously unused N83 points for a separate online panel."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
COLD = ROOT / "research/n83_cold_factor_base_20261016/generate.py"
EXPECTED = {
    "research/n83_cold_factor_base_20261016/generate.py":
        "1807dadb4dc423850f32742faa03beb5774dacb4d3c17618502b6debb623ca58",
    "research/n83_cold_factor_base_20261016/prior_manifest.json":
        "f0b29cb6e7d0081711ccb4b5ff54455f06198a2a3562d4654fbe7e81d07b1c8b",
    "research/n83_cold_factor_base_20261016/targets.jsonl":
        "2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d",
    "research/n83_cold_factor_base_20261016/reserve_targets.jsonl":
        "6ab08db41c5e86c96162f538a0c45d53331c2f9f823a66eea9bf0b458dd3554d",
    "research/n83_cold_factor_base_20261016/answers.jsonl":
        "b71c5a8bfdae603931fa6c03fc6f6faa3b261288cf787ac682bf8a3e3afe262d",
    "research/n83_cold_factor_base_20261016/reserve_answers.jsonl":
        "b21d3a1d54c701255a079a3f1d88940f46563a2f857303b4e3bb98ac40483f39",
}
DOMAIN = b"n83-tier1-online-five-fresh-20261016-v1"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for relative, expected in EXPECTED.items():
    if sha256(ROOT / relative) != expected:
        raise ValueError(f"prior N83 source or target changed: {relative}")

spec = importlib.util.spec_from_file_location("n83_cold_fixture", COLD)
cold = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cold)
group = cold.group


def frozen_prior():
    prior = cold.prior_points()
    for relative in (
        "research/n83_cold_factor_base_20261016/targets.jsonl",
        "research/n83_cold_factor_base_20261016/reserve_targets.jsonl",
    ):
        for line in (ROOT / relative).read_text().splitlines():
            encoded = json.loads(line)
            point = group.point_from_json(encoded)
            if not group.on_curve(point) or group.multiply(point, group.R) is not None:
                raise ValueError(f"prior target fails subgroup replay: {relative}")
            prior.add(tuple(encoded))
    return prior


def main():
    if (not group.on_curve(group.G)
            or group.multiply(group.G, group.R) is not None):
        raise ValueError("N83 subgroup generator or order changed")
    prior = frozen_prior()
    targets = []
    answers = []
    for index in range(5):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if (tuple(encoded) in prior or not group.on_curve(point)
                or group.multiply(point, group.R) is not None):
            raise ValueError(f"fresh target {index} repeats a prior point or fails replay")
        prior.add(tuple(encoded))
        targets.append(encoded)
        answers.append({"fixture_index": index, "point": encoded,
                        "scalar": str(scalar), "seed_sha256": seed.hex()})
    result = {
        "target_count": len(targets),
        "targets_sha256": cold.previous.freeze(HERE / "targets.jsonl", targets),
        "answers_sha256": cold.previous.freeze(HERE / "answers.jsonl", answers),
        "excluded_source_hashes": EXPECTED,
        "domain": DOMAIN.decode(),
    }
    manifest = json.dumps(result, sort_keys=True, indent=2) + "\n"
    path = HERE / "fixture_manifest.json"
    if path.exists() and path.read_text() != manifest:
        raise ValueError("frozen Tier 1 fixture manifest changed")
    path.write_text(manifest)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
