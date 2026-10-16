#!/usr/bin/env python3
"""Freeze five N83 cold-factor-base holdout points before the K-grid result."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
COLD = ROOT / "research/n83_cold_factor_base_20261016/generate.py"
PRIOR_FILES = {
    "research/n83_cold_factor_base_20261016/targets.jsonl":
        "2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d",
    "research/n83_cold_factor_base_20261016/reserve_targets.jsonl":
        "6ab08db41c5e86c96162f538a0c45d53331c2f9f823a66eea9bf0b458dd3554d",
    "research/n83_tier1_online_20261016/targets.jsonl":
        "4976679ac6c9e3858cd26312ec2e50b8eb30971041ea38286c9582dd9217a1c1",
}
DOMAIN = b"n83-e2e-cold-holdout-five-20261016-v1"

spec = importlib.util.spec_from_file_location("n83_cold_fixture", COLD)
cold = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cold)
group = cold.group


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    prior = cold.prior_points()
    for relative, expected in PRIOR_FILES.items():
        path = ROOT / relative
        if sha256(path) != expected:
            raise ValueError(f"prior target file changed: {relative}")
        for line in path.read_text().splitlines():
            encoded = json.loads(line)
            point = group.point_from_json(encoded)
            if not group.on_curve(point) or group.multiply(point, group.R) is not None:
                raise ValueError(f"prior point fails subgroup replay: {relative}")
            prior.add(tuple(encoded))
    if (not group.on_curve(group.G)
            or group.multiply(group.G, group.R) is not None):
        raise ValueError("N83 subgroup generator or order changed")
    targets, answers = [], []
    for index in range(5):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if (tuple(encoded) in prior or not group.on_curve(point)
                or group.multiply(point, group.R) is not None):
            raise ValueError(f"holdout point {index} repeats or fails replay")
        prior.add(tuple(encoded))
        targets.append(encoded)
        answers.append({"fixture_index": index, "point": encoded,
                        "scalar": str(scalar), "seed_sha256": seed.hex()})
    manifest = {
        "domain": DOMAIN.decode(), "target_count": len(targets),
        "prior_target_hashes": PRIOR_FILES,
        "targets_sha256": cold.previous.freeze(HERE / "targets.jsonl", targets),
        "answers_sha256": cold.previous.freeze(HERE / "answers.jsonl", answers),
    }
    blob = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    path = HERE / "fixture_manifest.json"
    if path.exists() and path.read_text() != blob:
        raise ValueError("holdout manifest differs")
    path.write_text(blob)
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    main()
