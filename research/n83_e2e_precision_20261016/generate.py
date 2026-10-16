#!/usr/bin/env python3
"""Freeze and independently replay sixteen new N83 public points."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HOLDOUT = ROOT / "research/n83_e2e_holdout_20261016/generate.py"
DOMAIN = b"n83-e2e-cold-precision-sixteen-20261016-v1"
COUNT = 16
PRIOR_FILES = {
    "research/n83_cold_factor_base_20261016/targets.jsonl":
        "2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d",
    "research/n83_cold_factor_base_20261016/reserve_targets.jsonl":
        "6ab08db41c5e86c96162f538a0c45d53331c2f9f823a66eea9bf0b458dd3554d",
    "research/n83_tier1_online_20261016/targets.jsonl":
        "4976679ac6c9e3858cd26312ec2e50b8eb30971041ea38286c9582dd9217a1c1",
    "research/n83_e2e_holdout_20261016/targets.jsonl":
        "722bc60cd010e533cd715a9bfb321fba52d950352c99684cd422b23af7f93d26",
}

spec = importlib.util.spec_from_file_location("n83_holdout_fixture", HOLDOUT)
holdout = importlib.util.module_from_spec(spec)
spec.loader.exec_module(holdout)
group = holdout.group


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_points():
    points = set()
    for relative, expected in PRIOR_FILES.items():
        path = ROOT / relative
        if digest(path) != expected:
            raise ValueError(f"prior N83 point file differs: {relative}")
        for line in path.read_text().splitlines():
            encoded = json.loads(line)
            point = group.point_from_json(encoded)
            if not group.on_curve(point) or group.multiply(point, group.R) is not None:
                raise ValueError(f"prior point fails subgroup replay: {relative}")
            points.add(tuple(encoded))
    return points


def expected_rows():
    prior = prior_points()
    if not group.on_curve(group.G) or group.multiply(group.G, group.R) is not None:
        raise ValueError("N83 subgroup generator or order differs")
    targets, answers = [], []
    for index in range(COUNT):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if (tuple(encoded) in prior or not group.on_curve(point)
                or group.multiply(point, group.R) is not None):
            raise ValueError(f"precision point {index} repeats or fails replay")
        prior.add(tuple(encoded))
        targets.append(encoded)
        answers.append({"fixture_index": index, "point": encoded,
                        "scalar": str(scalar), "seed_sha256": seed.hex()})
    return targets, answers


def canonical_lines(rows):
    return "".join(json.dumps(row, separators=(",", ":")) + "\n" for row in rows)


def main(verify):
    targets, answers = expected_rows()
    target_blob, answer_blob = canonical_lines(targets), canonical_lines(answers)
    manifest = {
        "schema": "n83.e2e-precision-fixture/v1",
        "domain": DOMAIN.decode(), "target_count": COUNT,
        "prior_target_hashes": PRIOR_FILES,
        "targets_sha256": hashlib.sha256(target_blob.encode()).hexdigest(),
        "answers_sha256": hashlib.sha256(answer_blob.encode()).hexdigest(),
    }
    files = {HERE / "targets.jsonl": target_blob,
             HERE / "answers.jsonl": answer_blob,
             HERE / "fixture_manifest.json": json.dumps(manifest, indent=2,
                                                         sort_keys=True) + "\n"}
    for path, blob in files.items():
        if path.exists():
            if path.read_text() != blob:
                raise ValueError(f"frozen precision fixture differs: {path.name}")
        elif verify:
            raise ValueError(f"missing committed precision fixture: {path.name}")
        else:
            path.write_text(blob)
    print(json.dumps({"mode": "sealed_replay" if verify else "freeze",
                      **manifest}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-sealed", action="store_true")
    args = parser.parse_args()
    main(args.verify_sealed)
