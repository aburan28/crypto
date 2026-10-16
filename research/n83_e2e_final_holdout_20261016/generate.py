#!/usr/bin/env python3
"""Freeze sixteen N83 public points before the finalist is selected."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "research/n83_e2e_precision_20261016/generate.py"
PARENT_SHA256 = "3f2ddb2b0f56c74bc030a32a138cf6230310d050cd7d200f13f9782aa679b4b9"
PRIOR_TARGETS = ROOT / "research/n83_e2e_precision_20261016/targets.jsonl"
PRIOR_TARGETS_SHA256 = "266baaec03569f6e6502a50c20829b9fc0b81b185abc2aceb507700904b52b42"
PRIOR_ANSWERS = ROOT / "research/n83_e2e_precision_20261016/answers.jsonl"
PRIOR_ANSWERS_SHA256 = "3b1af7adec31bf394230551c3c4ffea9af01951dd911362b3f1a0eea7c694ee1"
DOMAIN = b"n83-e2e-final-holdout-sixteen-20261016-v1"
COUNT = 16


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def precision_module():
    if digest(PARENT) != PARENT_SHA256:
        raise ValueError("precision fixture generator differs")
    spec = importlib.util.spec_from_file_location("n83_final_prior_fixture", PARENT)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def expected_rows():
    parent = precision_module()
    group = parent.group
    if (digest(PRIOR_TARGETS) != PRIOR_TARGETS_SHA256
            or digest(PRIOR_ANSWERS) != PRIOR_ANSWERS_SHA256):
        raise ValueError("sealed precision target or answer file differs")
    prior = parent.prior_points()
    precision_points, precision_answers = parent.expected_rows()
    if (parent.canonical_lines(precision_points) != PRIOR_TARGETS.read_text()
            or parent.canonical_lines(precision_answers) != PRIOR_ANSWERS.read_text()):
        raise ValueError("precision target or answer replay differs")
    for encoded in precision_points:
        prior.add(tuple(encoded))
    targets, answers = [], []
    for index in range(COUNT):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if (tuple(encoded) in prior or not group.on_curve(point)
                or group.multiply(point, group.R) is not None):
            raise ValueError(f"final holdout point {index} repeats or fails replay")
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
        "schema": "n83.e2e-final-holdout-fixture/v1",
        "domain": DOMAIN.decode(), "target_count": COUNT,
        "parent_generator_sha256": PARENT_SHA256,
        "prior_targets_sha256": PRIOR_TARGETS_SHA256,
        "prior_answers_sha256": PRIOR_ANSWERS_SHA256,
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
                raise ValueError(f"frozen finalist fixture differs: {path.name}")
        elif verify:
            raise ValueError(f"missing committed finalist fixture: {path.name}")
        else:
            path.write_text(blob)
    print(json.dumps({"mode": "sealed_replay" if verify else "freeze",
                      **manifest}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-sealed", action="store_true")
    args = parser.parse_args()
    main(args.verify_sealed)
