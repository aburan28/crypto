#!/usr/bin/env python3
"""Freeze a new n53 scalar/point corpus disjoint from both earlier corpora."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parent / "autolab_shared_log_n53_20260925"
OLD_POINTS = HERE.parent / "autolab_combined_l384_coldbase_n53_20260925/points_L384.jsonl"
MATH_PATH = SHARED / "generate_targets.py"
MANIFEST = SHARED / "TARGET_MANIFEST.json"
DOMAIN = b"ECC2K53-TARGET-CYCLIC-RANK-20260925-v1/"
COUNT = 512


def load_math():
    spec = importlib.util.spec_from_file_location("n53_independent_group", MATH_PATH)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def generate() -> tuple[bytes, bytes, dict]:
    math = load_math()
    manifest = json.loads(MANIFEST.read_text())
    old_training = math.training_scalars()
    old_q_scalars = {
        scalar for block in manifest["blocks"] for scalar in block["validator_scalars"]
    }
    old_q_points = {tuple(json.loads(line)) for line in OLD_POINTS.read_text().splitlines()}
    assert len(old_training) == 512
    assert len(old_q_scalars) == len(old_q_points) == 384
    assert not old_training & old_q_scalars
    excluded = old_training | old_q_scalars
    values = []
    points = []
    seen = set()
    counter = 0
    while len(values) < COUNT:
        candidate = 1 + int.from_bytes(
            hashlib.sha256(DOMAIN + str(counter).encode()).digest(), "big"
        ) % (math.ORDER - 1)
        counter += 1
        if candidate in excluded or candidate in seen:
            continue
        point = math.scalar(math.GENERATOR, candidate)
        assert point is not None and math.on_curve(point)
        if point in old_q_points:
            continue
        seen.add(candidate)
        values.append(candidate)
        points.append(point)
    assert len(set(points)) == len(points) == COUNT
    assert not set(points) & old_q_points
    assert not set(values) & old_training
    scalar_bytes = b"".join(f"{value}\n".encode() for value in values)
    point_bytes = b"".join((json.dumps(point, separators=(",", ":")) + "\n").encode() for point in points)
    report = {
        "domain": DOMAIN.decode(),
        "target_count": COUNT,
        "counters_consumed": counter,
        "old_training_count": len(old_training),
        "old_public_q_count": len(old_q_points),
        "old_training_scalar_overlap": len(set(values) & old_training),
        "old_public_q_scalar_overlap": len(set(values) & old_q_scalars),
        "old_public_q_point_overlap": len(set(points) & old_q_points),
        "scalar_sha256": hashlib.sha256(scalar_bytes).hexdigest(),
        "point_sha256": hashlib.sha256(point_bytes).hexdigest(),
    }
    return scalar_bytes, point_bytes, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    scalar_bytes, point_bytes, report = generate()
    for path, data in ((HERE / "target_scalars.txt", scalar_bytes),
                       (HERE / "target_points.jsonl", point_bytes)):
        if args.write:
            assert not path.exists(), f"refuse to overwrite {path}"
            path.write_bytes(data)
        else:
            assert path.read_bytes() == data, f"target schedule drift: {path}"
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
