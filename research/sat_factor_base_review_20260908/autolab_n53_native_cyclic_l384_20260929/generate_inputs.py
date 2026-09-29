#!/usr/bin/env python3
"""Rebuild the frozen A/B training and public-Q schedules without n53 producers."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SHARED = ROOT / "autolab_shared_log_n53_20260925"
ROTATION = ROOT / "autolab_n53_rank_rotation_20260925"
OLD_Q = ROOT / "autolab_combined_l384_certified_coldbase_n53_20260925/points_L384.jsonl"
DOMAIN_B = b"ECC2K53-NATIVE-CYCLIC-L384-20260929-TRAIN-B-v1/"
DOMAIN_Q = b"ECC2K53-NATIVE-CYCLIC-L384-20260929-Q-v1/"


def load_math():
    spec = importlib.util.spec_from_file_location("native_l384_independent_math", SHARED / "generate_targets.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def scalar_bytes(values: list[int]) -> bytes:
    return b"".join(f"{value}\n".encode() for value in values)


def point_bytes(points: list[tuple[int, int]]) -> bytes:
    return b"".join((json.dumps(point, separators=(",", ":")) + "\n").encode() for point in points)


def draw(math, domain: bytes, count: int, excluded: set[int], old_q_points: set[tuple[int, int]]):
    values: list[int] = []
    points: list[tuple[int, int]] = []
    counter = 0
    while len(values) < count:
        candidate = 1 + int.from_bytes(hashlib.sha256(domain + str(counter).encode()).digest(), "big") % (math.ORDER - 1)
        counter += 1
        if candidate in excluded:
            continue
        point = math.scalar(math.GENERATOR, candidate)
        assert point is not None and math.on_curve(point)
        if point in old_q_points:
            continue
        excluded.add(candidate)
        values.append(candidate)
        points.append(point)
    assert len(set(values)) == len(set(points)) == count
    return values, points, counter


def generate():
    math = load_math()
    a_values = [int(line) for line in (ROTATION / "target_scalars.txt").read_text().splitlines()]
    a_points = [tuple(json.loads(line)) for line in (ROTATION / "target_points.jsonl").read_text().splitlines()]
    assert len(a_values) == len(set(a_values)) == len(a_points) == len(set(a_points)) == 512
    assert all(math.scalar(math.GENERATOR, scalar) == point for scalar, point in zip(a_values, a_points))
    old_train = math.training_scalars()
    old_manifest = json.loads((SHARED / "TARGET_MANIFEST.json").read_text())
    old_q_values = {scalar for block in old_manifest["blocks"] for scalar in block["validator_scalars"]}
    old_q_points = {tuple(json.loads(line)) for line in OLD_Q.read_text().splitlines()}
    assert len(old_train) == 512 and len(old_q_values) == len(old_q_points) == 384
    assert not set(a_values) & (old_train | old_q_values)
    assert not set(a_points) & old_q_points
    excluded = old_train | old_q_values | set(a_values)
    b_values, b_points, b_counter = draw(math, DOMAIN_B, 512, excluded, old_q_points)
    q_values, q_points, q_counter = draw(math, DOMAIN_Q, 384, excluded, old_q_points)
    assert not (set(a_values) & set(b_values) or set(a_values) & set(q_values)
                or set(b_values) & set(q_values))
    assert not (set(a_points) & set(b_points) or set(a_points) & set(q_points)
                or set(b_points) & set(q_points))
    files = {
        "training_A_scalars.txt": scalar_bytes(a_values),
        "training_B_scalars.txt": scalar_bytes(b_values),
        "training_B_points.jsonl": point_bytes(b_points),
        "points_L384.jsonl": point_bytes(q_points),
        "validator_scalars_L384.txt": scalar_bytes(q_values),
    }
    assert files["training_A_scalars.txt"] == (ROTATION / "target_scalars.txt").read_bytes()
    report = {
        "schema": "native_cyclic_l384_input_proof_v1",
        "domains": {"B": DOMAIN_B.decode(), "Q": DOMAIN_Q.decode()},
        "counters_consumed": {"B": b_counter, "Q": q_counter},
        "counts": {"A": 512, "B": 512, "Q": 384, "old_training": 512, "old_Q": 384},
        "overlaps": {"A_B": 0, "A_Q": 0, "B_Q": 0,
                     "A_old_training": 0, "B_old_training": 0, "Q_old_training": 0,
                     "A_old_Q": 0, "B_old_Q": 0, "Q_old_Q": 0},
        "sha256": {name: sha(data) for name, data in files.items()},
    }
    return files, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    files, report = generate()
    for name, data in files.items():
        path = HERE / name
        if args.write:
            assert not path.exists(), f"refuse overwrite: {path}"
            path.write_bytes(data)
        else:
            assert path.read_bytes() == data, f"input drift: {name}"
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
