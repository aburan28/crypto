#!/usr/bin/env python3
"""Freeze one N83 cold-sweep point and four untouched confirmation points."""

import hashlib
import importlib.util
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREVIOUS = ROOT / "research/n83_compact4_crosslayout_20261016/generate.py"
PREVIOUS_SHA256 = "638624dc0d9ed678f1c0a716c31701c8ed81eba328b139916cf222bc3ce43df6"
ADMISSION = ROOT / "research/n83_compact4_crosslayout_results_20261016/selection.json"
ADMISSION_SHA256 = "75a8f9f7638839d1d3d37e91576a490597f97fb0dfa5cab8b7526d38176415be"
CURRENT_TARGETS = ROOT / "research/n83_compact4_crosslayout_20261016/targets.jsonl"
CURRENT_TARGETS_SHA256 = "6c5558613d4f1115dc7bdb09b74d569645c0e85328ddad6555a2484317803a08"
DOMAIN = b"n83-cold-factor-base-one-plus-four-20261016-v1"
PRIOR_MANIFEST = HERE / "prior_manifest.json"
PRIOR_ROOT = Path(os.environ.get("N83_COLD_PRIOR_ROOT", str(ROOT))).resolve()

if hashlib.sha256(PREVIOUS.read_bytes()).hexdigest() != PREVIOUS_SHA256:
    raise ValueError("previous target generator changed")
spec = importlib.util.spec_from_file_location("n83_previous_target_generator", PREVIOUS)
previous = importlib.util.module_from_spec(spec)
spec.loader.exec_module(previous)
group = previous.group


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_path(relative):
    local = ROOT / relative
    if local.exists():
        return local
    alternate = PRIOR_ROOT / relative
    if alternate.exists():
        return alternate
    raise FileNotFoundError(f"generate the pinned prior fixture before this sweep: {relative}")


def prior_points():
    if sha256(ADMISSION) != ADMISSION_SHA256 or sha256(CURRENT_TARGETS) != CURRENT_TARGETS_SHA256:
        raise ValueError("admission panel or its public points changed")
    decision = json.loads(ADMISSION.read_text())
    if (decision["status"] != "BOUNDED_UNKNOWN"
            or decision["descriptive_lowest_pipeline"] != "root11k3450wide3"):
        raise ValueError("cold-sweep admission result changed")
    files = set(previous.PRIOR) | set(previous.MORE_PRIOR)
    files.update((
        "research/n83_best_e2e_20261015/excluded_bloom_targets.jsonl",
        "research/n83_compact4_crosslayout_20261016/targets.jsonl",
    ))
    current = {relative: sha256(prior_path(relative)) for relative in sorted(files)}
    if PRIOR_MANIFEST.exists():
        if json.loads(PRIOR_MANIFEST.read_text()) != current:
            raise ValueError("prior target manifest changed")
    else:
        PRIOR_MANIFEST.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n")
    points = set()
    for relative in current:
        points.update(tuple(json.loads(line)) for line in
                      prior_path(relative).read_text().splitlines() if line)
    wide3 = []
    for index in range(1_000):
        seed = hashlib.sha256(previous.WIDE3_DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        points.add(tuple(encoded))
        wide3.append(encoded)
    blob = "".join(json.dumps(point, separators=(",", ":")) + "\n" for point in wide3)
    if hashlib.sha256(blob.encode()).hexdigest() != previous.WIDE3_POINT_SHA256:
        raise ValueError("concurrent wide3 target exclusion changed")
    return points


def main():
    if (sha256(previous.ANALYZER) != previous.ANALYZER_SHA
            or not group.on_curve(group.G)
            or group.multiply(group.G, group.R) is not None):
        raise ValueError("independent N83 group implementation or subgroup changed")
    prior = prior_points()
    points, answers = [], []
    for index in range(5):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        encoded = [str(point[0]), str(point[1])]
        if tuple(encoded) in prior or not group.on_curve(point):
            raise ValueError(f"point {index} repeats a previous panel or leaves the curve")
        prior.add(tuple(encoded))
        points.append(encoded)
        answers.append({"fixture_index": index, "point": encoded,
                        "scalar": str(scalar), "seed_sha256": seed.hex()})
    digests = {
        "primary_targets_sha256": previous.freeze(HERE / "targets.jsonl", points[:1]),
        "primary_answers_sha256": previous.freeze(HERE / "answers.jsonl", answers[:1]),
        "reserve_targets_sha256": previous.freeze(HERE / "reserve_targets.jsonl", points[1:]),
        "reserve_answers_sha256": previous.freeze(HERE / "reserve_answers.jsonl", answers[1:]),
    }
    print(json.dumps({"primary_count": 1, "reserve_count": 4, **digests}, sort_keys=True))


if __name__ == "__main__":
    main()
