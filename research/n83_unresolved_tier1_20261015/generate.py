#!/usr/bin/env python3
"""Freeze one public N83 point excluded from every earlier local panel."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
ANALYZER = REPO / "research/n83_indexed5_20261013/analyze.py"
FINAL = REPO / "research/n83_final_tournament_20261015/generate.py"
DOMAIN = b"n83-unresolved-tier1-one-fresh-target-20261015-v1"
EXTRA = {
    "n83_final_tournament_20261015/targets.jsonl":
        "cfa943ae538fcc41d2af8bb8fd1c8fc0f29adc340b55644e90382da25b213424",
    "n83_edge_confirmation_20261015/targets.jsonl":
        "66c4683e252f1cd321b4b6c7c9d06744c82140222ddce160faa4e7a0e0984b7e",
    "n83_final_tier1_20261015/targets.jsonl":
        "34c530851eb135b7e348ac976caaae37c5883c52b52962ad496028c058b1dcda",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_once(path, value):
    blob = json.dumps(value, separators=(",", ":")).encode() + b"\n"
    if path.exists() and path.read_bytes() != blob:
        raise ValueError(f"frozen point file differs: {path}")
    if not path.exists():
        path.write_bytes(blob)
    return digest(path)


def main():
    fspec = importlib.util.spec_from_file_location("n83_diagnostic_prior", FINAL)
    final = importlib.util.module_from_spec(fspec)
    fspec.loader.exec_module(final)
    gspec = importlib.util.spec_from_file_location("n83_diagnostic_group", ANALYZER)
    group = importlib.util.module_from_spec(gspec)
    gspec.loader.exec_module(group)
    if not group.on_curve(group.G) or group.multiply(group.G, group.R) is not None:
        raise ValueError("frozen subgroup generator or order differs")
    seen = set()
    for relative, expected in (final.PRIOR_SHA256 | EXTRA).items():
        path = REPO / "research" / relative
        if digest(path) != expected:
            raise ValueError(f"prior public panel changed: {relative}")
        for line in path.read_text().splitlines():
            record = json.loads(line)
            point = record if isinstance(record, list) else record.get("target") or record.get("point")
            if isinstance(point, list) and len(point) == 2:
                seen.add(tuple(point))
    seed = hashlib.sha256(DOMAIN + (0).to_bytes(4, "little")).digest()
    scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
    point = group.multiply(group.G, scalar)
    if point is None or not group.on_curve(point):
        raise ValueError("fresh public point is invalid")
    encoded = [str(point[0]), str(point[1])]
    if tuple(encoded) in seen:
        raise ValueError("fresh public point repeats a prior panel")
    targets_sha = write_once(HERE / "targets.jsonl", encoded)
    answers_sha = write_once(HERE / "answers.jsonl", {
        "fixture_index": 0, "point": encoded, "scalar": str(scalar),
        "seed_sha256": seed.hex()})
    print(json.dumps({"count": 1, "targets_sha256": targets_sha,
                      "answers_sha256": answers_sha}, sort_keys=True))


if __name__ == "__main__":
    main()
