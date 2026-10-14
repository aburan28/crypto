#!/usr/bin/env python3
"""Deterministically freeze five new N83 public points before panel selection."""

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
ANALYZER = ROOT / "research/n83_indexed5_20261013/analyze.py"
ANALYZER_SHA256 = "e34d06b3275d86dcb3a659b79adf87dcdc078cebb6f510c4b5c42403d83a1804"
HOLDOUT = ROOT / "research/n83_evaluation_20261012/optimization/fixtures/holdout_5_targets.jsonl"
PRIMARY = ROOT / "research/n83_indexed5_20261013/primary_target.jsonl"
DOMAIN = b"n83-tier1-final-20261014-v1"

spec = importlib.util.spec_from_file_location("n83_fixture_group", ANALYZER)
group = importlib.util.module_from_spec(spec)
spec.loader.exec_module(group)


def canonical_line(value):
    return json.dumps(value, separators=(",", ":"), ensure_ascii=False) + "\n"


def write_once(path, content):
    if path.exists():
        if path.read_text() != content:
            raise ValueError(f"frozen fixture differs: {path}")
    else:
        path.write_text(content)


def main():
    if hashlib.sha256(ANALYZER.read_bytes()).hexdigest() != ANALYZER_SHA256:
        raise ValueError("independent group implementation changed")
    previous = {tuple(json.loads(line)) for line in HOLDOUT.read_text().splitlines()}
    previous.add(tuple(json.loads(PRIMARY.read_text())))
    targets = []
    answers = []
    for index in range(5):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, "little")).digest()
        scalar = 1 + int.from_bytes(seed, "little") % (group.R - 1)
        point = group.multiply(group.G, scalar)
        if (point is None or not group.on_curve(point)
                or group.multiply(point, group.R) is not None):
            raise ValueError(f"fixture {index}: group verification failed")
        encoded = [str(point[0]), str(point[1])]
        if tuple(encoded) in previous:
            raise ValueError(f"fixture {index}: point was used by selection")
        previous.add(tuple(encoded))
        targets.append(encoded)
        answers.append({"fixture_index": index, "scalar": str(scalar),
                        "point": encoded, "seed_sha256": seed.hex()})
    write_once(HERE / "targets.jsonl", "".join(map(canonical_line, targets)))
    write_once(HERE / "answers.jsonl", "".join(map(canonical_line, answers)))
    print(json.dumps({"target_count": len(targets),
                      "targets_sha256": hashlib.sha256((HERE / "targets.jsonl").read_bytes()).hexdigest(),
                      "answers_sha256": hashlib.sha256((HERE / "answers.jsonl").read_bytes()).hexdigest()},
                     sort_keys=True))


if __name__ == "__main__":
    main()
