#!/usr/bin/env python3
"""Freeze 4,000 public N83 points before the root11 edge result is known."""

import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "research/n83_final_tournament_20261015/generate.py"

spec = importlib.util.spec_from_file_location("n83_final_fixture_generator", SOURCE)
final = importlib.util.module_from_spec(spec)
spec.loader.exec_module(final)
final.HERE = HERE
final.DOMAIN = b"n83-edge-confirm-fresh-4000-20261015-v1"
final.COUNT = 4000
final.PRIOR_SHA256.update({
    "n83_final_tournament_20261015/targets.jsonl":
        "cfa943ae538fcc41d2af8bb8fd1c8fc0f29adc340b55644e90382da25b213424",
    "n83_root11_edgek_20261015/targets.jsonl":
        "c75cbb92ef1a23057272cf3481cc649b12495135c75affd1c2610205a42b7aef",
    "n83_final_tier1_20261015/targets.jsonl":
        "34c530851eb135b7e348ac976caaae37c5883c52b52962ad496028c058b1dcda",
})


if __name__ == "__main__":
    final.main()
