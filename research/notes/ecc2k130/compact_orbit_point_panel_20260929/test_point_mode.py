#!/usr/bin/env python3
"""Exercise generated and public-point rho modes on an n13 known answer."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import tempfile
import sys


def call(rho: Path, corpus: str, points: Path | None = None,
         generate: bool = False) -> list[dict]:
    env = dict(os.environ, KIC_RHO_BATCH_CORPUS=corpus, KIC_RHO_DP_BITS="4")
    if points is not None:
        env["KIC_RHO_POINT_INPUT"] = str(points)
    if generate:
        env["KIC_RHO_GENERATE_ONLY"] = "1"
    output = subprocess.check_output([str(rho), "13", "0", "signed_frobenius",
                                      "1", "531310"], env=env, text=True)
    return [json.loads(line) for line in output.splitlines()]


def main() -> None:
    rho = Path(sys.argv[1]).resolve()
    corpus = "point-panel-n13-control-v1"
    fixture, = call(rho, corpus, generate=True)
    assert fixture["kind"] == "rho_ks_public_fixture"
    assert fixture["published_fixture_scalar"] > 0
    legacy_fixture, legacy_summary = call(rho, corpus)
    assert legacy_fixture["published_q"] == fixture["published_q"]
    assert legacy_fixture["published_fixture_scalar"] == fixture["published_fixture_scalar"]
    assert legacy_summary["target_source"] == "generated_fixture"
    with tempfile.TemporaryDirectory() as temp:
        points = Path(temp) / "points.jsonl"
        points.write_text(json.dumps(fixture["published_q"]) + "\n")
        point_fixture, point_summary = call(rho, corpus, points=points)
        assert point_fixture["published_q"] == fixture["published_q"]
        assert point_fixture["published_fixture_scalar"] is None
        assert point_fixture["recovered_fixture_scalar"] == fixture["published_fixture_scalar"]
        assert point_fixture["target_source"] == "public_point_jsonl"
        assert point_summary["target_source"] == "public_point_jsonl"
        assert point_summary["all_verified"] is True
        assert (legacy_summary["charges"]["scalar_multiplications"]
                == point_summary["charges"]["scalar_multiplications"] + 1)
        points.write_text(json.dumps([1 << 13, 0]) + "\n")
        invalid = subprocess.run(
            [str(rho), "13", "0", "signed_frobenius", "1", "531310"],
            env=dict(os.environ, KIC_RHO_BATCH_CORPUS=corpus,
                     KIC_RHO_POINT_INPUT=str(points)), capture_output=True)
        assert invalid.returncode != 0
    print(json.dumps({"status": "PASS", "n": 13,
                      "published_q": fixture["published_q"],
                      "scalar": fixture["published_fixture_scalar"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
