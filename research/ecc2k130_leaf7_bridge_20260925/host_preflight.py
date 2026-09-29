#!/usr/bin/env python3
"""Check the frozen Sage host and hard resource cap without building a bridge."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
import resource

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    spec = json.loads((HERE / "FROZEN.json").read_text())
    assert spec["schema"] == "ecc2k130-leaf7-bridge-frozen-v1"
    assert spec["release_gate"] == "hold_host_preflight_and_review"
    assert spec["release_main_head"] is not None
    assert spec["host_image"]["manifest_sha256"] == (
        "2401ffa8e9fc85c7ea17d3649bde5958b4dbf0858b3e504098c4102720151711")
    assert platform.system() == "Linux" and platform.machine() == "x86_64"
    for relative, expected in (spec["input_sha256"] | spec["host_refusal_sha256"]).items():
        assert sha(REPO / relative) == expected, relative
    for name, expected in spec["implementation_sha256"].items():
        assert sha(HERE / name) == expected, name
    limit = int(spec["caps"]["child_peak_rss_bytes"])
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    assert resource.getrlimit(resource.RLIMIT_AS) == (limit, limit)
    from sage.version import version  # noqa: PLC0415
    assert str(version).startswith(spec["sage_version"]), version
    from sage.all import GF, EllipticCurve  # noqa: PLC0415
    field = GF(2)
    curve = EllipticCurve(field, [field.one(), field.zero(), field.zero(),
                                  field.zero(), field.one()])
    assert hasattr(curve, "division_polynomial") and hasattr(curve, "isogeny")
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    assert peak < limit
    print(json.dumps({
        "decision": "HOST_PREFLIGHT_PASS",
        "structural_status": "UNMEASURED",
        "freeze_sha256": sha(HERE / "FROZEN.json"),
        "image_manifest_sha256": spec["host_image"]["manifest_sha256"],
        "sage_version": str(version),
        "platform": platform.platform(),
        "rlimit_as_bytes": limit,
        "peak_rss_bytes": peak,
        "release_main_head": spec["release_main_head"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
