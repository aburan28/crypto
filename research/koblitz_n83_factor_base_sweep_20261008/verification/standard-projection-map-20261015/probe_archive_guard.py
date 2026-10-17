#!/usr/bin/env python3
"""Exercise the frozen-byte upload admission without calling S3."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ARCHIVER = HERE / "archive_maps.py"


def load_archiver(name: str):
    spec = importlib.util.spec_from_file_location(name, ARCHIVER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: probe_archive_guard.py SCRATCH BINARY BASE_ROOT")
    arguments = [str(ARCHIVER), *sys.argv[1:], "--upload"]
    prepared = (HERE / "local-prepared.json").read_bytes()

    positive = load_archiver("n83_archive_guard_positive")
    positive.put = lambda *args, **kwargs: (_ for _ in ()).throw(
        RuntimeError("NETWORK_BOUNDARY_AFTER_FROZEN_MATCH"))
    sys.argv = arguments
    try:
        positive.main()
    except RuntimeError as error:
        if str(error) != "NETWORK_BOUNDARY_AFTER_FROZEN_MATCH":
            raise
    else:
        raise AssertionError("positive probe did not reach the network boundary")

    negative = load_archiver("n83_archive_guard_negative")
    negative.put = lambda *args, **kwargs: (_ for _ in ()).throw(
        AssertionError("forged plan reached network"))
    original_read_text = Path.read_text

    def forged_read_text(path: Path, *args, **kwargs):
        content = original_read_text(path, *args, **kwargs)
        if path == HERE / "local-prepared.json":
            plan = json.loads(content)
            plan["rows"][0]["object_sha256"] = "0" * 64
            return json.dumps(plan)
        return content

    Path.read_text = forged_read_text
    try:
        negative.main()
    except ValueError as error:
        if str(error) != "dimension 12 differs from prepared upload":
            raise
    else:
        raise AssertionError("forged plan was accepted")
    finally:
        Path.read_text = original_read_text

    receipt = {
        "schema": "n83.standard-projection-map-upload-guard-probe/v1",
        "status": "PASS_FROZEN_MATCH_AND_TAMPER_REJECTION_BEFORE_NETWORK",
        "archiver_sha256": hashlib.sha256(ARCHIVER.read_bytes()).hexdigest(),
        "prepared_sha256": hashlib.sha256(prepared).hexdigest(),
        "positive": "REACHED_NETWORK_BOUNDARY_WITH_MATCHING_PREPARED_BYTES",
        "negative": "FORGED_D12_OBJECT_HASH_REJECTED_BEFORE_NETWORK",
        "network_calls": 0,
    }
    with (HERE / "upload-guard-probe.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
