#!/usr/bin/env python3
"""Non-scoring custody preflight for the frozen PDP difficulty protocol."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tarfile

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / "disjoint_cold_v2_outcome_20261001/evidence_run_36803331080"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    config_bytes = (HERE / "CONFIG.json").read_bytes()
    config = json.loads(config_bytes)
    manifest_bytes = (ARCHIVE / "MANIFEST.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    assert sha(manifest_bytes) == config["manifest_sha256"]
    assert manifest["source_head"] == config["source_head"]
    receipt = {"schema": "ecc2k130-pdp-difficulty-input-check-v1",
               "status": "PASS", "config_sha256": sha(config_bytes),
               "manifest_sha256": sha(manifest_bytes), "cells": {}}
    for cell, spec in config["cells"].items():
        entry = manifest["cases"][cell]
        assert entry["status"] == "PASS" and entry["second_host_matches_hosted"]
        assert spec["raw_sha256"] == entry["raw_sha256"]
        path = ARCHIVE / entry["raw_path"]
        data = path.read_bytes()
        assert len(data) == entry["raw_bytes"] and sha(data) == spec["raw_sha256"]
        prefix = f"disjoint-cold-v2-{cell}/"
        seen = set()
        with tarfile.open(path, "r:gz") as archive:
            for member in archive:
                assert member.isfile(), member.name
                assert member.name.startswith(prefix), member.name
                relative = member.name[len(prefix):]
                assert relative not in seen and relative in entry["member_sha256"]
                seen.add(relative)
                assert sha(archive.extractfile(member).read()) == entry["member_sha256"][relative], relative
            header_name = f"{prefix}{cell}/b00_ic_a.base.jsonl"
            header = json.loads(archive.extractfile(archive.getmember(header_name)).read().splitlines()[0])
            assert header["n"] == spec["n"] and header["orbit_columns"] == spec["K"]
            assert int(header["subgroup_order"]) == spec["subgroup_order"]
            assert header["field_modulus_low_terms"] == spec["field_modulus_low_terms"]
        assert seen == set(entry["member_sha256"])
        receipt["cells"][cell] = {"raw_sha256": sha(data), "members_verified": len(seen),
                                  "source_head": manifest["source_head"]}
    receipt["checker_sha256"] = sha(Path(__file__).read_bytes())
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({cell: row["members_verified"] for cell, row in receipt["cells"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
