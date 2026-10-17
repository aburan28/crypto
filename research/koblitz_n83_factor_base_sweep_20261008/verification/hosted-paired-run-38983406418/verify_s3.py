#!/usr/bin/env python3
"""Check fresh S3 GETs against the three hosted regenerated N83 archives."""

from __future__ import annotations

import argparse
import filecmp
import hashlib
import json
from pathlib import Path

from audit_replay import ARCHIVES, digest, require


PREFIX = ("s3://crypto-autoresearcher/factor-bases/icv1/etc/"
          "koblitz_n83_factor_base_sweep_20261008/a0/objects/")


def verify(downloaded: Path, hosted: Path) -> dict:
    rows = []
    for object_hash, (_, _, expected_sha) in ARCHIVES.items():
        name = f"{object_hash}.jsonl.gz"
        s3_copy = downloaded / name
        hosted_copy = hosted / name
        require(s3_copy.stat().st_size == hosted_copy.stat().st_size
                and digest(s3_copy, hashlib.sha256) == expected_sha
                and digest(hosted_copy, hashlib.sha256) == expected_sha
                and filecmp.cmp(s3_copy, hosted_copy, shallow=False),
                f"S3 and hosted archive bytes differ: {name}")
        rows.append({"s3_uri": PREFIX + name, "sha256": expected_sha,
                     "bytes": s3_copy.stat().st_size,
                     "hosted_byte_identical": True})
    return {"schema": "n83.hosted-paired-s3-roundtrip/v1",
            "status": "PASS_S3_HOSTED_BYTE_IDENTITY", "objects": rows,
            "selected_best_total_runtime": None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--downloaded", type=Path, required=True)
    parser.add_argument("--hosted", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.downloaded, args.hosted)
    with args.receipt.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write("\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
