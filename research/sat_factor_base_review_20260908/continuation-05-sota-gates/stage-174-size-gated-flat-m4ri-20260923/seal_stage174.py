#!/usr/bin/env python3
"""Write the immutable Stage 174 file inventory."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "result-seal.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory() -> list[dict]:
    files = []
    for path in sorted(HERE.rglob("*")):
        if (
            not path.is_file()
            or path == OUTPUT
            or "__pycache__" in path.parts
            or path.suffix == ".pyc"
        ):
            continue
        files.append(
            {
                "path": path.relative_to(HERE).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return files


def main() -> None:
    files = inventory()
    result = HERE / "result.json"
    verification = HERE / "verification.json"
    if sys.argv[1:] == ["--verify"]:
        seal = json.loads(OUTPUT.read_text())
        assert seal["schema"] == "koblitz_stage174_result_seal.v1"
        assert seal["inventory_entries"] == len(files)
        assert seal["inventory"] == files
        assert seal["inventory_sha256"] == hashlib.sha256(
            json.dumps(files, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        assert seal["result_sha256"] == sha256(result)
        assert seal["verification_sha256"] == sha256(verification)
        print(json.dumps({"entries": len(files), "ok": True, "seal_sha256": sha256(OUTPUT)}))
        return
    if sys.argv[1:]:
        raise SystemExit("usage: seal_stage174.py [--verify]")
    seal = {
        "schema": "koblitz_stage174_result_seal.v1",
        "source_commit": "e51efb219edd4511a08d545de21184928e5e771c",
        "result_sha256": sha256(result),
        "verification_sha256": sha256(verification),
        "inventory_entries": len(files),
        "inventory": files,
        "inventory_sha256": hashlib.sha256(
            json.dumps(files, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "claim_boundary": "same-target native-F4 engineering; not SOTA",
    }
    OUTPUT.write_text(json.dumps(seal, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"entries": len(files), "result_sha256": seal["result_sha256"]}))


if __name__ == "__main__":
    main()
