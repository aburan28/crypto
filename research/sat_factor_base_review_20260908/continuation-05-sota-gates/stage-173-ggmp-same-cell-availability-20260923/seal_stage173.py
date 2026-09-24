#!/usr/bin/env python3
"""Seal every retained Stage 173 artifact."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
SEAL = HERE / "result-seal.json"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    inventory = []
    for path in sorted(path for path in HERE.rglob("*") if path.is_file() and path != SEAL):
        data = path.read_bytes()
        inventory.append(
            {"path": path.relative_to(HERE).as_posix(), "bytes": len(data), "sha256": digest(data)}
        )
    inventory_sha256 = digest(json.dumps(inventory, sort_keys=True, separators=(",", ":")).encode())
    payload = {
        "schema": "koblitz_stage173_ggmp_same_cell_availability_seal.v1",
        "status": "stage173_ggmp_same_cell_availability_result_frozen",
        "result_sha256": digest((HERE / "result.json").read_bytes()),
        "inventory_sha256": inventory_sha256,
    }
    seal = {
        **payload,
        "inventory": inventory,
        "seal_payload_sha256": digest(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()),
    }
    SEAL.write_text(json.dumps(seal, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in seal.items() if key != "inventory"}, indent=2))


if __name__ == "__main__":
    main()
