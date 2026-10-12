#!/usr/bin/env python3
"""Hash every committed Stage 171 artifact into one deterministic seal."""

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
    for path in sorted(
        p
        for p in HERE.rglob("*")
        if p.is_file()
        and p != SEAL
        and "__pycache__" not in p.parts
        and p.suffix != ".pyc"
        and "development/clean-t14-b64/clean-t14-b64/" not in p.relative_to(HERE).as_posix()
    ):
        data = path.read_bytes()
        inventory.append(
            {
                "path": path.relative_to(HERE).as_posix(),
                "bytes": len(data),
                "sha256": digest(data),
            }
        )
    inventory_sha256 = digest(
        json.dumps(inventory, sort_keys=True, separators=(",", ":")).encode()
    )
    result_sha256 = digest((HERE / "result.json").read_bytes())
    payload = {
        "schema": "koblitz_stage171_dense_symbolic_sets_seal.v1",
        "status": "stage171_dense_symbolic_sets_result_frozen",
        "result_sha256": result_sha256,
        "inventory_sha256": inventory_sha256,
    }
    seal = {
        **payload,
        "inventory": inventory,
        "seal_payload_sha256": digest(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        ),
    }
    SEAL.write_text(json.dumps(seal, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in seal.items() if k != "inventory"}, indent=2))


if __name__ == "__main__":
    main()
