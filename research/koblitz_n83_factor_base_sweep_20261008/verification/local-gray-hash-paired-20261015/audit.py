#!/usr/bin/env python3
"""Independently replay the local N83 prelaunch rejection; launch nothing."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
PLAN = HERE.parent / "k600-gray-hash-paired-protocol-20261014" / "plan.json"
IMAGE = "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31"
OBJECTS = {
    "3888e2597fd3519807c9d722f49f83329c93ce5cb0929501cf7cb7cbc3cab1a0":
        "8463aaf30365d6f6f5d6b189aff4e8afcb6c8305ec27c216700776cbdccba29d",
    "966a0493c96674f0af6d184615ed05309d08fba5c7342957b4186cd50366dc47":
        "df02963d1706fa7a64fd7558d88db79b23ea34801ac58062c4c500113e8c57fe",
    "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e":
        "3e85c9e419fcddd3c9faed3c84470bf96fcad9178e48e1d73629369335841ecf",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def audit(archives: Path) -> dict:
    receipt = json.loads((HERE / "prelaunch.json").read_text())
    require(receipt.get("schema") == "n83.local-paired-prelaunch/v1", "schema")
    require(receipt.get("status") == "HOST_BUSY_NO_SOLVER_LAUNCH", "disposition")
    require(receipt.get("n83_solver_launched") is False
            and receipt.get("conservative_charge_seconds") == 0, "no-launch charge")
    require(receipt.get("plan_sha256") == sha256(PLAN), "plan hash")
    require(receipt.get("authorization_sha256") == sha256(HERE / "authorization.json"),
            "authorization hash")
    expected_archives = {f"{object_hash}.jsonl.gz": digest for object_hash, digest in OBJECTS.items()}
    require(receipt.get("archive_sha256") == expected_archives, "archive receipt hashes")
    require(all(sha256(archives / name) == digest for name, digest in expected_archives.items()),
            "local archive bytes")
    probes = receipt.get("probes")
    require(isinstance(probes, dict) and all(
        isinstance(probes.get(name), dict) and probes[name].get("exit_code") == 0
        for name in ("top", "load", "cores", "containers", "image", "cgroup")),
        "probe execution")
    require(probes["image"].get("stdout", "").strip() == f"{IMAGE} arm64 linux",
            "image identity")
    require(probes["cgroup"].get("stdout", "").splitlines() == ["4294967296", "0"],
            "cgroup limits")
    require(probes["containers"].get("stdout", "").splitlines() == ["conductor-db"],
            "container isolation")
    top = probes["top"].get("stdout", "")
    match = re.search(r"CPU usage: .*?([0-9]+(?:\.[0-9]+)?)% idle", top)
    require(match is not None, "CPU idle parse")
    load = re.findall(r"[0-9]+(?:\.[0-9]+)?", probes["load"].get("stdout", ""))
    require(bool(load), "load parse")
    idle = float(match.group(1))
    one_minute = float(load[0])
    cpus = int(probes["cores"].get("stdout", "0").strip())
    measured = receipt.get("measured")
    require(isinstance(measured, dict) and measured.get("idle_percent") == idle
            and measured.get("load_1m") == one_minute
            and measured.get("logical_cpus") == cpus
            and measured.get("other_containers") == []
            and measured.get("system_free_bytes", 0) >= 10 * (1 << 30),
            "measured fields")
    require(idle < 30.0 or one_minute > cpus, "busy-host rejection condition")
    return {"schema": "n83.local-paired-prelaunch-audit/v1",
            "status": "PASS_PRELAUNCH_REJECTION_ONLY",
            "prelaunch_sha256": sha256(HERE / "prelaunch.json"),
            "authorization_sha256": sha256(HERE / "authorization.json"),
            "plan_sha256": sha256(PLAN),
            "idle_percent": idle, "load_1m": one_minute,
            "logical_cpus": cpus, "n83_solver_launched": False,
            "selected_best_total_runtime": None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archives", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.archives)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
