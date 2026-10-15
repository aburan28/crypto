#!/usr/bin/env python3
"""Independently check the retained 900-second mixed-orbit import receipt."""

import hashlib
import json
from pathlib import Path
import subprocess


EVIDENCE = Path(__file__).resolve().parent
ROOT = EVIDENCE.parents[3]
COMMIT = "d68cee586ae4b59b18ac5885ef064fe3fe709772"
IMAGE = "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31"
OBJECTS = {
    "small": ("eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380",
              "76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92"),
    "envelope": ("aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a",
                 "82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851"),
}


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, label: str) -> None:
    if not condition:
        raise ValueError(label)


def check() -> dict:
    auth = json.loads((EVIDENCE / "authorization.json").read_text())
    outer = json.loads((EVIDENCE / "outer.json").read_text())
    budget = json.loads((EVIDENCE / "budget.json").read_text())
    imported = json.loads((EVIDENCE / "verify.stdout").read_text())
    require(outer["schema"] == "n83.mixed-orbit-import-900-outer/v1", "outer schema")
    require(outer["status"] == "PASS_import_only" and outer["import_passed"] is True,
            "import status")
    require(outer["source_commit"] == COMMIT and outer["container_image_id"] == IMAGE,
            "source or image")
    require(outer["authorization_sha256"] == sha(EVIDENCE / "authorization.json"),
            "authorization bytes")
    require(auth["schema"] == "n83.additional-local-budget/v1"
            and auth["purpose"] == "mixed_orbit_build_import_900_gate"
            and auth["prior_pilot_charged_seconds"] == "11700"
            and auth["additional_authorized_seconds"] == "900"
            and "900 additional local seconds" in auth["user_approval_reference"],
            "authorization terms")
    require(budget["schema"] == "n83.mixed-orbit-import-900-budget/v1"
            and budget["authorization_sha256"] == outer["authorization_sha256"]
            and budget["conservative_charge_seconds"] == 900
            and budget["remaining_additional_seconds"] == "0"
            and budget["status"] == outer["status"]
            and budget["selected_best_total_runtime"] is None,
            "budget terms")
    require(0 <= outer["elapsed_supervisor_seconds_L0"] <= 900
            and outer["conservative_charge_seconds"] == 900
            and outer["total_index_calculus_runtime_ms"] is None
            and outer["selected_best_total_runtime"] is None
            and outer["frozen_inputs_unchanged"] is True,
            "outer accounting")
    require(len(outer["stages"]) == 2
            and [item["name"] for item in outer["stages"]] == ["build", "verify"],
            "stage set")
    for item, wall, cpus, memory in zip(outer["stages"], (450, 240), ("4", "1"),
                                       ("8192m", "4096m")):
        command = item["command"]
        require(item["exit_code"] == 0 and item["timed_out"] is False
                and item["cleanup_verified"] is True and item["launch_error"] is None
                and item["container_state"]["ExitCode"] == 0
                and item["container_state"]["OOMKilled"] is False
                and item["wall_elapsed_seconds_L0"] <= wall + 30,
                f"{item['name']} disposition")
        require(command[command.index("--network") + 1] == "none"
                and command[command.index("--pull") + 1] == "never"
                and command[command.index("--platform") + 1] == "linux/arm64"
                and command[command.index("--cpus") + 1] == cpus
                and command[command.index("--memory") + 1] == memory
                and command[command.index("--memory-swap") + 1] == memory
                and IMAGE in command,
                f"{item['name']} isolation")
        for suffix in ("stdout", "stderr"):
            path = EVIDENCE / f"{item['name']}.{suffix}"
            require(sha(path) == item[f"{suffix}_sha256"], f"{item['name']} {suffix} bytes")
    require(outer["stages"][0]["command"][-2:] ==
            ["--example", "koblitz_n83_mixed_orbit_worker"], "build target")
    worker_args = outer["stages"][1]["command"][outer["stages"][1]["command"].index("/worker") + 1:]
    require(worker_args == ["verify", "/small-manifest.json", "/small.jsonl.gz",
                            "/envelope-manifest.json", "/envelope.jsonl.gz"], "verify only")
    for relative, digest in outer["source_sha256"].items():
        blob = subprocess.run(["git", "show", f"{COMMIT}:{relative}"], cwd=ROOT,
                              check=True, capture_output=True).stdout
        require(hashlib.sha256(blob).hexdigest() == digest, f"committed source {relative}")
    for path, digest in outer["input_sha256"].items():
        key = "small" if OBJECTS["small"][1] in path else (
            "envelope" if OBJECTS["envelope"][1] in path else None)
        if key:
            require(digest == OBJECTS[key][0], f"{key} input digest")
        elif Path(path).name == "authorization.json":
            require(digest == outer["authorization_sha256"], "authorization input digest")
        else:
            raise ValueError(f"unexpected input {path}")
        if Path(path).exists():
            require(sha(Path(path)) == digest, f"current input bytes {path}")
    require(imported == {
        "schema": "n83.mixed-orbit-import/v1", "status": "PASS",
        "curve_slug": "icv1-f2m83-tm6151469093347-debefd74",
        "small_object_blake3": OBJECTS["small"][1],
        "envelope_object_blake3": OBJECTS["envelope"][1],
        "small_points": 196212, "small_signed_orbits": 1182,
        "residual_points": 143756, "residual_signed_orbits": 866,
        "solver_search_executed": False,
    }, "worker import receipt")
    return {
        "schema": "n83.mixed-orbit-import-900-audit/v1", "status": "PASS",
        "source_commit": COMMIT, "stage_count": 2,
        "small_object_sha256": OBJECTS["small"][0],
        "envelope_object_sha256": OBJECTS["envelope"][0],
        "small_signed_orbits": 1182, "residual_signed_orbits": 866,
        "conservative_charge_seconds": 900,
        "elapsed_supervisor_seconds_L0": outer["elapsed_supervisor_seconds_L0"],
        "solver_search_executed": False,
        "total_index_calculus_runtime_ms": None,
    }


if __name__ == "__main__":
    print(json.dumps(check(), indent=2, sort_keys=True))
