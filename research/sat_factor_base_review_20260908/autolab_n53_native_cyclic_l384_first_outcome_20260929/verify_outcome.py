#!/usr/bin/env python3
"""Replay the one archived n53 result without changing its measured freeze."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CAMPAIGN = "research/sat_factor_base_review_20260908/autolab_n53_native_cyclic_l384_20260929"
SOURCE_HEAD = "3c9232e0dc0f71aa7b08c5c844bd1377c95f4257"
EVIDENCE_HEAD = "cc8eef5c8f00705e8df9cab3c84f9ae2aba616c0"
RUN_ID = 36532278977
ARTIFACT_ID = 11018400513
ARTIFACT_SHA256 = "1ec46f66868eff353c640cb763751ed7ed64a8845ea9e57d467740c82a833bcb"
FROZEN_SHA256 = "1d30ca19ad2859dfbc605a205f7faa2e243337a79929fa07957af247da346fce"
FILES = {
    "archive_manifest.json": "70d4da30e319df065e497fe17281cb4011cfa76d543b889ce595351568e696a9",
    "evidence.tar.gz": "4591db23e9e295cadcd1f27434a84d7e1cdd87be392d59019d4599c0dc2ec8a6",
    "panel.json": "aa02885af461ae7c0204f6f1b07787e4567db9e9772c1d8c8aad6dbe9b9af369",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git(*args: str) -> bytes:
    return subprocess.check_output(["git", *args], cwd=ROOT)


def require(ok: bool, message: str) -> None:
    if not ok:
        raise AssertionError(message)


def main() -> None:
    require(subprocess.run(["git", "merge-base", "--is-ancestor", SOURCE_HEAD,
                            EVIDENCE_HEAD], cwd=ROOT).returncode == 0,
            "measured source does not precede bot evidence")
    require(subprocess.run(["git", "merge-base", "--is-ancestor", EVIDENCE_HEAD,
                            "HEAD"], cwd=ROOT).returncode == 0,
            "bot evidence is not in this checkout's history")
    frozen_path = ROOT / CAMPAIGN / "FROZEN.json"
    require(sha(frozen_path.read_bytes()) == FROZEN_SHA256 and
            git("show", SOURCE_HEAD + ":" + CAMPAIGN + "/FROZEN.json") ==
            frozen_path.read_bytes(), "measured freeze changed")
    require(not (ROOT / CAMPAIGN / "evidence/archive_manifest.json").exists(),
            "pre-outcome namespace was repopulated")
    manifest = json.loads((HERE / "archive_manifest.json").read_text())
    require(manifest["classification"] == "COMPLETE_PRIMARY_B_NO_CROSSOVER" and
            manifest["files"] == 64 and
            manifest["archive_sha256"] == FILES["evidence.tar.gz"] and
            manifest["panel_sha256"] == FILES["panel.json"],
            "first archive manifest drift")
    for name, digest in FILES.items():
        current = (HERE / name).read_bytes()
        require(sha(current) == digest, f"{name} digest drift")
        original = git("show", EVIDENCE_HEAD + ":" + CAMPAIGN + "/evidence/" + name)
        require(current == original, f"{name} differs from bot's first archive")
    zip_path = HERE / f"actions_artifact_{ARTIFACT_ID}.zip"
    require(sha(zip_path.read_bytes()) == ARTIFACT_SHA256,
            "first Actions artifact ZIP digest drift")
    # The original Actions ZIP and bot's deterministic tar must contain the
    # same 64 raw files, byte for byte; the tar's SHA256SUMS is checked again
    # by the frozen archive verifier below.
    with zipfile.ZipFile(zip_path) as zipped, tarfile.open(HERE / "evidence.tar.gz", "r:gz") as tar:
        zip_names = [info.filename for info in zipped.infolist() if not info.is_dir()]
        tar_names = [member.name.removeprefix("panel/") for member in tar.getmembers()
                     if member.name.startswith("panel/")]
        require(len(zip_names) == len(set(zip_names)) == len(tar_names) == 64 and
                set(zip_names) == set(tar_names), "raw artifact member set drift")
        for name in zip_names:
            path = PurePosixPath(name)
            require(not path.is_absolute() and ".." not in path.parts,
                    "unsafe raw artifact member")
            member = tar.getmember("panel/" + name)
            require(member.isfile() and zipped.read(name) == tar.extractfile(member).read(),
                    f"raw artifact member differs from first tar: {name}")
    panel = json.loads((HERE / "panel.json").read_text())
    require(panel["github_run_id"] == str(RUN_ID) and
            panel["checkout_head"] == SOURCE_HEAD and
            panel["classification"] == "COMPLETE_PRIMARY_B_NO_CROSSOVER",
            "measured run/head/outcome identity drift")
    command = [sys.executable, str(ROOT / CAMPAIGN / "verify_archive.py"),
               "--bundle", str(HERE)]
    checked = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                             timeout=1240)
    require(checked.returncode == 0,
            "frozen archive replay failed: " + checked.stderr[-2000:])
    report = json.loads(checked.stdout.strip().splitlines()[-1])
    require(report == {
        "archive_sha256": FILES["evidence.tar.gz"],
        "audit_classification": "INDEPENDENT_FULL_REPLAY",
        "classification": "COMPLETE_PRIMARY_B_NO_CROSSOVER",
        "verdict": "PASS_COMPLETE_INDEPENDENT_REPLAY",
    }, "frozen semantic replay verdict drift")
    print(json.dumps({"decision": "PASS_ARCHIVED_FIRST_OUTCOME",
                      "run_id": RUN_ID, "artifact_id": ARTIFACT_ID,
                      "artifact_sha256": ARTIFACT_SHA256, **report}, sort_keys=True))


if __name__ == "__main__":
    main()
