"""Tamper controls for the hosted 4+2 artifact audit and S3 custody."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parent


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


custody = load("custody")
replay = load("replay")


class EvidenceControls(unittest.TestCase):
    def test_custody_is_deterministic_and_detects_changed_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as location:
            root = Path(location)
            artifact = root / "artifact"
            (artifact / "n83-m6-setup").mkdir(parents=True)
            (artifact / "n83-m6-setup/source-commit.txt").write_text("17f03c5\n")
            (artifact / "base.jsonl.gz").write_bytes(b"public base\x00")
            first, second = root / "first.tar.gz", root / "second.tar.gz"
            receipt = custody.create(artifact, first)
            self.assertEqual(custody.verify(artifact, first, receipt)["status"],
                             "PASS_BYTE_IDENTICAL_CUSTODY")
            self.assertEqual(custody.create(artifact, second)["archive_sha256"],
                             receipt["archive_sha256"])
            (artifact / "base.jsonl.gz").write_bytes(b"wrong base\x00")
            with self.assertRaisesRegex(ValueError, "custody receipt differs"):
                custody.verify(artifact, first, receipt)

    def test_mount_replay_rejects_wrong_pool_source(self) -> None:
        arm = replay.plan.build()["arms"][0]
        mounts = [
            ("/worker", "/tmp/n83-m6-setup/worker", True),
            ("/base.jsonl.gz", "/tmp/n83-archives/" +
             arm["base_object_blake3"] + ".jsonl.gz", True),
            ("/pool.jsonl.gz", "/tmp/n83-archives/" +
             arm["pool_object_blake3"] + ".jsonl.gz", True),
            ("/pilot-manifest.json", "/tmp/research/koblitz_n83_factor_base_sweep_20261008/"
             "pilot-01/manifest.json", True),
            ("/corpus.json", "/tmp/research/koblitz_n83_factor_base_sweep_20261008/"
             "pilot-01/probe-corpus.json", True),
            ("/out", "/tmp/n83-m6-pair/gray_pool", False),
        ]
        command = ["docker", "run"]
        inspected = {"Mounts": []}
        for destination, source, readonly in mounts:
            command += ["--mount", f"type=bind,src={source},dst={destination}" +
                        (",readonly" if readonly else "")]
            inspected["Mounts"].append({
                "Type": "bind", "Destination": destination,
                "Source": source, "RW": not readonly,
            })
        replay.audit_mounts(command, inspected, arm, "search")
        inspected["Mounts"][2]["Source"] = "/tmp/unpinned-pool.jsonl.gz"
        with self.assertRaisesRegex(ValueError, "inspected bind differs"):
            replay.audit_mounts(command, inspected, arm, "search")


if __name__ == "__main__":
    unittest.main()
