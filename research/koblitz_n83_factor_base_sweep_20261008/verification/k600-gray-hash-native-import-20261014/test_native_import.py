"""No-launch checks for the pinned native import commands."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import native_import


class NativeImportCommandTest(unittest.TestCase):
    def test_frozen_pairs_keep_hashes_and_archive_paths_in_worker_slots(self):
        plan = json.loads(native_import.preflight.PLAN.read_text())
        with TemporaryDirectory() as directory:
            archives = Path(directory)
            binary = archives / "worker"
            manifest = str(native_import.preflight.PILOT / "manifest.json")
            for arm in plan["arms"]:
                with self.subTest(role=arm["role"]):
                    self.assertEqual(native_import.verify_argv(binary, arm, archives), [
                        str(binary), "verify-cross-pilot",
                        arm["base_object_blake3"], manifest,
                        str(archives / f"{arm['base_object_blake3']}.jsonl.gz"),
                        arm["pool_object_blake3"], manifest,
                        str(archives / f"{arm['pool_object_blake3']}.jsonl.gz"),
                    ])


if __name__ == "__main__":
    unittest.main()
