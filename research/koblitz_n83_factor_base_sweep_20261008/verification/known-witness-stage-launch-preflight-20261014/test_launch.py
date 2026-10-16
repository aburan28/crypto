"""Check stage-launch custody, watchdog cleanup and authorization rejection."""

from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory
import unittest

import launch


EVIDENCE = launch.HERE / "selftest-evidence"


class StageLaunchTests(unittest.TestCase):
    def test_retained_marker_only_evidence_replays(self):
        row = launch.audit(EVIDENCE)
        self.assertEqual((row["status"], row["mode"],
                          row["trace"]["marker_count"]),
                         ("PASS_marker_transport_only", "self_test", 2))
        self.assertEqual(row["conservative_new_charge_seconds"], 0)
        self.assertEqual(row["cumulative_charged_seconds"], 22200)
        self.assertIsNone(row["selected_best_total_runtime"])

    def test_raw_marker_mutation_fails_audit(self):
        with TemporaryDirectory() as directory:
            folder = Path(directory)
            for name in ("launch.json", "exit.json", "stdout.raw", "stderr.raw"):
                shutil.copyfile(EVIDENCE / name, folder / name)
            with (folder / "stderr.raw").open("ab") as output:
                output.write(b"altered\n")
            with self.assertRaisesRegex(ValueError, "raw output differs"):
                launch.audit(folder)

    def test_watchdog_preserves_raw_exit_and_cleans_process_group(self):
        with TemporaryDirectory() as directory:
            folder = Path(directory)
            row = launch.launch_child(
                [sys.executable, "-c", "import time; time.sleep(30)"],
                0.05, folder)
            self.assertTrue(row["timed_out"])
            self.assertEqual(row["watchdog_exit_code"], 124)
            self.assertTrue(row["process_group_cleanup_verified"])
            self.assertIsNone(row["cleanup_error"])
            self.assertEqual(row["stdout_sha256"],
                             launch.sha256(folder / "stdout.raw"))
            self.assertEqual(row["stderr_sha256"],
                             launch.sha256(folder / "stderr.raw"))

    def test_full_domain_refuses_missing_or_invalid_grant_before_launch(self):
        with TemporaryDirectory() as directory:
            parent = Path(directory)
            output = parent / "unlaunched"
            with self.assertRaisesRegex(ValueError, "needs authorization"):
                launch.run("full_domain", output)
            self.assertFalse(output.exists())
            invalid = parent / "invalid.json"
            invalid.write_text("{}\n")
            with self.assertRaisesRegex(ValueError, "authorization differs"):
                launch.run("full_domain", output, invalid)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
