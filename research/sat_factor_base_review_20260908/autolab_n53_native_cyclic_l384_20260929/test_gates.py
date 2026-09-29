"""Meaningful held controls for native base admission and dispatch one-shot."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import check_protocol
import cold_base
import dispatch_gate
import run_panel

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent / "autolab_n53_rank_rotation_20260925/evidence/evidence.tar.gz"


class NativeBaseGate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with tarfile.open(ARCHIVE, "r:gz") as tar:
            stream = tar.extractfile("panel/native_cyclic/producer.stdout.jsonl")
            assert stream is not None
            cls.raw = json.loads(stream.read())
        # The later merged producer reports its explicit default selection mode.
        cls.raw["compact_orbit_base_header"]["selection_mode"] = "ascending_x_v1"

    def materialize(self, raw):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw_path = root / "raw.jsonl"
            raw_path.write_text(json.dumps(raw) + "\n")
            output = root / "header.jsonl"
            receipt = root / "receipt.json"
            result = cold_base.materialize(raw_path, output, receipt)
            return result, json.loads(output.read_text())

    def test_exact_native_arrays_are_admitted_from_fresh_raw(self):
        result, header = self.materialize(self.raw)
        self.assertEqual(result["classification"], "FRESH_NATIVE_BASE_MATCHES_PINNED_ORDER")
        self.assertEqual(header["base_hash"], cold_base.BASE_ARRAYS_SHA256)
        self.assertEqual(result["scanned_x"], 465)

    def test_changed_point_label_is_rejected(self):
        raw = json.loads(json.dumps(self.raw))
        raw["compact_orbit_base_header"]["factor_base_point_labels"][0][0] ^= 1
        with self.assertRaises(AssertionError):
            self.materialize(raw)

    def test_loaded_base_is_rejected(self):
        raw = json.loads(json.dumps(self.raw))
        raw["factor_base_input_hash"] = "pretend"
        with self.assertRaises(AssertionError):
            self.materialize(raw)


class OneShotGate(unittest.TestCase):
    def event(self, root):
        path = root / "event.json"
        path.write_text(json.dumps({
            "action": "labeled", "label": {"name": dispatch_gate.LABEL},
            "pull_request": {"number": 999, "head": {
                "sha": "a" * 40, "ref": "codex/n53-native-cyclic-l384-20260929",
                "repo": {"full_name": "aburan28/crypto"}}},
        }))
        return path

    def exercise(self, job_conclusion):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            event_path = self.event(root)
            frozen = {"expected_branch": "codex/n53-native-cyclic-l384-20260929",
                      "release_main_head": "b" * 40}
            def fake_api(path, token):
                if "/pulls/999" in path:
                    return {"state": "open", "head": {"sha": "a" * 40}}
                if "/workflows/" in path:
                    return {"total_count": 1, "workflow_runs": [
                        {"id": 77, "head_branch": frozen["expected_branch"]}]}
                if "/runs/77/jobs" in path:
                    return {"total_count": 1, "jobs": [
                        {"name": "outcome", "id": 88, "status": "completed",
                         "conclusion": job_conclusion}]}
                raise AssertionError(path)
            env = {"GITHUB_REPOSITORY": "aburan28/crypto",
                   "GITHUB_RUN_ID": "100", "GITHUB_RUN_ATTEMPT": "1",
                   "GITHUB_EVENT_PATH": str(event_path), "GITHUB_TOKEN": "unit-test-only",
                   "KIC_NATIVE_L384_EXPECTED_HEAD": "a" * 40}
            panel = root / "panel"
            with patch.dict("os.environ", env), patch.object(check_protocol, "preflight", return_value=frozen),                  patch.object(dispatch_gate, "api", side_effect=fake_api):
                if job_conclusion == "skipped":
                    dispatch_gate.run(panel)
                else:
                    with self.assertRaises(AssertionError):
                        dispatch_gate.run(panel)
            return json.loads((panel / "predispatch.json").read_text())["status"]

    def test_prior_dispatched_job_blocks_relabel(self):
        self.assertEqual(self.exercise("failure"), "REFUSED")

    def test_prior_skipped_job_does_not_consume_attempt(self):
        self.assertEqual(self.exercise("skipped"), "ADMITTED")


class FailedAttemptArchive(unittest.TestCase):
    def verify(self, admitted: bool):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            panel = root / "panel"
            panel.mkdir()
            frozen = json.loads((HERE / "FROZEN.json").read_text())
            gate = {"schema": "n53_native_cyclic_l384_dispatch_gate_v1",
                    "status": "ADMITTED" if admitted else "REFUSED",
                    "label": frozen["one_shot_pr_label"], "run_attempt": "1",
                    "run_id": "unit-test"}
            (panel / "predispatch.json").write_text(json.dumps(gate) + "\n")
            if admitted:
                source = frozen["source_sha256"]
                build = {"schema": "n53_native_cyclic_l384_build_v1",
                         "command": frozen["build_command"], "status": "FAILED",
                         "cargo_lock_sha256": source["cargo_lock"],
                         "cargo_toml_sha256": source["cargo_toml"],
                         "ic_source_sha256": source["ic"],
                         "rho_source_sha256": source["rho"]}
                (panel / "build_receipt.json").write_text(json.dumps(build) + "\n")
            bundle = root / "bundle"
            subprocess.run([sys.executable, str(HERE / "archive.py"), "--panel", str(panel),
                            "--out", str(bundle)], check=True, capture_output=True, text=True)
            verified = subprocess.run([sys.executable, str(HERE / "verify_archive.py"),
                                       "--bundle", str(bundle)], check=True,
                                      capture_output=True, text=True)
            self.assertEqual(json.loads(verified.stdout)["verdict"],
                             "PASS_PREDISPATCH_OR_BUILD_RAW_INTEGRITY_ONLY")

    def test_refused_gate_is_preserved(self):
        self.verify(False)

    def test_failed_build_is_preserved(self):
        self.verify(True)


class ResourceBudgetControls(unittest.TestCase):
    def test_audit_slot_survives_exhausted_operational_budget(self):
        spent = {"rho": {"wall_ms": 600_000}, "B_training": {"wall_ms": 2_400_000}}
        self.assertEqual(run_panel.stage_timeout("A_training", 360, spent), 0)
        self.assertEqual(run_panel.stage_timeout("independent_audit", 1200, spent), 1200)
        self.assertEqual(run_panel.MAX_OPERATIONAL_WALL + run_panel.MAX_AUDIT_WALL,
                         run_panel.MAX_WALL)

    def test_late_exception_attempts_one_reserved_prefix_audit(self):
        with tempfile.TemporaryDirectory() as temp:
            panel = Path(temp); rho = panel / "rho"; rho.mkdir()
            receipt = {"returncode": 0, "timed_out": False, "rss_gate": False,
                       "peak_rss_bytes": 1, "observed_group_peak_rss_bytes": 1}
            (rho / "resource_receipt.json").write_text(json.dumps(receipt) + "\n")
            (rho / "rho.stdout.jsonl").write_text("complete-prefix\n")
            replay = {"classification": "INDEPENDENT_PARTIAL_REPLAY"}
            (panel / "audit.json").write_text(json.dumps(replay) + "\n")
            summary = {"classification": "INVALID_RUNNER_EXCEPTION",
                       "stage_sequence": ["rho"], "steps": {"rho": receipt}}
            with patch.object(run_panel, "measure", return_value=receipt) as measured:
                result = run_panel.audit_after_exception(panel, summary)
                self.assertEqual(result["stage_sequence"], ["rho", "independent_audit"])
                self.assertEqual(result["replay"], replay)
                self.assertEqual(measured.call_args.args[4], run_panel.MAX_AUDIT_WALL)
                run_panel.audit_after_exception(panel, result)
                measured.assert_called_once()

    @unittest.skipUnless(sys.platform.startswith("linux"), "Linux /proc resource control")
    def test_linux_wall_cap_kills_child_group(self):
        with tempfile.TemporaryDirectory() as temp:
            receipt = run_panel.measure(
                [sys.executable, "-c", "import time; time.sleep(5)"],
                {"PATH": os.environ.get("PATH", "")}, Path(temp) / "wall",
                "wall", 0.15, {"test_source": Path(__file__)})
            self.assertTrue(receipt["timed_out"])
            self.assertFalse(receipt["rss_gate"])
            self.assertLess(receipt["wall_ms"], 3000)

    @unittest.skipUnless(sys.platform.startswith("linux"), "Linux /proc resource control")
    def test_linux_group_rss_cap_kills_grandchild(self):
        child = "import time; data=bytearray(b'X' * (96 * 1024 * 1024)); time.sleep(10)"
        parent = ("import subprocess,sys,time; "
                  f"subprocess.Popen([sys.executable, '-c', {child!r}]); time.sleep(10)")
        limit = 32 * 1024**2
        with tempfile.TemporaryDirectory() as temp, patch.object(run_panel, "MAX_RSS", limit):
            receipt = run_panel.measure(
                [sys.executable, "-c", parent], {"PATH": os.environ.get("PATH", "")},
                Path(temp) / "rss", "rss", 6.0, {"test_source": Path(__file__)})
            self.assertTrue(receipt["rss_gate"])
            self.assertFalse(receipt["timed_out"])
            self.assertGreaterEqual(receipt["observed_group_peak_rss_bytes"], limit)


class ExceptionArchiveStatus(unittest.TestCase):
    def test_unreplayed_runner_failure_is_raw_integrity_only(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); panel = root / "panel"; panel.mkdir()
            frozen = json.loads((HERE / "FROZEN.json").read_text())
            gate = {"schema": "n53_native_cyclic_l384_dispatch_gate_v1",
                    "status": "ADMITTED", "label": frozen["one_shot_pr_label"],
                    "run_attempt": "1", "run_id": "unit-test", "reviewed_head": "a" * 40}
            (panel / "predispatch.json").write_text(json.dumps(gate) + "\n")
            failure_path = panel / "failure.json"
            failure_path.write_text('{"error":"synthetic late exception"}\n')
            summary = {"classification": "INVALID_RUNNER_EXCEPTION",
                       "failure_sha256": run_panel.sha(failure_path)}
            (panel / "panel.json").write_text(json.dumps(summary) + "\n")
            # An unverified audit-shaped file cannot promote this to a semantic result.
            (panel / "audit.json").write_text('{"classification":"INDEPENDENT_FULL_REPLAY"}\n')
            bundle = root / "bundle"
            subprocess.run([sys.executable, str(HERE / "archive.py"), "--panel", str(panel),
                            "--out", str(bundle)], check=True, capture_output=True, text=True)
            verified = subprocess.run([sys.executable, str(HERE / "verify_archive.py"),
                                       "--bundle", str(bundle)], check=True,
                                      capture_output=True, text=True)
            self.assertEqual(json.loads(verified.stdout)["verdict"],
                             "PASS_RUNNER_FAILURE_RAW_INTEGRITY_ONLY")


if __name__ == "__main__":
    unittest.main()
