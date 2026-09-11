#!/usr/bin/env python3
"""Exercise build receipt custody using tiny compiler fixtures, never PDP inputs."""

from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import build_koblitz_phase_b_tools as builder
import run_koblitz_blind_pdp_phase_b as phase_b
import score_koblitz_blind_pdp_phase_b as scorer


class BuildReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="phase-b-build-receipt-")
        cls.root = Path(cls.temporary.name)
        cls.source = cls.root / "source-repository"
        cls.source.mkdir()
        (cls.source / "src").mkdir()
        (cls.source / "examples").mkdir()
        (cls.source / "Cargo.toml").write_text('[package]\nname = "build-receipt-fixture"\nversion = "0.1.0"\nedition = "2021"\n')
        lock = cls.source / builder.RUST_LOCK_PATH
        lock.parent.mkdir(parents=True)
        lock.write_text('version = 3\n\n[[package]]\nname = "build-receipt-fixture"\nversion = "0.1.0"\n')
        (cls.source / "src/lib.rs").write_text('pub fn fixture() -> u32 { 7 }\n')
        for name in ("koblitz_pdp_export", "koblitz_pdp_backend"):
            (cls.source / f"examples/{name}.rs").write_text('fn main() { println!("compiler fixture only"); }\n')
        for arguments in (["init", "-q"], ["add", "Cargo.toml", builder.RUST_LOCK_PATH, "src", "examples"], ["-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-qm", "Tiny build fixture"]):
            subprocess.run(["git", *arguments], cwd=cls.source, check=True, capture_output=True)
        cls.output = cls.root / "measured-build"
        cls.receipt = builder.build("rust", cls.source, cls.output, jobs=2, timeout=120)
        cls.identities = cls.receipt["binaries"]
        cls.objects = builder.rust_source_objects(cls.source)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def setUp(self):
        self.capsule = Path(tempfile.mkdtemp(prefix="case-", dir=self.root)) / "capsule"
        builder.copy_capsule(self.output / "receipt.json", self.capsule)

    def validate(self, identities=None, objects=None):
        return builder.validate_receipt(self.capsule / "receipt.json", "rust", identities or self.identities, objects or self.objects)

    def change_receipt(self, change, *, rehash=True):
        value = json.loads((self.capsule / "receipt.json").read_text())
        change(value)
        if rehash:
            value.pop("receipt_payload_sha256")
            value["receipt_payload_sha256"] = phase_b.canonical_sha256(value)
        (self.capsule / "receipt.json").write_bytes(phase_b.pretty_bytes(value))

    def test_parallel_build_and_relocated_capsule(self):
        validated = self.validate()
        self.assertEqual(validated["receipt"]["requested_parallel_jobs"], 2)
        self.assertGreater(validated["receipt"]["resources"]["total_core_seconds"], 0)
        self.assertIsNone(validated["receipt"]["resources"]["single_core_elapsed_seconds"])
        self.assertIsNone(validated["receipt"]["resources"]["aggregate_process_tree_peak_rss_bytes"])
        self.assertFalse(validated["receipt"]["accounting_boundary"]["full_cost_gate_passed"])
        compilation = next(row for row in validated["receipt"]["build_processes"] if row["role"] == "build")
        self.assertIn("--jobs", compilation["command"])
        self.assertEqual(Path(self.identities["exporter"]["path"]).stat().st_nlink, 1)
        self.assertEqual(builder.git(self.source, "ls-files", "Cargo.lock"), "")
        self.assertEqual((self.output / "source/Cargo.lock").read_bytes(), (self.source / builder.RUST_LOCK_PATH).read_bytes())

    def test_build_command_substitution_with_rehashed_receipt(self):
        def change(value):
            next(row for row in value["build_processes"] if row["role"] == "build")["command"] = ["/usr/bin/true"]
        self.change_receipt(change)
        with self.assertRaisesRegex(phase_b.PhaseBError, "frozen build recipe"):
            self.validate()

    def test_incompatible_frozen_lock_stops_before_compilation(self):
        source = self.root / "incompatible-lock-source"
        shutil.copytree(self.source, source)
        manifest = source / "Cargo.toml"
        manifest.write_text(manifest.read_text().replace('version = "0.1.0"', 'version = "0.2.0"'))
        subprocess.run(["git", "add", "Cargo.toml"], cwd=source, check=True, capture_output=True)
        subprocess.run(["git", "-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-qm", "Incompatible lock fixture"], cwd=source, check=True, capture_output=True)
        output = self.root / "incompatible-lock-build"
        with self.assertRaisesRegex(phase_b.PhaseBError, "vendor failed"):
            builder.build("rust", source, output, jobs=1, timeout=120)
        self.assertTrue((output / "evidence/vendor.metrics.json").exists())
        self.assertFalse((output / "receipt.json").exists())
        self.assertFalse((output / "evidence/build.intent.json").exists())

    def test_only_fixed_unused_cms_submodules_are_excluded(self):
        def command(source, *args):
            if args == ("rev-parse", "HEAD"):
                return builder.CMS_COMMIT
            if args[0] == "status":
                return ""
            if args == ("ls-files", "--stage"):
                return "\n".join(f"160000 {commit} 0\t{path}" for path, commit in builder.CMS_UNUSED_SUBMODULES.items())
            self.fail(args)
        with patch.object(builder, "git", side_effect=command):
            self.assertEqual(builder.clean_commit(self.source, builder.CMS_COMMIT), builder.CMS_COMMIT)
            with self.assertRaisesRegex(phase_b.PhaseBError, "fixed unused-test"):
                builder.clean_commit(self.source)

    def test_scorer_corroborates_rust_inputs_with_recorded_revision(self):
        commit = builder.git(self.source, "rev-parse", "HEAD")
        plan = {
            "additional_tool_builds": {"rust": {}},
            "source_revision": {"commit": commit, "dirty": False, "porcelain": []},
            "rust_build_source_objects": self.objects,
            "evidence_class": "operational_smoke", "allow_dirty_requested": True,
        }
        protocol = {"phase_a_binding": {"source_revision": commit}}
        with patch.object(phase_b, "REPO", self.source):
            self.assertEqual(scorer.validate_build_provenance(plan, protocol, allow_smoke=True), {"rust": {}})
            bad = deepcopy(plan)
            bad["rust_build_source_objects"]["src"] = "0" * 40
            with self.assertRaisesRegex(phase_b.PhaseBError, "recorded implementation revision"):
                scorer.validate_build_provenance(bad, protocol, allow_smoke=True)
            for label in ("made_up_measurement", "internal_measurement_pending_outer_and_tool_build_receipts"):
                bad = deepcopy(plan)
                bad["evidence_class"] = label
                with self.assertRaises(phase_b.PhaseBError):
                    scorer.validate_build_provenance(bad, protocol, allow_smoke=True)
            with self.assertRaisesRegex(phase_b.PhaseBError, "production scoring requires"):
                scorer.validate_build_provenance(plan, protocol, allow_smoke=False)

    def test_receipt_self_hash(self):
        self.change_receipt(lambda value: value.update(requested_parallel_jobs=9), rehash=False)
        with self.assertRaisesRegex(phase_b.PhaseBError, "payload hash"):
            self.validate()

    def test_rehashed_build_driver_substitution(self):
        self.change_receipt(lambda value: value.update(build_driver_sha256="0" * 64))
        with self.assertRaisesRegex(
            phase_b.PhaseBError, "build driver differs from the trusted current implementation"
        ):
            self.validate()

    def test_rehashed_process_meter_substitution(self):
        self.change_receipt(lambda value: value.update(process_meter_sha256="0" * 64))
        with self.assertRaisesRegex(
            phase_b.PhaseBError, "process meter differs from the trusted current implementation"
        ):
            self.validate()

    def test_rehashed_phase_b_helper_substitution(self):
        self.change_receipt(lambda value: value.update(phase_b_helper_sha256="0" * 64))
        with self.assertRaisesRegex(
            phase_b.PhaseBError, "phase b helper differs from the trusted current implementation"
        ):
            self.validate()

    def test_binary_substitution(self):
        wrong = deepcopy(self.identities)
        wrong["exporter"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(phase_b.PhaseBError, "requested executable"):
            self.validate(identities=wrong)

    def test_changed_rust_source(self):
        wrong = dict(self.objects)
        wrong["src"] = "0" * 40
        with self.assertRaisesRegex(phase_b.PhaseBError, "Rust build inputs"):
            self.validate(objects=wrong)

    def test_raw_metrics_and_output_tampering(self):
        for name in ("build.metrics.json", "build.stdout"):
            path = self.capsule / "evidence" / name
            original = path.read_bytes()
            path.write_bytes(original + b"tamper")
            with self.assertRaisesRegex(phase_b.PhaseBError, "evidence inventory"):
                self.validate()
            path.write_bytes(original)

    def test_rehashed_summary_disagrees_with_raw_metrics(self):
        self.change_receipt(lambda value: value["build_processes"][-1]["metrics"].update(total_core_seconds=0))
        with self.assertRaisesRegex(phase_b.PhaseBError, "raw receipt or intent"):
            self.validate()

    def test_missing_configure_or_compile_stage(self):
        self.change_receipt(lambda value: value["build_processes"].pop())
        with self.assertRaisesRegex(phase_b.PhaseBError, "required separately metered stages"):
            self.validate()

    def test_refuses_unknown_parallel_count_and_full_cost_promotion(self):
        self.change_receipt(lambda value: value["accounting_boundary"].update(full_cost_gate_passed=True))
        with self.assertRaisesRegex(phase_b.PhaseBError, "accounting boundary"):
            self.validate()

    def test_unexpected_file_or_symlink(self):
        unexpected = self.capsule / "evidence/extra"
        unexpected.write_text("extra")
        with self.assertRaisesRegex(phase_b.PhaseBError, "evidence inventory"):
            self.validate()
        unexpected.unlink()
        unexpected.symlink_to(self.output / "receipt.json")
        with self.assertRaises(phase_b.PhaseBError):
            self.validate()

    def test_no_build_resume(self):
        with self.assertRaisesRegex(phase_b.PhaseBError, "must be new"):
            builder.build("rust", self.source, self.output, 1, 120)

    def test_dirty_source_rejected_before_output_creation(self):
        dirty = self.source / "uncommitted"
        dirty.write_text("fixture")
        output = self.root / "dirty-build"
        try:
            with self.assertRaisesRegex(phase_b.PhaseBError, "must be clean"):
                builder.build("rust", self.source, output, 1, 120)
            self.assertFalse(output.exists())
        finally:
            dirty.unlink()

    def test_cryptominisat_commit_is_pinned_before_build(self):
        output = self.root / "wrong-cms-source"
        with self.assertRaisesRegex(phase_b.PhaseBError, "must be clean"):
            builder.build("cryptominisat", self.source, output, 1, 120)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
