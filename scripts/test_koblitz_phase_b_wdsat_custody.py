#!/usr/bin/env python3
"""Test WDSat build-capsule custody with a tiny, non-PDP fixture."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import build_koblitz_phase_b_wdsat as builder
import run_koblitz_blind_pdp_phase_b as phase_b


CLEAN_IMPLEMENTATION = {"commit": "1" * 40, "dirty": False, "porcelain": []}


def config_text(limit: int = 128) -> str:
    return "".join(
        f"#define {name} {limit}\n"
        for name in (
            "__MAX_ANF_ID__", "__MAX_DEGREE__", "__MAX_ID__",
            "__MAX_BUFFER_SIZE__", "__MAX_EQ__", "__MAX_EQ_SIZE__",
            "__MAX_XEQ__", "__MAX_XEQ_SIZE__",
        )
    )


class WDSatCapsuleCustodyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory(prefix="wdsat-custody-")
        cls.root = Path(cls.temporary.name)
        cls.source = cls.root / "source"
        source_dir = cls.source / "src"
        source_dir.mkdir(parents=True)
        (source_dir / "makefile").write_text(
            "all:\n\tcp solver.sh ../wdsat_solver\n\tchmod 755 ../wdsat_solver\n"
            "clean:\n\trm -f *.o ../wdsat_solver\n"
        )
        (cls.source / ".gitignore").write_text("*.o\n")
        solver = source_dir / "solver.sh"
        solver.write_text("#!/bin/sh\nprintf 'fixture only\\n'\n")
        solver.chmod(0o755)
        for command in (
            ["git", "init", "-q"],
            ["git", "add", "."],
            [
                "git", "-c", "user.name=Fixture", "-c",
                "user.email=fixture@example.invalid", "commit", "-qm", "fixture",
            ],
        ):
            subprocess.run(command, cwd=cls.source, check=True, capture_output=True)
        cls.source_commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=cls.source, check=True,
            text=True, capture_output=True,
        ).stdout.strip()
        (source_dir / "stale.o").write_bytes(b"ignored stale object")
        cls.config = cls.root / "config.h"
        cls.config.write_text(config_text())
        cls.protocol, _ = phase_b.read_json(phase_b.DEFAULT_PROTOCOL, "test protocol base")
        cls.protocol = deepcopy(cls.protocol)
        cls.protocol["wdsat_build"] = {
            "source_commit": cls.source_commit,
            "frozen_config_path": "config.h",
            "frozen_config_sha256": phase_b.sha256_file(cls.config, "test config"),
            "limits": {
                "max_anf_id": 128, "max_degree": 128, "max_id": 128,
                "max_buffer_size": 128, "max_eq": 128, "max_eq_size": 128,
                "max_xeq": 128, "max_xeq_size": 128,
            },
        }
        cls.protocol_path = cls.root / "protocol.json"
        phase_b.write_json_new(cls.protocol_path, cls.protocol)
        cls.output = cls.root / "original-capsule"
        with patch.object(builder, "current_implementation_state", return_value=deepcopy(CLEAN_IMPLEMENTATION)):
            builder.build(
                cls.protocol_path, cls.source, cls.config, cls.output,
                phase_b.DEFAULT_METER,
            )
        cls.identity = phase_b.executable_identity(
            cls.output / "wdsat_solver", "test WDSat binary"
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def setUp(self) -> None:
        self.case_root = Path(tempfile.mkdtemp(prefix="case-", dir=self.root))
        self.capsule = self.case_root / "capsule"
        shutil.copytree(self.output, self.capsule)

    def validate(self, *, identity=None, state=None, clean=True):
        return builder.validate_capsule(
            self.capsule / "build-seal.json", self.protocol,
            identity or self.identity, state or CLEAN_IMPLEMENTATION,
            require_clean_implementation=clean,
        )

    def reseal(self, *, receipt_changed: bool = False) -> None:
        receipt_path = self.capsule / "receipt.json"
        if receipt_changed:
            receipt = json.loads(receipt_path.read_text())
            receipt.pop("receipt_payload_sha256", None)
            receipt["receipt_payload_sha256"] = phase_b.canonical_sha256(receipt)
            receipt_path.write_bytes(phase_b.pretty_bytes(receipt))
        seal_path = self.capsule / "build-seal.json"
        seal = json.loads(seal_path.read_text())
        seal["receipt_sha256"] = phase_b.sha256_file(receipt_path, "test receipt")
        seal["inventory"] = phase_b.all_regular_inventory(self.capsule, {"build-seal.json"})
        seal["inventory_sha256"] = phase_b.canonical_sha256(seal["inventory"])
        seal.pop("seal_payload_sha256", None)
        seal["seal_payload_sha256"] = phase_b.canonical_sha256(seal)
        seal_path.write_bytes(phase_b.pretty_bytes(seal))

    def test_valid_capsule_and_relocation(self) -> None:
        normalized = self.validate()
        self.assertEqual(normalized["binary_sha256"], self.identity["sha256"])
        self.assertEqual(
            [row["role"] for row in normalized["build_processes"]],
            builder.EXPECTED_ROLES,
        )
        destination = self.case_root / "relocated"
        copied = builder.copy_capsule(
            self.capsule / "build-seal.json", destination, self.protocol,
            self.identity, CLEAN_IMPLEMENTATION,
        )
        self.assertEqual(normalized, copied)

    def test_receipt_and_seal_self_hashes(self) -> None:
        receipt_path = self.capsule / "receipt.json"
        receipt = json.loads(receipt_path.read_text())
        receipt["status"] = "forged"
        receipt_path.write_bytes(phase_b.pretty_bytes(receipt))
        self.reseal()
        with self.assertRaisesRegex(phase_b.PhaseBError, "receipt payload self-hash"):
            self.validate()

        shutil.rmtree(self.capsule)
        shutil.copytree(self.output, self.capsule)
        seal_path = self.capsule / "build-seal.json"
        seal = json.loads(seal_path.read_text())
        seal["status"] = "forged"
        seal_path.write_bytes(phase_b.pretty_bytes(seal))
        with self.assertRaisesRegex(phase_b.PhaseBError, "seal payload self-hash|schema or status"):
            self.validate()

    def test_rehashed_recipe_and_raw_metrics_tampering(self) -> None:
        receipt_path = self.capsule / "receipt.json"
        receipt = json.loads(receipt_path.read_text())
        receipt["build_processes"][0]["command"] = ["/usr/bin/true"]
        receipt_path.write_bytes(phase_b.pretty_bytes(receipt))
        self.reseal(receipt_changed=True)
        with self.assertRaisesRegex(phase_b.PhaseBError, "frozen WDSat build recipe"):
            self.validate()

        shutil.rmtree(self.capsule)
        shutil.copytree(self.output, self.capsule)
        metrics_path = self.capsule / "build.metrics.json"
        metrics = json.loads(metrics_path.read_text())
        metrics["metrics"]["total_core_seconds"] = 0
        metrics_path.write_bytes(phase_b.pretty_bytes(metrics))
        self.reseal()
        with self.assertRaisesRegex(phase_b.PhaseBError, "exact raw metrics"):
            self.validate()

    def test_inventory_binary_and_implementation_guards(self) -> None:
        (self.capsule / "unexpected").write_text("tamper")
        with self.assertRaisesRegex(phase_b.PhaseBError, "inventory"):
            self.validate()
        (self.capsule / "unexpected").unlink()

        wrong = dict(self.identity)
        wrong["sha256"] = "0" * 64
        with self.assertRaisesRegex(phase_b.PhaseBError, "requested executable"):
            self.validate(identity=wrong)
        wrong_state = {"commit": "2" * 40, "dirty": False, "porcelain": []}
        with self.assertRaisesRegex(phase_b.PhaseBError, "expected implementation state"):
            self.validate(state=wrong_state)

    def test_driver_meter_and_source_inventory_substitution(self) -> None:
        receipt_path = self.capsule / "receipt.json"
        receipt = json.loads(receipt_path.read_text())
        receipt["build_driver_identity"]["sha256"] = "0" * 64
        receipt_path.write_bytes(phase_b.pretty_bytes(receipt))
        self.reseal(receipt_changed=True)
        with self.assertRaisesRegex(phase_b.PhaseBError, "trusted current implementation"):
            self.validate()

        shutil.rmtree(self.capsule)
        shutil.copytree(self.output, self.capsule)
        inventory_path = self.capsule / "configured-source-inventory.json"
        inventory = json.loads(inventory_path.read_text())
        inventory["inventory"][0]["bytes"] += 1
        inventory["inventory_sha256"] = phase_b.canonical_sha256(inventory["inventory"])
        inventory.pop("record_payload_sha256")
        inventory["record_payload_sha256"] = phase_b.canonical_sha256(inventory)
        inventory_path.write_bytes(phase_b.pretty_bytes(inventory))
        self.reseal()
        with self.assertRaisesRegex(phase_b.PhaseBError, "source input changed"):
            self.validate()

    def test_copy_validates_source_before_destination_creation(self) -> None:
        (self.capsule / "unexpected").write_text("tamper")
        destination = self.case_root / "must-not-exist"
        with self.assertRaisesRegex(phase_b.PhaseBError, "inventory"):
            builder.copy_capsule(
                self.capsule / "build-seal.json", destination, self.protocol,
                self.identity, CLEAN_IMPLEMENTATION,
            )
        self.assertFalse(destination.exists())

    def test_prebuild_inputs_allow_generated_objects_but_reject_source_mutation(self) -> None:
        receipt = json.loads((self.capsule / "receipt.json").read_text())
        inventory_path = self.capsule / "configured-source-inventory.json"
        generated = self.capsule / "build/src/generated.o"
        generated.write_bytes(b"generated build output")
        builder.validate_source_inventory(inventory_path, receipt)
        makefile = self.capsule / "build/src/makefile"
        makefile.write_bytes(makefile.read_bytes() + b"\n# changed\n")
        with self.assertRaisesRegex(phase_b.PhaseBError, "source input changed"):
            builder.validate_source_inventory(inventory_path, receipt)


if __name__ == "__main__":
    unittest.main()
