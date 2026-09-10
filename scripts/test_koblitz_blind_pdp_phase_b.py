#!/usr/bin/env python3
"""Adversarial and tiny end-to-end checks for blinded PDP Phase B."""

from __future__ import annotations

import argparse
from collections import Counter
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

import run_koblitz_blind_pdp_phase_b as phase_b
import build_koblitz_phase_b_wdsat as wdsat_builder
import score_koblitz_blind_pdp_phase_b as phase_b_score


PHASE_A_BASE = "47235e51b74a6fa8f3c8dc85d68bf886e20a1e88"
TEST_WDSAT_MAXIMUM = 100000000
TEST_WDSAT_MACROS = (
    "__MAX_ANF_ID__",
    "__MAX_DEGREE__",
    "__MAX_ID__",
    "__MAX_BUFFER_SIZE__",
    "__MAX_EQ__",
    "__MAX_EQ_SIZE__",
    "__MAX_XEQ__",
    "__MAX_XEQ_SIZE__",
)


def test_wdsat_config(maximum: int = TEST_WDSAT_MAXIMUM) -> str:
    return "\n".join(f"#define {name} {maximum}" for name in TEST_WDSAT_MACROS) + "\n"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str], *, env: dict[str, str] | None = None, expected: int = 0) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(command, text=True, capture_output=True, env=env, check=False)
    if completed.returncode != expected:
        raise AssertionError(
            f"command returned {completed.returncode}, expected {expected}: {command}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return completed


def tiny_phase_a(preparer: Path, root: Path) -> tuple[Path, Path, Path]:
    protocol = root / "tiny-phase-a-protocol.json"
    protocol.write_text(
        json.dumps(
            {
                "schema": "koblitz_balanced_pdp_phase_a_protocol.v1",
                "status": "frozen_before_preparation",
                "master_seed_hex": "00112233445566778899aabbccddeeff00112233445566778899aabbccddeeff",
                "master_seed_anchor_hex": None,
                "quotas_per_cell": {"decomposable": 2, "nondecomposable": 2},
                "max_unconditional_draws_per_cell": 10000,
                "cells": [
                    {
                        "id": "tiny-n7-standard",
                        "n": 7,
                        "ell": 2,
                        "m": 3,
                        "basis": "standard",
                        "curve_a": 0,
                        "factor_index": 0,
                    }
                ],
            },
            indent=2,
        )
        + "\n"
    )
    prepared = root / "tiny-phase-a"
    run(
        [
            str(preparer),
            "prepare",
            "--protocol",
            str(protocol),
            "--output",
            str(prepared),
            "--max-canonical-triples",
            "1000",
        ]
    )
    seal_path = prepared / "seal.json"
    seal = json.loads(seal_path.read_text())
    # This is a deliberately synthetic integration fixture.  Local tests run
    # before their implementation commit exists, so bind the fixture to the
    # clean Phase-A base rather than mislabelling the dirty build checkout.
    seal["implementation"]["source_revision"] = {
        "available": True,
        "commit": PHASE_A_BASE,
        "dirty": False,
        "porcelain": [],
    }
    seal_path.write_text(json.dumps(seal, indent=2) + "\n")
    return protocol, prepared / "blind/bundle.json", seal_path


def tiny_phase_b_protocol(
    path: Path,
    phase_a_protocol: Path,
    bundle_path: Path,
    seal_path: Path,
) -> dict:
    bundle = json.loads(bundle_path.read_text())
    seal = json.loads(seal_path.read_text())
    instance = bundle["instances"][0]
    protocol = {
        "schema": phase_b.PROTOCOL_SCHEMA,
        "status": "frozen_before_solver_execution",
        "phase_a_binding": {
            "source_revision": PHASE_A_BASE,
            "protocol_sha256": sha256(phase_a_protocol),
            "seal_sha256": sha256(seal_path),
            "blind_bundle_sha256": sha256(bundle_path),
            "blind_bundle_schema": phase_b.BLIND_BUNDLE_SCHEMA,
            "blind_instance_schema": phase_b.BLIND_INSTANCE_SCHEMA,
            "instance_count": 4,
        },
        "solver_root": {
            "exact_files": ["binding.json", "blind/bundle.json"],
            "symlinks_allowed": False,
            "unexpected_files_allowed": False,
            "binding_schema": phase_b.BINDING_SCHEMA,
        },
        "cells": [
            {
                "id": instance["cell_id"],
                "n": instance["n"],
                "ell": instance["ell"],
                "m": instance["m"],
                "basis": instance["basis"],
                "curve_a": instance["curve_a"],
                "factor_index": instance["factor_index"],
                "expected_instances": 4,
            }
        ],
        "execution": {
            "backend_order": list(phase_b.BACKENDS),
            "expected_backend_runs": 12,
            "one_process_at_a_time": True,
            "per_process_watchdog_seconds": 120,
            "native_conflict_budget": 100000,
            "single_thread_requested": True,
            "unknown_or_timeout_is_unsat": False,
        },
        "wdsat_build": {
            "source_commit": "1" * 40,
            "frozen_config_path": "synthetic-test-config.h",
            "frozen_config_sha256": hashlib.sha256(test_wdsat_config().encode()).hexdigest(),
            "limits": {
                "max_anf_id": TEST_WDSAT_MAXIMUM,
                "max_degree": TEST_WDSAT_MAXIMUM,
                "max_id": TEST_WDSAT_MAXIMUM,
                "max_buffer_size": TEST_WDSAT_MAXIMUM,
                "max_eq": TEST_WDSAT_MAXIMUM,
                "max_eq_size": TEST_WDSAT_MAXIMUM,
                "max_xeq": TEST_WDSAT_MAXIMUM,
                "max_xeq_size": TEST_WDSAT_MAXIMUM,
            },
        },
        "tool_build_accounting": {
            "rust_release_build": {
                "status": "not_run",
                "required_for_full_cost_gate": True,
            },
            "cryptominisat_build": {
                "status": "not_run",
                "source_commit": "3" * 40,
                "required_for_full_cost_gate": True,
            },
            "full_cost_gate_passed": False,
        },
        "source_contract": {
            "forbidden_solver_material": sorted(phase_b.DEFAULT_FORBIDDEN_MATERIAL),
        },
    }
    path.write_text(json.dumps(protocol, indent=2) + "\n")
    # The Phase-A seal binds the Phase-A protocol bytes, not the Phase-B file.
    assert seal["protocol_sha256"] == protocol["phase_a_binding"]["protocol_sha256"]
    return protocol


def fake_unknown_solver(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import json
import os
import sys

forbidden_names = {"ORACLE_PATH", "TARGET_CLASS", "WITNESS_INDICES", "PHASE_B_SECRET"}
leaked = {key: value for key, value in os.environ.items() if key in forbidden_names}
open_fds = []
for descriptor in range(3, 128):
    try:
        os.fstat(descriptor)
    except OSError:
        continue
    open_fds.append(descriptor)
stdin_byte = sys.stdin.buffer.read(1)
print("PHASE_B_FAKE_AUDIT=" + json.dumps({
    "leaked": leaked,
    "open_fds": open_fds,
    "stdin_hex": stdin_byte.hex(),
    "argv": sys.argv[1:],
}, sort_keys=True))
print("UNKNOWN")
raise SystemExit(0)
"""
    )
    path.chmod(0o755)


def fake_wdsat_receipt(path: Path, binary: Path, *, maximum: int = TEST_WDSAT_MAXIMUM) -> None:
    zero_metrics = {
        "wall_seconds": 0.0,
        "user_seconds": 0.0,
        "system_seconds": 0.0,
        "total_core_seconds": 0.0,
        "single_core_seconds": 0.0,
        "peak_rss_bytes": 0,
        "meter": "fresh-process getrusage(RUSAGE_CHILDREN)",
    }
    config = test_wdsat_config(maximum)
    receipt = {
        "schema": "koblitz_pdp_phase_b_wdsat_build_receipt.v1",
        "status": "completed",
        "source_commit": "1" * 40,
        "source_dirty": False,
        "config": config,
        "config_sha256": hashlib.sha256(config.encode()).hexdigest(),
        "binary_sha256": sha256(binary),
        "limits": {
            "max_anf_id": maximum,
            "max_degree": maximum,
            "max_id": maximum,
            "max_buffer_size": maximum,
            "max_eq": maximum,
            "max_eq_size": maximum,
            "max_xeq": maximum,
            "max_xeq_size": maximum,
        },
        "build_processes": [
            {
                "role": role,
                "command": ["/synthetic/" + role],
                "returncode": 0,
                "timed_out": False,
                "orphan_group_terminated": False,
                "metrics": zero_metrics,
                "stdout_sha256": hashlib.sha256(b"").hexdigest(),
                "stderr_sha256": hashlib.sha256(b"").hexdigest(),
            }
            for role in ("source-copy", "config-install", "clean", "build", "binary-copy")
        ],
        "claim_boundary": "Synthetic test build receipt with separately visible zero-cost fixture metrics",
    }
    path.write_text(json.dumps(receipt, indent=2) + "\n")


class ContractUnitTests(unittest.TestCase):
    def test_frozen_production_protocol_is_internally_consistent(self) -> None:
        protocol, _ = phase_b.read_json(phase_b.DEFAULT_PROTOCOL, "production Phase-B protocol")
        phase_b.validate_protocol(protocol)
        self.assertEqual(protocol["phase_a_binding"]["instance_count"], 160)
        self.assertEqual(protocol["execution"]["expected_backend_runs"], 480)
        config_path = phase_b.REPO / protocol["wdsat_build"]["frozen_config_path"]
        self.assertEqual(sha256(config_path), protocol["wdsat_build"]["frozen_config_sha256"])
        self.assertEqual(
            wdsat_builder.parse_config(config_path.read_text()),
            protocol["wdsat_build"]["limits"],
        )

    def test_duplicate_json_key_rejected(self) -> None:
        with self.assertRaisesRegex(phase_b.PhaseBError, "duplicate JSON key"):
            phase_b.parse_json_bytes(b'{"a":1,"a":2}', "duplicate fixture")

    def test_ambient_secrets_are_not_in_child_environment(self) -> None:
        previous = dict(os.environ)
        try:
            os.environ["ORACLE_PATH"] = "/secret/oracle-ledger.json"
            os.environ["TARGET_CLASS"] = "decomposable"
            os.environ["PHASE_B_SECRET"] = "sentinel-value"
            child = phase_b.safe_child_environment()
        finally:
            os.environ.clear()
            os.environ.update(previous)
        self.assertNotIn("ORACLE_PATH", child)
        self.assertNotIn("TARGET_CLASS", child)
        self.assertNotIn("PHASE_B_SECRET", child)
        self.assertEqual(child["OMP_NUM_THREADS"], "1")
        self.assertEqual(child["RAYON_NUM_THREADS"], "1")

    def test_launch_boundary_rejects_paths_labels_and_inputs(self) -> None:
        markers = set(phase_b.DEFAULT_FORBIDDEN_MATERIAL)
        with self.assertRaises(phase_b.PhaseBError):
            phase_b.assert_launch_boundary(["solver", "/tmp/oracle-ledger.json"], {}, [], markers)
        with self.assertRaises(phase_b.PhaseBError):
            phase_b.assert_launch_boundary([], {"TARGET_CLASS": "decomposable"}, [], markers)
        with self.assertRaises(phase_b.PhaseBError):
            phase_b.assert_launch_boundary([], {}, [{"target_class": "decomposable"}], markers)

    def test_wdsat_capacity_fails_closed(self) -> None:
        receipt = {"limits": {name: 10 for name in (
            "max_anf_id", "max_degree", "max_id", "max_buffer_size",
            "max_eq", "max_eq_size", "max_xeq", "max_xeq_size"
        )}}
        requirements = [{"cell_id": "tiny", "requirements": dict(receipt["limits"], max_id=11)}]
        with self.assertRaisesRegex(phase_b.PhaseBError, "undersized"):
            phase_b.require_wdsat_capacity(receipt, requirements)

    def test_scorer_rejects_sat_without_exact_witness_validation(self) -> None:
        row = {
            "solver": "wdsat",
            "status": "sat",
            "source_model_valid": True,
            "source_witness_valid": False,
            "process": {
                "command": ["/solver"],
                "returncode": 0,
                "timed_out": False,
                "orphan_group_terminated": False,
                "metrics": {
                    "wall_seconds": 0.0,
                    "user_seconds": 0.0,
                    "system_seconds": 0.0,
                    "total_core_seconds": 0.0,
                    "single_core_seconds": 0.0,
                    "peak_rss_bytes": 0,
                    "meter": "fresh-process getrusage(RUSAGE_CHILDREN)",
                },
            },
        }
        with self.assertRaisesRegex(phase_b.PhaseBError, "lacks exact model and witness"):
            phase_b_score.validate_backend_admission(row, "b-" + "0" * 64)
        native_unsat = copy.deepcopy(row)
        native_unsat.update(
            {
                "solver": "native-xor",
                "status": "unsat",
                "source_model_valid": None,
                "source_witness_valid": None,
            }
        )
        native_unsat["process"]["returncode"] = 7
        with self.assertRaisesRegex(phase_b.PhaseBError, "native UNSAT"):
            phase_b_score.validate_backend_admission(
                native_unsat, "b-" + "1" * 64
            )
        with TemporaryDirectory() as directory:
            root = Path(directory)
            stdout = b"s UNSATISFIABLE\n"
            (root / "cms.stdout").write_bytes(stdout)
            cms_unsat = copy.deepcopy(row)
            cms_unsat.update(
                {
                    "solver": "cryptominisat",
                    "status": "unsat",
                    "source_model_valid": None,
                    "source_witness_valid": None,
                }
            )
            cms_unsat["process"].update(
                {
                    "command": ["/cms", "--threads", "1"],
                    "returncode": 0,
                    "stdout_path": "cms.stdout",
                    "stdout_bytes": len(stdout),
                    "stdout_sha256": hashlib.sha256(stdout).hexdigest(),
                }
            )
            with self.assertRaisesRegex(phase_b.PhaseBError, "CryptoMiniSat UNSAT"):
                phase_b_score.validate_backend_admission(
                    cms_unsat,
                    "b-" + "2" * 64,
                    {"wdsat_requirements": {"max_anf_id": 5}},
                    root,
                )

    def test_timeout_is_explicit_and_attempt_files_cannot_be_overwritten(self) -> None:
        with TemporaryDirectory(prefix="phase-b-timeout-") as directory:
            root = Path(directory)
            record = phase_b.run_metered(
                role="timeout-fixture",
                command=[sys.executable, "-c", "import time; time.sleep(10)"],
                cwd=root,
                task_root=root,
                input_paths=[],
                timeout=0.1,
                meter=phase_b.DEFAULT_METER,
                environment=phase_b.safe_child_environment(),
                markers=set(phase_b.DEFAULT_FORBIDDEN_MATERIAL),
            )
            self.assertTrue(record["timed_out"])
            self.assertNotEqual(record["returncode"], 0)
            with self.assertRaisesRegex(phase_b.PhaseBError, "already exists"):
                phase_b.run_metered(
                    role="timeout-fixture",
                    command=[sys.executable, "-c", "raise SystemExit(0)"],
                    cwd=root,
                    task_root=root,
                    input_paths=[],
                    timeout=0.1,
                    meter=phase_b.DEFAULT_METER,
                    environment=phase_b.safe_child_environment(),
                    markers=set(phase_b.DEFAULT_FORBIDDEN_MATERIAL),
                )

    def test_leftover_descendant_invalidates_attempt(self) -> None:
        with TemporaryDirectory(prefix="phase-b-orphan-") as directory:
            root = Path(directory)
            with self.assertRaisesRegex(phase_b.PhaseBError, "left descendant processes"):
                phase_b.run_metered(
                    role="orphan-fixture",
                    command=[
                        sys.executable,
                        "-c",
                        'import subprocess; subprocess.Popen(["/bin/sleep", "2"])',
                    ],
                    cwd=root,
                    task_root=root,
                    input_paths=[],
                    timeout=5,
                    meter=phase_b.DEFAULT_METER,
                    environment=phase_b.safe_child_environment(),
                    markers=set(phase_b.DEFAULT_FORBIDDEN_MATERIAL),
                )

    def test_misleading_or_conflicting_solver_terminals_are_errors(self) -> None:
        base_record = {
            "stdout_text": "NOT UNSAT\n",
            "stderr_text": "",
            "returncode": 0,
            "timed_out": False,
            "metrics": {},
            "command": [],
        }
        manifest = {"source_variables": 4, "exports": {"cryptominisat_xor_dimacs": {"variables": 4}}}
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "instance.anf").write_text("p cnf 4 0\n")
            wdsat = phase_b.external_result(
                solver="wdsat",
                record=base_record,
                manifest=manifest,
                manifest_path=root / "manifest.json",
                backend=root / "backend",
                task_root=root,
                instance_root=root,
                timeout=120,
                meter=root / "meter",
                environment={},
                markers=set(phase_b.DEFAULT_FORBIDDEN_MATERIAL),
                backend_sha256="0" * 64,
            )
            self.assertEqual(wdsat["status"], "solver_error")
            conflicting = dict(base_record, stdout_text="s SATISFIABLE\ns UNSATISFIABLE\n")
            cms = phase_b.external_result(
                solver="cryptominisat",
                record=conflicting,
                manifest=manifest,
                manifest_path=root / "manifest.json",
                backend=root / "backend",
                task_root=root,
                instance_root=root,
                timeout=120,
                meter=root / "meter",
                environment={},
                markers=set(phase_b.DEFAULT_FORBIDDEN_MATERIAL),
                backend_sha256="0" * 64,
            )
            self.assertEqual(cms["status"], "solver_error")
            for stdout, returncode in (("0000\nUNSAT\n", 0), ("0000\n", 7)):
                contradictory_model = dict(
                    base_record, stdout_text=stdout, returncode=returncode
                )
                result = phase_b.external_result(
                    solver="wdsat",
                    record=contradictory_model,
                    manifest=manifest,
                    manifest_path=root / "manifest.json",
                    backend=root / "backend",
                    task_root=root,
                    instance_root=root,
                    timeout=120,
                    meter=root / "meter",
                    environment={},
                    markers=set(phase_b.DEFAULT_FORBIDDEN_MATERIAL),
                    backend_sha256="0" * 64,
                )
                self.assertEqual(result["status"], "solver_error")

    def test_explicit_manifest_must_be_v2_and_exact(self) -> None:
        protocol = {"source_contract": {"forbidden_solver_material": sorted(phase_b.DEFAULT_FORBIDDEN_MATERIAL)}}
        instance = {
            "n": 7,
            "ell": 2,
            "m": 3,
            "basis": "standard",
            "factor_index": 0,
            "source_nonce": 3,
            "blind_instance_id": "b-" + "0" * 64,
            "curve_a": 0,
            "target": {"x": "1", "y": "2"},
        }
        manifest = {
            "kind": "binary_koblitz_pdp_cross_solver_instance",
            "n": 7,
            "ell": 2,
            "m": 3,
            "seed": 3,
            "blind_instance_id": instance["blind_instance_id"],
            "target_mode": "explicit_affine",
            "curve_a": 0,
            "target": instance["target"],
            "factor_base_predicate": {
                "kind": "polynomial_subspace",
                "enumerates_target_subgroup": False,
                "uses_discrete_log_labels": False,
            },
            "factor_base_basis_bitmasks": ["1", "2"],
            "native_sat": {"status": "not_run_in_export_process"},
            "direct_meet_in_the_middle": {"status": "not_run_in_export_process"},
            "source_instance": {"schema": "koblitz_pdp_source_instance.v1"},
        }
        with TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "manifest.json").write_text("{}\n")
            (root / "instance.anf").write_text("x\n")
            (root / "instance.xor.cnf").write_text("x\n")
            (root / "instance.magma").write_text("x\n")
            unsafe = copy.deepcopy(manifest)
            unsafe["factor_base_predicate"]["enumerates_target_subgroup"] = True
            with self.assertRaisesRegex(phase_b.PhaseBError, "public algebraic"):
                phase_b.validate_export_manifest(unsafe, instance, root, protocol)
            with self.assertRaisesRegex(phase_b.PhaseBError, "source identity v2"):
                phase_b.validate_export_manifest(manifest, instance, root, protocol)

    def test_wdsat_builder_emits_runnable_authenticated_receipt(self) -> None:
        with TemporaryDirectory(prefix="phase-b-builder-") as directory:
            root = Path(directory)
            source = root / "wdsat-source"
            (source / "src").mkdir(parents=True)
            (source / "src/config.h").write_text("#define PLACEHOLDER 1\n")
            (source / "src/main.c").write_text("int main(void) { return 0; }\n")
            (source / "src/makefile").write_text(
                "SHELL=/bin/sh\n"
                "CC?=cc\n"
                "all: ../wdsat_solver\n"
                "../wdsat_solver: main.c config.h\n"
                "\t$(CC) -O2 main.c -o ../wdsat_solver\n"
                "clean:\n"
                "\trm -f *.o ../wdsat_solver\n"
            )
            subprocess.run(["git", "init", "-q"], cwd=source, check=True)
            subprocess.run(["git", "add", "src"], cwd=source, check=True)
            subprocess.run(
                [
                    "git",
                    "-c",
                    "user.name=Phase B Test",
                    "-c",
                    "user.email=phase-b-test@example.invalid",
                    "commit",
                    "-q",
                    "-m",
                    "fixture",
                ],
                cwd=source,
                check=True,
            )
            commit = subprocess.run(
                ["git", "rev-parse", "HEAD"], cwd=source, text=True, capture_output=True, check=True
            ).stdout.strip()
            config = root / "config.h"
            config.write_text(test_wdsat_config())
            protocol, _ = phase_b.read_json(phase_b.DEFAULT_PROTOCOL, "production protocol fixture")
            protocol = copy.deepcopy(protocol)
            protocol["wdsat_build"] = {
                "source_commit": commit,
                "frozen_config_path": "config.h",
                "frozen_config_sha256": sha256(config),
                "limits": {
                    "max_anf_id": TEST_WDSAT_MAXIMUM,
                    "max_degree": TEST_WDSAT_MAXIMUM,
                    "max_id": TEST_WDSAT_MAXIMUM,
                    "max_buffer_size": TEST_WDSAT_MAXIMUM,
                    "max_eq": TEST_WDSAT_MAXIMUM,
                    "max_eq_size": TEST_WDSAT_MAXIMUM,
                    "max_xeq": TEST_WDSAT_MAXIMUM,
                    "max_xeq_size": TEST_WDSAT_MAXIMUM,
                },
            }
            protocol_path = root / "protocol.json"
            protocol_path.write_text(json.dumps(protocol, indent=2) + "\n")
            output = root / "build-output"
            seal = wdsat_builder.build(
                protocol_path, source, config, output, phase_b.DEFAULT_METER
            )
            self.assertEqual(seal["status"], "build_frozen")
            self.assertTrue((output / "wdsat_solver").is_file())
            receipt, receipt_bytes = phase_b.read_json(output / "receipt.json", "test WDSat receipt")
            normalized = phase_b.validate_wdsat_build_receipt(
                receipt,
                receipt_bytes,
                phase_b.executable_identity(output / "wdsat_solver", "test WDSat binary"),
                protocol,
            )
            self.assertEqual(
                [item["role"] for item in normalized["build_processes"]],
                ["source-copy", "config-install", "clean", "build", "binary-copy"],
            )


class TinyIntegrationTests(unittest.TestCase):
    preparer: Path
    exporter: Path
    backend: Path

    def test_blind_run_isolated_then_scored(self) -> None:
        with TemporaryDirectory(prefix="phase-b-test-") as directory:
            root = Path(directory)
            phase_a_protocol, bundle, phase_a_seal = tiny_phase_a(self.preparer, root)
            phase_b_protocol = root / "tiny-phase-b-protocol.json"
            tiny_phase_b_protocol(phase_b_protocol, phase_a_protocol, bundle, phase_a_seal)
            solver_root = root / "solver-input"
            stage = run(
                [
                    sys.executable,
                    str(Path(phase_b.__file__)),
                    "stage",
                    "--protocol",
                    str(phase_b_protocol),
                    "--phase-a-seal",
                    str(phase_a_seal),
                    "--blind-bundle",
                    str(bundle),
                    "--output",
                    str(solver_root),
                ]
            )
            self.assertEqual(json.loads(stage.stdout)["blind_instance_count"], 4)
            self.assertEqual(
                sorted(path.relative_to(solver_root).as_posix() for path in solver_root.rglob("*") if path.is_file()),
                ["binding.json", "blind/bundle.json"],
            )
            self.assertFalse((solver_root / "seal.json").exists())
            with self.assertRaises(phase_b.PhaseBError):
                phase_b.stage_solver_root(
                    phase_b_protocol, phase_a_seal, bundle, solver_root
                )

            extra_root = root / "solver-extra"
            shutil.copytree(solver_root, extra_root)
            (extra_root / "unexpected.json").write_text("{}\n")
            protocol_value, _ = phase_b.read_json(phase_b_protocol, "test protocol")
            with self.assertRaises(phase_b.PhaseBError):
                phase_b.load_solver_root(extra_root, protocol_value)

            linked_root = root / "solver-linked"
            shutil.copytree(solver_root, linked_root)
            (linked_root / "blind/bundle.json").unlink()
            (linked_root / "blind/bundle.json").symlink_to(bundle)
            with self.assertRaises(phase_b.PhaseBError):
                phase_b.load_solver_root(linked_root, protocol_value)

            wdsat = root / "fake-wdsat"
            cms = root / "fake-cms"
            fake_unknown_solver(wdsat)
            fake_unknown_solver(cms)
            build_receipt = root / "wdsat-build.json"
            fake_wdsat_receipt(build_receipt, wdsat)
            run_root = root / "solver-run"
            hostile_env = dict(os.environ)
            hostile_env.update(
                {
                    "ORACLE_PATH": "/secret/oracle-ledger.json",
                    "TARGET_CLASS": "decomposable",
                    "WITNESS_INDICES": "0,1,2",
                    "PHASE_B_SECRET": "sentinel-value",
                }
            )
            completed = run(
                [
                    sys.executable,
                    str(Path(phase_b.__file__)),
                    "run",
                    "--protocol",
                    str(phase_b_protocol),
                    "--solver-root",
                    str(solver_root),
                    "--output",
                    str(run_root),
                    "--exporter",
                    str(self.exporter),
                    "--backend",
                    str(self.backend),
                    "--wdsat",
                    str(wdsat),
                    "--wdsat-build-receipt",
                    str(build_receipt),
                    "--cryptominisat",
                    str(cms),
                    "--allow-dirty",
                ],
                env=hostile_env,
            )
            run_seal = json.loads(completed.stdout)
            self.assertTrue(run_seal["full_panel_complete"])
            self.assertEqual(run_seal["backend_outcomes"], 12)
            export_inventory = json.loads((run_root / "export-inventory.json").read_text())
            self.assertEqual(export_inventory["selected_instances"], 4)
            self.assertTrue(export_inventory["wdsat_capacity"]["capacity_verified"])
            manifests = sorted(run_root.glob("tasks/*/instance/manifest.json"))
            self.assertEqual(len(manifests), 4)
            for manifest_path in manifests:
                manifest = json.loads(manifest_path.read_text())
                self.assertEqual(manifest["target_mode"], "explicit_affine")
                self.assertEqual(manifest["source_instance"]["schema"], "koblitz_pdp_source_instance.v2")
                self.assertNotIn("planted_points", manifest)
            intents = sorted(run_root.glob("tasks/**/*.intent.json"))
            self.assertEqual(len([path for path in intents if path.name == "export.intent.json"]), 4)
            for path in intents:
                intent = json.loads(path.read_text())
                encoded = json.dumps(intent).lower()
                for forbidden in phase_b.DEFAULT_FORBIDDEN_MATERIAL:
                    self.assertNotIn(forbidden, encoded)
                self.assertEqual(intent["stdin"], "devnull")
                self.assertTrue(intent["close_fds"])
                self.assertEqual(intent["environment"]["OMP_NUM_THREADS"], "1")
                self.assertNotIn("ORACLE_PATH", intent["environment"])
                self.assertNotIn("TARGET_CLASS", intent["environment"])
                self.assertNotIn("PHASE_B_SECRET", intent["environment"])
            fake_outputs = list(run_root.glob("tasks/*/wdsat.stdout")) + list(
                run_root.glob("tasks/*/cryptominisat.stdout")
            )
            self.assertEqual(len(fake_outputs), 8)
            for path in fake_outputs:
                audit_line = next(
                    line for line in path.read_text().splitlines() if line.startswith("PHASE_B_FAKE_AUDIT=")
                )
                audit = json.loads(audit_line.split("=", 1)[1])
                self.assertEqual(audit["leaked"], {})
                self.assertEqual(audit["open_fds"], [])
                self.assertEqual(audit["stdin_hex"], "")
            task_results = [json.loads(path.read_text()) for path in run_root.glob("tasks/*/task-result.json")]
            self.assertEqual(len(task_results), 4)
            for task in task_results:
                self.assertEqual([row["solver"] for row in task["backends"]], list(phase_b.BACKENDS))
                self.assertEqual(task["backends"][1]["status"], "unknown_inconclusive")
                self.assertEqual(task["backends"][2]["status"], "unknown_inconclusive")

            score_root = root / "score"
            scorer = Path(phase_b.__file__).with_name("score_koblitz_blind_pdp_phase_b.py")
            score_completed = run(
                [
                    sys.executable,
                    str(scorer),
                    "--protocol",
                    str(phase_b_protocol),
                    "--solver-root",
                    str(solver_root),
                    "--run-root",
                    str(run_root),
                    "--phase-a-seal",
                    str(phase_a_seal),
                    "--oracle-ledger",
                    str(phase_a_seal.parent / "sealed-oracle/oracle-ledger.json"),
                    "--output",
                    str(score_root),
                    "--allow-smoke",
                ]
            )
            score_seal = json.loads(score_completed.stdout)
            score = json.loads((score_root / "score.json").read_text())
            self.assertEqual(score_seal["status"], "post_run_truth_scoring_complete")
            self.assertEqual(score["selected_instances"], 4)
            self.assertEqual(score["backend_rows"], 12)
            self.assertEqual(
                Counter(row["target_class"] for row in score["rows"]),
                Counter({"decomposable": 6, "nondecomposable": 6}),
            )

            fake_outer = root / "fake-outer.metrics.json"
            fake_outer.write_text(
                json.dumps(
                    {
                        "command": [
                            "/usr/bin/true",
                            "run",
                            "--output",
                            str(run_root),
                            "--solver-root",
                            str(solver_root),
                        ],
                        "returncode": 0,
                        "watchdog_seconds": 999,
                        "timed_out": False,
                        "orphan_group_terminated": False,
                        "metrics": {
                            "wall_seconds": 0.0,
                            "user_seconds": 0.0,
                            "system_seconds": 0.0,
                            "total_core_seconds": 0.0,
                            "single_core_seconds": 0.0,
                            "peak_rss_bytes": 0,
                            "meter": "fresh-process getrusage(RUSAGE_CHILDREN)",
                        },
                    },
                    indent=2,
                )
                + "\n"
            )
            fake_outer_output = root / "fake-outer-score"
            rejected_outer = run(
                [
                    sys.executable,
                    str(scorer),
                    "--protocol",
                    str(phase_b_protocol),
                    "--solver-root",
                    str(solver_root),
                    "--run-root",
                    str(run_root),
                    "--phase-a-seal",
                    str(phase_a_seal),
                    "--oracle-ledger",
                    str(phase_a_seal.parent / "sealed-oracle/oracle-ledger.json"),
                    "--output",
                    str(fake_outer_output),
                    "--outer-metrics",
                    str(fake_outer),
                    "--allow-smoke",
                ],
                expected=2,
            )
            self.assertIn("differs from the frozen execution command", rejected_outer.stderr)
            self.assertFalse(fake_outer_output.exists())

            run_seal_path = run_root / "run-seal.json"
            original_run_seal = run_seal_path.read_bytes()
            forged_run_seal = json.loads(original_run_seal)
            forged_run_seal["full_panel_complete"] = False
            run_seal_path.write_text(json.dumps(forged_run_seal, indent=2) + "\n")
            forged_output = root / "forged-score"
            forged = run(
                [
                    sys.executable,
                    str(scorer),
                    "--protocol",
                    str(phase_b_protocol),
                    "--solver-root",
                    str(solver_root),
                    "--run-root",
                    str(run_root),
                    "--phase-a-seal",
                    str(phase_a_seal),
                    "--oracle-ledger",
                    str(phase_a_seal.parent / "sealed-oracle/oracle-ledger.json"),
                    "--output",
                    str(forged_output),
                    "--allow-smoke",
                ],
                expected=2,
            )
            self.assertIn("self-hash", forged.stderr)
            self.assertFalse(forged_output.exists())
            run_seal_path.write_bytes(original_run_seal)

            tampered_task = sorted(run_root.glob("tasks/*/task-result.json"))[0]
            original = tampered_task.read_bytes()
            tampered_task.write_bytes(original.replace(b"terminal_outputs_recorded", b"terminal_outputs_tampered", 1))
            tamper_output = root / "tampered-score"
            failed = run(
                [
                    sys.executable,
                    str(scorer),
                    "--protocol",
                    str(phase_b_protocol),
                    "--solver-root",
                    str(solver_root),
                    "--run-root",
                    str(run_root),
                    "--phase-a-seal",
                    str(phase_a_seal),
                    "--oracle-ledger",
                    str(phase_a_seal.parent / "sealed-oracle/oracle-ledger.json"),
                    "--output",
                    str(tamper_output),
                    "--allow-smoke",
                ],
                expected=2,
            )
            self.assertIn("inventory differs", failed.stderr)
            self.assertFalse(tamper_output.exists())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparer", type=Path, required=True)
    parser.add_argument("--exporter", type=Path, required=True)
    parser.add_argument("--backend", type=Path, required=True)
    args, unittest_args = parser.parse_known_args()
    TinyIntegrationTests.preparer = args.preparer.resolve()
    TinyIntegrationTests.exporter = args.exporter.resolve()
    TinyIntegrationTests.backend = args.backend.resolve()
    unittest.main(argv=[sys.argv[0], *unittest_args])


if __name__ == "__main__":
    main()
