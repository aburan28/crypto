#!/usr/bin/env python3
"""Focused parser/acceptance tests for the Koblitz matched-matrix runner."""

from __future__ import annotations

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import run_koblitz_pdp_matrix as matrix


def run_record(stdout: str, returncode: int = 0) -> dict:
    return {
        "stdout": stdout,
        "stderr": "",
        "returncode": returncode,
        "timed_out": False,
        "metrics": {},
        "command": [],
    }


class MagmaTerminalTests(unittest.TestCase):
    def test_parses_complete_direct_f4_sat_record(self) -> None:
        output = """\
KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1
KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse
KOBLITZ_MAGMA_STATUS=SAT
KOBLITZ_MAGMA_F4_DEGREES=[ 3, 4, 5 ]
KOBLITZ_MAGMA_BASIS_SIZE=42
KOBLITZ_MAGMA_CPU_SECONDS=1.25
KOBLITZ_MAGMA_WALL_SECONDS=1.5
"""
        parsed = matrix.parse_magma_terminal(output)
        self.assertIsNotNone(parsed)
        self.assertEqual(parsed["terminal_status"], "sat")
        self.assertEqual(parsed["f4_step_degrees"], [3, 4, 5])
        row = matrix.solver_status(run_record(output), "magma-f4", Path("."), {})
        self.assertEqual(row["status"], "sat_basis_certificate_unverified_model")

    def test_parses_unit_basis_as_unsat(self) -> None:
        output = """\
KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1
KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse
KOBLITZ_MAGMA_STATUS=UNSAT
KOBLITZ_MAGMA_F4_DEGREES=[]
KOBLITZ_MAGMA_BASIS_SIZE=1
KOBLITZ_MAGMA_CPU_SECONDS=0.25
KOBLITZ_MAGMA_WALL_SECONDS=0.5
"""
        parsed = matrix.parse_magma_terminal(output)
        self.assertIsNotNone(parsed)
        row = matrix.solver_status(run_record(output), "magma-f4", Path("."), {})
        self.assertEqual(row["status"], "unsat")

    def test_accepts_single_nonunit_basis_as_unverified_sat(self) -> None:
        output = """\
KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1
KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse
KOBLITZ_MAGMA_STATUS=SAT
KOBLITZ_MAGMA_F4_DEGREES=[ 2 ]
KOBLITZ_MAGMA_BASIS_SIZE=1
KOBLITZ_MAGMA_CPU_SECONDS=0.1
KOBLITZ_MAGMA_WALL_SECONDS=0.2
"""
        parsed = matrix.parse_magma_terminal(output)
        self.assertIsNotNone(parsed)
        row = matrix.solver_status(run_record(output), "magma-f4", Path("."), {})
        self.assertEqual(row["status"], "sat_basis_certificate_unverified_model")

    def test_rejects_missing_duplicate_and_inconsistent_markers(self) -> None:
        complete = """\
KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1
KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse
KOBLITZ_MAGMA_STATUS=SAT
KOBLITZ_MAGMA_F4_DEGREES=[ 3 ]
KOBLITZ_MAGMA_BASIS_SIZE=2
KOBLITZ_MAGMA_CPU_SECONDS=0.1
KOBLITZ_MAGMA_WALL_SECONDS=0.2
"""
        self.assertIsNone(matrix.parse_magma_terminal(complete.replace(
            "KOBLITZ_MAGMA_WALL_SECONDS=0.2\n", ""
        )))
        self.assertIsNone(matrix.parse_magma_terminal(
            complete + "KOBLITZ_MAGMA_STATUS=UNSAT\n"
        ))
        self.assertIsNone(matrix.parse_magma_terminal(
            complete.replace("KOBLITZ_MAGMA_STATUS=SAT", "KOBLITZ_MAGMA_STATUS=UNSAT")
        ))
        self.assertIsNone(matrix.parse_magma_terminal(
            complete.replace("KOBLITZ_MAGMA_CPU_SECONDS=0.1", "KOBLITZ_MAGMA_CPU_SECONDS=nan")
        ))
        row = matrix.solver_status(run_record("Magma exited normally\n"), "magma-f4", Path("."), {})
        self.assertEqual(row["status"], "terminal_certificate_missing")


class SourceArtifactCustodyTests(unittest.TestCase):
    def test_sha256_snapshots_detect_same_length_mutation(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            exports = {}
            for name, filename, data in [
                ("wdsat_anf", "instance.anf", b"abcd"),
                ("cryptominisat_xor_dimacs", "instance.xor.cnf", b"efgh"),
                ("magma_boolean_f4", "instance.magma", b"ijkl"),
            ]:
                (root / filename).write_bytes(data)
                exports[name] = {"path": filename, "bytes": len(data), "blake3": "bound"}
            manifest = {"exports": exports}
            before = matrix.source_artifact_snapshot(root, manifest)
            (root / "instance.anf").write_bytes(b"wxyz")
            after = matrix.source_artifact_snapshot(root, manifest)
            self.assertNotEqual(before, after)
            self.assertNotEqual(
                before["wdsat_anf"]["sha256"], after["wdsat_anf"]["sha256"]
            )

    def test_snapshot_rejects_path_escape(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = {
                "exports": {
                    "wdsat_anf": {"path": "../outside", "bytes": 1, "blake3": "x"},
                    "cryptominisat_xor_dimacs": {
                        "path": "instance.xor.cnf",
                        "bytes": 1,
                        "blake3": "x",
                    },
                    "magma_boolean_f4": {
                        "path": "instance.magma",
                        "bytes": 1,
                        "blake3": "x",
                    },
                }
            }
            with self.assertRaises(ValueError):
                matrix.source_artifact_snapshot(root, manifest)


class AssignmentValidationTests(unittest.TestCase):
    def test_accepts_exact_source_bound_assignment_and_rejects_mismatch(self) -> None:
        import json

        with TemporaryDirectory() as directory:
            assignment_path = Path(directory) / "model.json"
            assignment_path.write_text("[true,false]\n")
            report = {
                "schema": "koblitz_pdp_assignment_validation.v1",
                "status": "valid_point_witness",
                "source_instance_id": "source-id",
                "source_instance_verified": True,
                "regenerated_source_exact": True,
                "assignment_values": 2,
                "assignment_blake3": "a" * 64,
                "source_assignment": [True, False],
                "source_model_valid": True,
                "source_witness_valid": True,
            }
            manifest = {
                "source_variables": 2,
                "source_instance": {"id_blake3": "source-id"},
            }
            accepted = matrix.assignment_validation_status(
                run_record(json.dumps(report)), manifest, assignment_path
            )
            self.assertEqual(accepted["status"], "valid_point_witness")
            report["source_assignment"] = [False, True]
            rejected = matrix.assignment_validation_status(
                run_record(json.dumps(report)), manifest, assignment_path
            )
            self.assertEqual(rejected["status"], "validation_contract_error")


class IsolatedBackendTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = {"source_instance": {"id_blake3": "source-id"}}
        self.base = {
            "schema": "koblitz_pdp_isolated_backend.v1",
            "backend": "native-sat",
            "status": "sat",
            "source_instance_id": "source-id",
            "source_instance_verified": True,
            "regenerated_source_exact": True,
            "source_artifacts": {
                "wdsat_anf": {"valid": True},
                "cryptominisat_xor_dimacs": {"valid": True},
                "magma_boolean_f4": {"valid": True},
            },
            "source_model_valid": True,
            "source_witness_valid": True,
            "stats": {"conflicts": 7},
        }

    def test_accepts_source_bound_valid_model(self) -> None:
        import json

        row = matrix.isolated_backend_status(
            run_record(json.dumps(self.base)), "native-sat", self.manifest
        )
        self.assertEqual(row["status"], "sat")
        self.assertEqual(row["conflicts"], 7)

    def test_rejects_invalid_model_or_wrong_source(self) -> None:
        import json

        invalid_model = dict(self.base, source_model_valid=False)
        invalid_witness = dict(self.base, source_witness_valid=False)
        wrong_source = dict(self.base, source_instance_id="different")
        self.assertEqual(
            matrix.isolated_backend_status(
                run_record(json.dumps(invalid_witness)), "native-sat", self.manifest
            )["status"],
            "backend_contract_error",
        )
        self.assertEqual(
            matrix.isolated_backend_status(
                run_record(json.dumps(invalid_model)), "native-sat", self.manifest
            )["status"],
            "backend_contract_error",
        )
        self.assertEqual(
            matrix.isolated_backend_status(
                run_record(json.dumps(wrong_source)), "native-sat", self.manifest
            )["status"],
            "backend_contract_error",
        )

    def test_accepts_model_cap_only_as_inconclusive(self) -> None:
        import json

        capped = dict(
            self.base,
            status="model_cap_inconclusive",
            source_model_valid=True,
            source_witness_valid=False,
        )
        row = matrix.isolated_backend_status(
            run_record(json.dumps(capped)), "native-sat", self.manifest
        )
        self.assertEqual(row["status"], "model_cap_inconclusive")


if __name__ == "__main__":
    unittest.main()
