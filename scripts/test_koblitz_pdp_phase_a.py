#!/usr/bin/env python3
"""Tiny end-to-end checks for balanced PDP Phase-A preparation and export."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory


SEED = "00112233445566778899aabbccddeeff00112233445566778899aabbccddeeff"
FORBIDDEN_BLIND_TEXT = (
    "decomposable",
    "nondecomposable",
    "target_class",
    "witness_indices",
    "selection_priority",
    "unconditional_draw_index",
    "planted",
)


def run(command: list[str], *, expect: int = 0) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != expect:
        raise AssertionError(
            f"command returned {completed.returncode}, expected {expect}: {command}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return completed


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparer", type=Path, required=True)
    parser.add_argument("--exporter", type=Path, required=True)
    parser.add_argument("--backend", type=Path, required=True)
    parser.add_argument("--reference-exporter", type=Path)
    args = parser.parse_args()
    preparer = str(args.preparer.resolve())
    exporter = str(args.exporter.resolve())
    backend = str(args.backend.resolve())

    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        protocol = root / "tiny-protocol.json"
        protocol.write_text(
            json.dumps(
                {
                    "schema": "koblitz_balanced_pdp_phase_a_protocol.v1",
                    "status": "frozen_before_preparation",
                    "master_seed_hex": SEED,
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

        plan = json.loads(run([preparer, "plan", "--protocol", str(protocol)]).stdout)
        assert plan["census_executed"] is False
        assert plan["requested_max_canonical_triples"] is None
        assert plan["selected_cell_ids"] == ["tiny-n7-standard"]
        assert plan["total_selected_instances"] == 4
        assert plan["total_planned_solver_runs"] == 12

        guarded = run(
            [
                preparer,
                "prepare",
                "--protocol",
                str(protocol),
                "--output",
                str(root / "guarded"),
                "--max-canonical-triples",
                "0",
            ],
            expect=2,
        )
        assert "above explicit ceiling" in guarded.stderr
        assert not (root / "guarded").exists()

        prepared = root / "prepared"
        run(
            [
                preparer,
                "prepare",
                "--protocol",
                str(protocol),
                "--output",
                str(prepared),
                "--max-canonical-triples",
                "1000",
            ]
        )
        blind_path = prepared / "blind/bundle.json"
        oracle_path = prepared / "sealed-oracle/oracle-ledger.json"
        blind_text = blind_path.read_text()
        for forbidden in FORBIDDEN_BLIND_TEXT:
            assert forbidden not in blind_text, f"blind bundle leaked {forbidden}"
        blind = json.loads(blind_text)
        oracle = json.loads(oracle_path.read_text())
        assert blind["instance_count"] == 4
        assert len({entry["blind_instance_id"] for entry in blind["instances"]}) == 4
        labels = [entry["target_class"] for entry in oracle["cells"][0]["entries"]]
        assert labels.count("decomposable") == 2
        assert labels.count("nondecomposable") == 2

        seal = json.loads((prepared / "seal.json").read_text())
        assert seal["blind_bundle_sha256"] == sha256(blind_path)
        assert seal["oracle_ledger_sha256"] == sha256(oracle_path)
        assert seal["requested_max_canonical_triples"] == 1000
        assert seal["selected_cell_ids"] == ["tiny-n7-standard"]
        assert len(seal["implementation"]["preparer_binary_sha256"]) == 64
        assert seal["execution_status"] == "not_run"

        instance = blind["instances"][0]
        target = instance["target"]
        unsafe_output = root / "unsafe-explicit-export"
        unsafe = run(
            [
                exporter,
                "7",
                "2",
                "standard",
                str(instance["source_nonce"]),
                "100000",
                str(unsafe_output),
                "0",
                "0",
                "--target-x",
                target["x"],
                "--target-y",
                target["y"],
                "--blind-instance-id",
                instance["blind_instance_id"],
            ],
            expect=101,
        )
        assert "require --export-only" in unsafe.stderr
        assert not unsafe_output.exists()

        invalid_output = root / "invalid-explicit-export"
        invalid = run(
            [
                exporter,
                "7",
                "2",
                "standard",
                str(instance["source_nonce"]),
                "100000",
                str(invalid_output),
                "0",
                "0",
                "--target-x",
                "0",
                "--target-y",
                "0",
                "--blind-instance-id",
                instance["blind_instance_id"],
                "--export-only",
            ],
            expect=101,
        )
        assert "not on the frozen curve" in invalid.stderr
        assert not invalid_output.exists()

        explicit_output = root / "explicit-export"
        run(
            [
                exporter,
                "7",
                "2",
                "standard",
                str(instance["source_nonce"]),
                "100000",
                str(explicit_output),
                "0",
                "0",
                "--target-x",
                target["x"],
                "--target-y",
                target["y"],
                "--blind-instance-id",
                instance["blind_instance_id"],
                "--export-only",
            ]
        )
        explicit = json.loads((explicit_output / "manifest.json").read_text())
        assert explicit["target"] == target
        assert explicit["seed"] == instance["source_nonce"]
        assert explicit["target_mode"] == "explicit_affine"
        assert explicit["blind_instance_id"] == instance["blind_instance_id"]
        assert explicit["source_instance"]["schema"] == "koblitz_pdp_source_instance.v2"
        assert (
            explicit["source_instance"]["identity"]["blind_instance_id"]
            == instance["blind_instance_id"]
        )
        assert "planted_points" not in explicit
        assert "planted_target_construction" not in explicit["timing_ns"]
        assert explicit["native_sat"]["status"] == "not_run_in_export_process"
        assert explicit["direct_meet_in_the_middle"]["status"] == "not_run_in_export_process"
        oracle_by_id = {
            entry["blind_instance_id"]: entry
            for entry in oracle["cells"][0]["entries"]
        }
        backend_receipt = json.loads(
            run([backend, "direct-mitm", str(explicit_output / "manifest.json")]).stdout
        )
        assert backend_receipt["source_instance_verified"] is True
        assert backend_receipt["regenerated_source_exact"] is True
        assert backend_receipt["source_instance_id"] == explicit["source_instance"]["id_blake3"]
        assert backend_receipt["status"] == (
            "sat"
            if oracle_by_id[instance["blind_instance_id"]]["target_class"] == "decomposable"
            else "unsat"
        )

        legacy_output = root / "legacy-export"
        run(
            [
                exporter,
                "7",
                "2",
                "standard",
                "1",
                "100000",
                str(legacy_output),
                "0",
                "0",
                "--export-only",
            ]
        )
        legacy = json.loads((legacy_output / "manifest.json").read_text())
        assert len(legacy["planted_points"]) == 3
        assert "blind_instance_id" not in legacy
        assert "planted_target_construction" in legacy["timing_ns"]
        legacy_backend = json.loads(
            run([backend, "direct-mitm", str(legacy_output / "manifest.json")]).stdout
        )
        assert legacy_backend["source_instance_verified"] is True
        assert legacy_backend["source_instance_id"] == legacy["source_instance"]["id_blake3"]

        if args.reference_exporter is not None:
            reference_output = root / "reference-legacy-export"
            run(
                [
                    str(args.reference_exporter.resolve()),
                    "7",
                    "2",
                    "standard",
                    "1",
                    "100000",
                    str(reference_output),
                    "0",
                    "0",
                    "--export-only",
                ]
            )
            reference = json.loads((reference_output / "manifest.json").read_text())
            legacy_without_timing = dict(legacy)
            reference_without_timing = dict(reference)
            legacy_without_timing.pop("timing_ns")
            reference_without_timing.pop("timing_ns")
            assert legacy_without_timing == reference_without_timing
            for filename in ("instance.anf", "instance.xor.cnf", "instance.magma"):
                assert (legacy_output / filename).read_bytes() == (
                    reference_output / filename
                ).read_bytes()

    print(
        json.dumps(
            {
                "schema": "koblitz_pdp_phase_a_tiny_test.v1",
                "status": "pass",
                "tiny_cells": 1,
                "prepared_instances": 4,
                "tiny_backend_validation_runs": 2,
                "legacy_reference_compared": args.reference_exporter is not None,
                "solver_runs": 0,
                "n59_census_runs": 0,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
