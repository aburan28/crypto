#!/usr/bin/env python3
"""Focused tests for the Stage-25 one-CPU receipt."""

from __future__ import annotations

from copy import deepcopy
import importlib.util
import math
from pathlib import Path
import unittest


SCRIPT = Path(__file__).with_name("run_koblitz_stage25_single_core.py")
SPEC = importlib.util.spec_from_file_location("stage25", SCRIPT)
assert SPEC and SPEC.loader
stage25 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage25)


class Stage25Tests(unittest.TestCase):
    def affinity(self) -> dict:
        return {
            "schema": stage25.AFFINITY_SCHEMA,
            "platform": "Linux",
            "uname": ["Linux", "host", "kernel", "version", "x86_64", ""],
            "initial_allowed_cpus": [4, 5],
            "selected_cpu": 4,
            "parent_effective_before": [4],
            "child_effective_before": [4],
            "parent_effective_after": [4],
            "child_effective_after": [4],
            "inheritance_contract": "Linux process affinity is inherited by fork/exec descendants unless explicitly changed",
        }

    def test_singleton_affinity_contract(self) -> None:
        self.assertEqual(stage25.validate_affinity(self.affinity())["selected_cpu"], 4)
        for field in (
            "parent_effective_before",
            "child_effective_before",
            "parent_effective_after",
            "child_effective_after",
        ):
            forged = self.affinity()
            forged[field] = [4, 5]
            with self.assertRaises(stage25.Stage25Error):
                stage25.validate_affinity(forged)

    def test_selected_cpu_must_be_initially_allowed(self) -> None:
        forged = self.affinity()
        forged["selected_cpu"] = 9
        with self.assertRaises(stage25.Stage25Error):
            stage25.validate_affinity(forged)

    def test_affinity_metadata_and_source_audit_fail_closed(self) -> None:
        for value in ([4, 4], [-1, 4], [5, 4]):
            forged = self.affinity()
            forged["initial_allowed_cpus"] = value
            with self.assertRaises(stage25.Stage25Error):
                stage25.validate_affinity(forged)
        forged = self.affinity()
        forged["uname"] = ["Linux"]
        with self.assertRaises(stage25.Stage25Error):
            stage25.validate_affinity(forged)
        self.assertEqual(stage25.affinity_reset_source_audit()["matches"], [])

    def test_result_rejects_claim_widening(self) -> None:
        result = {
            "schema": stage25.RESULT_SCHEMA,
            "status": "complete_single_cpu_affinity_control",
            "profile": "production",
            "source_commit": "0" * 40,
            "source_tree_sha256": "0" * 64,
            "affinity": self.affinity(),
            "stage23": {},
            "measurements": {
                "single_core_elapsed_seconds": 10.0,
                "outer_total_core_seconds": 9.0,
                "outer_user_seconds": 8.0,
                "outer_system_seconds": 1.0,
                "outer_peak_rss_bytes": 100,
                "stage23_child_total_core_seconds": 8.0,
                "stage23_summed_process_wall_seconds": 12.0,
                "stage23_peak_process_rss_bytes": 90,
                "ic_total_core_seconds": 7.0,
                "rho_total_core_seconds": 1.0,
                "online_ic_over_rho_core": 7.0,
                "setup_charged_ic_over_rho_core": 8.0,
                "affinity_wrapper_wall_seconds": 11.0,
                "affinity_wrapper_child_user_seconds": 8.0,
                "affinity_wrapper_child_system_seconds": 1.0,
            },
            "accounting_boundary": {
                "all_stage23_descendants_inherit_one_cpu_affinity": True,
                "single_core_elapsed_is_outer_wall_under_kernel_affinity": True,
                "simultaneous_process_tree_memory_measured": False,
                "source_dependency_toolchain_or_license_acquisition_charged": False,
                "covers_stage20_n31_n41_n59_backend_matrix": False,
                "licensed_magma_executed": False,
            },
            "retained_mathematical_witness_replay_completed": True,
            "independent_mathematical_payload_replay_completed": False,
            "scientific_measurement_admitted": False,
            "independent_external_reproduction_satisfied": False,
            "full_cost_gate_passed": False,
            "koblitz_index_calculus_sota": False,
        }
        stage25.validate_result(result)
        for field in (
            "independent_mathematical_payload_replay_completed",
            "scientific_measurement_admitted",
            "independent_external_reproduction_satisfied",
            "full_cost_gate_passed",
            "koblitz_index_calculus_sota",
        ):
            forged = deepcopy(result)
            forged[field] = True
            with self.assertRaises(stage25.Stage25Error):
                stage25.validate_result(forged)
        for invalid in (math.inf, math.nan):
            forged = deepcopy(result)
            forged["measurements"]["single_core_elapsed_seconds"] = invalid
            with self.assertRaises(stage25.Stage25Error):
                stage25.validate_result(forged)


if __name__ == "__main__":
    unittest.main()
