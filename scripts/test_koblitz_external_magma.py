#!/usr/bin/env python3
"""Focused custody and orchestration tests for the licensed-host Magma runner."""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import tempfile
import textwrap
import unittest

import run_koblitz_external_magma as external


FIRST_TASK = "seed-2026091301/n31-l5-m3-standard-a1-f0"


def executable(path: Path, source: str) -> Path:
    path.write_text(textwrap.dedent(source).lstrip())
    path.chmod(0o755)
    return path


def fake_tools(
    root: Path,
    *,
    always_nonlifting: bool = False,
    duplicate_models: bool = False,
    exhaust_after_nonlifting: bool = False,
):
    magma = executable(
        root / "magma",
        f"""
        #!/usr/bin/env python3
        import os
        import pathlib
        import re
        import sys

        if "--version" in sys.argv:
            print("Magma V2.fake")
            raise SystemExit(0)
        script = pathlib.Path(sys.argv[-1]).read_text()
        print("KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1")
        print("KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse")
        f4_status = os.environ.get("FAKE_F4_STATUS", "SAT")
        print("KOBLITZ_MAGMA_STATUS=" + f4_status)
        print("KOBLITZ_MAGMA_F4_DEGREES=" + ("[]" if f4_status == "UNSAT" else "[ 3, 4 ]"))
        print("KOBLITZ_MAGMA_BASIS_SIZE=" + ("1" if f4_status == "UNSAT" else "2"))
        print("KOBLITZ_MAGMA_CPU_SECONDS=0.01")
        print("KOBLITZ_MAGMA_WALL_SECONDS=0.02")
        if "KOBLITZ_MODEL_SCHEMA" in script:
            count = int(re.search(r"KOBLITZ_MODEL_EXCLUDED_COUNT=(\\d+)", script).group(1))
            variables = int(re.search(r"KOBLITZ_MODEL_VARIABLES=(\\d+)", script).group(1))
            number = 0 if {duplicate_models!r} else count
            bits = format(number, f"0{{variables}}b")[-variables:]
            print("KOBLITZ_MODEL_SCHEMA=koblitz_magma_f4_model.v1")
            exhausted = {exhaust_after_nonlifting!r} and count > 0
            print("KOBLITZ_MODEL_STATUS=" + ("UNSAT" if exhausted else "SAT"))
            print("KOBLITZ_MODEL_BITS=" + ("-" if exhausted else bits))
            print(f"KOBLITZ_MODEL_VARIABLES={{variables}}")
            print(f"KOBLITZ_MODEL_EXCLUDED_COUNT={{count}}")
            print("KOBLITZ_MODEL_CPU_SECONDS=0.003")
            print("KOBLITZ_MODEL_WALL_SECONDS=0.004")
        """,
    )
    backend = executable(
        root / "backend",
        f"""
        #!/usr/bin/env python3
        import json
        import pathlib
        import sys

        if "--version" in sys.argv:
            print("fake koblitz_pdp_backend 1")
            raise SystemExit(0)
        manifest = json.loads(pathlib.Path(sys.argv[2]).read_text())
        assignment_path = pathlib.Path(sys.argv[3])
        assignment = json.loads(assignment_path.read_text())
        valid = False if {always_nonlifting!r} else any(assignment)
        report = {{
            "schema": "koblitz_pdp_assignment_validation.v1",
            "status": "valid_point_witness" if valid else "nonlifting_source_model",
            "source_instance_id": manifest["source_instance"]["id_blake3"],
            "source_instance_verified": True,
            "regenerated_source_exact": True,
            "assignment_values": len(assignment),
            "assignment_blake3": "a" * 64,
            "source_assignment": assignment,
            "source_model_valid": True,
            "source_witness_valid": valid,
        }}
        print(json.dumps(report))
        raise SystemExit(0 if valid else 2)
        """,
    )
    minisat = executable(
        root / "minisat",
        """
        #!/usr/bin/env python3
        print("fake minisat 2")
        """,
    )
    return magma, backend, minisat


def arguments(output: Path, magma: Path, backend: Path, minisat: Path, model_cap: int = 3):
    return SimpleNamespace(
        protocol=external.DEFAULT_PROTOCOL,
        panel=external.DEFAULT_PANEL,
        meter=external.DEFAULT_METER,
        magma=str(magma),
        backend=backend,
        minisat=str(minisat),
        output=output,
        resume=False,
        only_task=[FIRST_TASK],
        f4_timeout=10.0,
        model_timeout=10.0,
        validation_timeout=10.0,
        model_cap=model_cap,
        synthetic_test_mode=True,
    )


class FrozenManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.packet = external.build_dry_run_manifest(
            external.DEFAULT_PROTOCOL,
            external.DEFAULT_PANEL,
            external.DEFAULT_METER,
        )

    def test_exact_twenty_task_seed_major_inventory(self) -> None:
        self.assertEqual(self.packet["schema"], external.MANIFEST_SCHEMA)
        self.assertEqual(self.packet["expected_tasks"], 20)
        self.assertEqual(len(self.packet["tasks"]), 20)
        self.assertEqual(self.packet["tasks"][0]["id"], FIRST_TASK)
        self.assertEqual(
            self.packet["tasks"][-1]["id"],
            "seed-2026091305/n59-l9-m3-standard-a1-f0",
        )
        self.assertNotEqual(
            self.packet["archive_commit"], self.packet["panel_source_revision"]
        )
        self.assertEqual(self.packet["archive_commit"], external.STAGE13_ARCHIVE_COMMIT)
        self.assertEqual(
            self.packet["panel_source_revision"],
            "67cf7f559cf6cf4c70fc16f22eb5dd68e4883d31",
        )
        self.assertTrue(
            self.packet["contract_pins"]["archive_tree_matches_pinned_commit"]
        )
        self.assertEqual(
            self.packet["contract_pins"]["stage13_verifier"]["sha256"],
            external.PINNED_STAGE13_VERIFIER_SHA256,
        )
        self.assertEqual(
            self.packet["contract_pins"]["matrix_contract"]["sha256"],
            external.PINNED_MATRIX_CONTRACT_SHA256,
        )
        self.assertEqual(len({task["source_instance_sha256"] for task in self.packet["tasks"]}), 20)

    def test_all_frozen_magma_inputs_have_exact_direct_f4_suffix(self) -> None:
        panel = Path(self.packet["archive_root"])
        for task in self.packet["tasks"]:
            path = external.resolve_manifest_path(
                panel, task["files"]["magma_boolean_f4"], "Magma input"
            )
            prefix = external.split_magma_source(path)
            self.assertTrue(prefix.endswith("I := ideal<R | F>;\n"))

    def test_changed_magma_template_is_rejected(self) -> None:
        first = self.packet["tasks"][0]
        panel = Path(self.packet["archive_root"])
        source = external.resolve_manifest_path(
            panel, first["files"]["magma_boolean_f4"], "Magma input"
        )
        with tempfile.TemporaryDirectory() as directory:
            changed = Path(directory) / "instance.magma"
            changed.write_bytes(source.read_bytes() + b"// changed\n")
            with self.assertRaises(external.ExternalMagmaError):
                external.split_magma_source(changed)


class ModelContractTests(unittest.TestCase):
    def test_orphan_group_cleanup_is_not_resource_complete(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            parent = executable(
                root / "leave-child",
                f"""
                #!{sys.executable}
                import subprocess
                import sys
                subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
                """,
            )
            receipt = external.run_metered(
                external.DEFAULT_METER,
                [str(parent)],
                root,
                5.0,
                root / "attempt",
            )
            self.assertTrue(receipt["process"]["orphan_group_terminated"])
            self.assertFalse(external.process_receipt_resource_complete(receipt))

    def test_binary_identity_falls_back_after_failed_nonempty_probe(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            tool = executable(
                Path(directory) / "solver",
                """
                #!/usr/bin/env python3
                import sys
                if "--version" in sys.argv:
                    print("unsupported --version")
                    raise SystemExit(2)
                if "-h" in sys.argv:
                    print("solver help version 1")
                    raise SystemExit(0)
                """,
            )
            identity = external.binary_identity(tool, [["--version"], ["-h"]])
            self.assertEqual(identity["version"], "solver help version 1")
            self.assertEqual(identity["version_command"][-1], "-h")
            self.assertEqual(identity["version_returncode"], 0)
            self.assertEqual(identity["failed_version_probes"][0]["returncode"], 2)

    def test_recognizes_only_documented_magma_version_forms(self) -> None:
        self.assertTrue(external.recognized_magma_version("Magma V2.29-10 Linux"))
        self.assertTrue(external.recognized_magma_version("2.29-10\n"))
        self.assertFalse(external.recognized_magma_version("Magma V2.fake"))
        self.assertFalse(external.recognized_magma_version("wrapper 2.29-10"))

    def test_strict_model_terminal(self) -> None:
        valid = "\n".join(
            [
                "KOBLITZ_MODEL_SCHEMA=koblitz_magma_f4_model.v1",
                "KOBLITZ_MODEL_STATUS=SAT",
                "KOBLITZ_MODEL_BITS=0101",
                "KOBLITZ_MODEL_VARIABLES=4",
                "KOBLITZ_MODEL_EXCLUDED_COUNT=2",
                "KOBLITZ_MODEL_CPU_SECONDS=1.0",
                "KOBLITZ_MODEL_WALL_SECONDS=2.0",
            ]
        )
        parsed = external.parse_model_terminal(valid, 4, 2)
        self.assertEqual(parsed["assignment"], [False, True, False, True])
        self.assertIsNone(external.parse_model_terminal(valid.replace("0101", "0121"), 4, 2))
        self.assertIsNone(external.parse_model_terminal(valid + "\nKOBLITZ_MODEL_STATUS=SAT", 4, 2))
        self.assertIsNone(external.parse_model_terminal(valid.replace("1.0", "nan"), 4, 2))
        self.assertIsNone(external.parse_model_terminal(valid, 4, 1))

    def test_witness_script_records_exclusions_and_recomputes_f4(self) -> None:
        packet = external.build_dry_run_manifest(
            external.DEFAULT_PROTOCOL, external.DEFAULT_PANEL, external.DEFAULT_METER
        )
        task = packet["tasks"][0]
        path = external.resolve_manifest_path(
            Path(packet["archive_root"]), task["files"]["magma_boolean_f4"], "Magma input"
        )
        excluded = [[False] * task["source_variables"]]
        script = external.render_witness_script(path, excluded, task["source_variables"])
        self.assertIn('GroebnerBasis(I : Al := "Direct"', script)
        self.assertIn("SAT(G : Exclude := Excluded, Verbose := false)", script)
        self.assertIn("KOBLITZ_MODEL_EXCLUDED_COUNT=1", script)


class FakeLicensedHostTests(unittest.TestCase):
    def run_fake(
        self,
        *,
        always_nonlifting=False,
        duplicate_models=False,
        exhaust_after_nonlifting=False,
        model_cap=3,
    ):
        temporary = tempfile.TemporaryDirectory()
        root = Path(temporary.name)
        tools = root / "tools"
        tools.mkdir()
        magma, backend, minisat = fake_tools(
            tools,
            always_nonlifting=always_nonlifting,
            duplicate_models=duplicate_models,
            exhaust_after_nonlifting=exhaust_after_nonlifting,
        )
        args = arguments(root / "run", magma, backend, minisat, model_cap=model_cap)
        summary = external.execute(args)
        packet = external.build_dry_run_manifest(
            external.DEFAULT_PROTOCOL, external.DEFAULT_PANEL, external.DEFAULT_METER
        )
        progress = external.read_json(args.output / "run.json")
        task = packet["tasks"][0]
        task_path = args.output / "tasks" / task["id"] / "task-result.json"
        result = external.read_json(task_path)
        return temporary, args, summary, packet, progress, task, task_path, result

    def test_nonlifting_model_is_excluded_then_valid_model_is_charged(self) -> None:
        temporary, args, summary, packet, progress, task, task_path, result = self.run_fake()
        self.addCleanup(temporary.cleanup)
        self.assertEqual(result["final_status"], "sat")
        self.assertEqual(result["valid_point_witness_ordinal"], 2)
        self.assertEqual(len(result["model_attempts"]), 2)
        self.assertEqual(
            result["model_attempts"][0]["validation"]["status"],
            "nonlifting_source_model",
        )
        self.assertEqual(
            result["model_attempts"][1]["validation"]["status"],
            "valid_point_witness",
        )
        self.assertFalse(summary["full_gate_passed"])
        self.assertEqual(summary["charged_totals"]["primary_f4"]["processes"], 1)
        self.assertEqual(summary["charged_totals"]["f4_plus_sat_g"]["processes"], 2)
        self.assertEqual(summary["charged_totals"]["validator"]["processes"], 2)
        self.assertEqual(summary["charged_totals"]["all_metered_processes"]["processes"], 5)
        self.assertEqual(len(summary["per_cell_distributions"]), 1)
        self.assertGreater(summary["charged_totals"]["all_metered_processes"]["total_core_seconds"], 0)
        tools = progress["tools"]
        policy = progress["policy"]
        validated = external.validate_completed_task(task_path, task, packet, tools, policy)
        self.assertEqual(validated["final_status"], "sat")

        synthetic_full_results = []
        for frozen_task in packet["tasks"]:
            row = copy.deepcopy(result)
            row["id"] = frozen_task["id"]
            row["seed"] = frozen_task["seed"]
            row["cell"] = frozen_task["cell"]
            synthetic_full_results.append(row)
        synthetic_full = external.summarize_run(
            packet, packet["tasks"], synthetic_full_results, policy
        )
        self.assertTrue(synthetic_full["full_selection"])
        self.assertTrue(synthetic_full["all_f4_terminals_resource_complete"])
        self.assertTrue(synthetic_full["all_sat_point_witnesses_validated"])
        self.assertTrue(synthetic_full["synthetic_test_mode"])
        self.assertFalse(synthetic_full["full_gate_passed"])
        dirty_real_policy = copy.deepcopy(policy)
        dirty_real_policy["synthetic_test_mode"] = False
        dirty_real_policy["relevant_checkout_clean"] = False
        dirty_full = external.summarize_run(
            packet, packet["tasks"], synthetic_full_results, dirty_real_policy
        )
        self.assertTrue(dirty_full["archive_contracts_pinned"])
        self.assertFalse(dirty_full["relevant_checkout_clean"])
        self.assertFalse(dirty_full["full_gate_passed"])

        args.resume = True
        before = external.sha256_file(task_path)
        external.execute(args)
        self.assertEqual(external.sha256_file(task_path), before)

    def test_fake_magma_banner_is_rejected_outside_synthetic_mode(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tools = root / "tools"
            tools.mkdir()
            magma, backend, minisat = fake_tools(tools)
            args = arguments(root / "run", magma, backend, minisat)
            args.synthetic_test_mode = False
            with self.assertRaisesRegex(
                external.ExternalMagmaError, "Magma --version must contain"
            ):
                external.execute(args)

    def test_output_inside_frozen_panel_is_rejected_before_creation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            tools = Path(directory) / "tools"
            tools.mkdir()
            magma, backend, minisat = fake_tools(tools)
            forbidden = external.DEFAULT_PANEL / "external-magma-output-must-not-exist"
            args = arguments(forbidden, magma, backend, minisat)
            with self.assertRaisesRegex(
                external.ExternalMagmaError, "inside the frozen Stage 13 panel"
            ):
                external.execute(args)
            self.assertFalse(forbidden.exists())

    def test_model_cap_is_inconclusive_after_every_attempt_is_metered(self) -> None:
        temporary, _, summary, _, _, _, _, result = self.run_fake(
            always_nonlifting=True, model_cap=2
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(result["final_status"], "model_cap_inconclusive")
        self.assertEqual(len(result["model_attempts"]), 2)
        self.assertEqual(summary["charged_totals"]["f4_plus_sat_g"]["processes"], 2)
        self.assertEqual(summary["charged_totals"]["validator"]["processes"], 2)
        self.assertFalse(summary["full_gate_passed"])

    def test_resume_recovers_and_charges_an_interrupted_process_receipt(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        tool_root = root / "tools"
        tool_root.mkdir()
        magma_path, backend_path, minisat_path = fake_tools(tool_root)
        args = arguments(root / "run", magma_path, backend_path, minisat_path)
        packet = external.build_dry_run_manifest(
            args.protocol, args.panel, args.meter
        )
        task = packet["tasks"][0]
        magma = external.binary_identity(args.magma, [["--version"]])
        magma["recognized_version_banner"] = False
        backend = external.binary_identity(args.backend, [["--version"]])
        minisat = external.binary_identity(args.minisat, [["--version"], ["-h"]])
        tools = {
            "magma": magma,
            "backend": backend,
            "minisat": minisat,
            "meter": {
                "path": str(args.meter.resolve()),
                "sha256": external.sha256_file(args.meter),
            },
            "runner": {
                "path": str(Path(external.__file__).resolve()),
                "sha256": external.sha256_file(Path(external.__file__)),
            },
            "matrix_contract": {
                "path": str((external.HERE / "run_koblitz_pdp_matrix.py").resolve()),
                "sha256": external.sha256_file(
                    external.HERE / "run_koblitz_pdp_matrix.py"
                ),
            },
            "stage13_verifier": {
                "path": str(
                    (external.STAGE / "verify_stage13_pdp_panel.py").resolve()
                ),
                "sha256": external.sha256_file(
                    external.STAGE / "verify_stage13_pdp_panel.py"
                ),
                "expected_sha256": external.PINNED_STAGE13_VERIFIER_SHA256,
            },
        }
        policy = {
            "parallel_workers": 1,
            "magma_threads": 1,
            "gpu_disabled_by_input": True,
            "f4_watchdog_seconds": args.f4_timeout,
            "model_watchdog_seconds": args.model_timeout,
            "validation_watchdog_seconds": args.validation_timeout,
            "model_cap": args.model_cap,
            "unknown_and_timeout_are_inconclusive": True,
            "witness_process_includes_f4_recomputation": True,
            "synthetic_test_mode": True,
            "tool_identity_bound": True,
            "license_entitlement_authenticated": False,
            "archive_contracts_pinned": True,
            "relevant_checkout_clean": packet["relevant_checkout"]["dirty"] is False,
            "relevant_checkout_state": packet["relevant_checkout"],
        }
        args.output.mkdir(parents=True)
        external.atomic_json(
            args.output / "run.json",
            {
                "schema": external.RUN_SCHEMA,
                "dry_run_manifest_sha256": packet["manifest_sha256"],
                "selection": [task["id"]],
                "tools": tools,
                "host": external.host_identity(),
                "policy": policy,
                "started_at": external.now(),
                "status": "running",
                "tasks": {},
            },
        )
        task_root = args.output / "tasks" / task["id"]
        task_root.mkdir(parents=True)
        frozen_magma = external.resolve_manifest_path(
            Path(packet["archive_root"]),
            task["files"]["magma_boolean_f4"],
            "Magma input",
        )
        external.run_metered(
            args.meter,
            [magma["path"], "-t", "1", "-b", str(frozen_magma)],
            task_root,
            args.f4_timeout,
            task_root / "f4-attempt-001",
        )
        args.resume = True
        summary = external.execute(args)
        result = external.read_json(task_root / "task-result.json")
        self.assertEqual(len(result["interrupted_processes"]), 1)
        self.assertEqual(
            result["interrupted_processes"][0]["kind"], "primary_f4"
        )
        self.assertEqual(
            summary["charged_totals"]["interrupted_prior"]["processes"], 1
        )
        self.assertEqual(
            summary["charged_totals"]["all_metered_processes"]["processes"], 6
        )
        self.assertFalse(result["incomplete_prior_attempts"])
        self.assertFalse(result["contradictions"])

    def test_interrupted_sat_model_without_validator_is_incomplete(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tools = root / "tools"
            tools.mkdir()
            magma_path, backend_path, _ = fake_tools(tools)
            magma = external.binary_identity(magma_path, [["--version"]])
            backend = external.binary_identity(backend_path, [["--version"]])
            packet = external.build_dry_run_manifest(
                external.DEFAULT_PROTOCOL,
                external.DEFAULT_PANEL,
                external.DEFAULT_METER,
            )
            task = packet["tasks"][0]
            panel = Path(packet["archive_root"])
            frozen_magma = external.resolve_manifest_path(
                panel, task["files"]["magma_boolean_f4"], "Magma input"
            )
            manifest = external.resolve_manifest_path(
                panel, task["files"]["manifest"], "manifest"
            )
            task_root = root / "task"
            attempt = task_root / "model-attempt-001"
            attempt.mkdir(parents=True)
            witness = attempt / "witness.magma"
            witness.write_text(
                external.render_witness_script(
                    frozen_magma, [], task["source_variables"]
                )
            )
            external.run_metered(
                external.DEFAULT_METER,
                [magma["path"], "-t", "1", "-b", str(witness.resolve())],
                attempt,
                10.0,
                attempt / "magma-process",
            )
            recovered, incomplete = external.collect_interrupted_processes(
                task_root,
                external.DEFAULT_METER,
                magma,
                backend,
                frozen_magma,
                manifest,
            )
            self.assertEqual(len(recovered), 1)
            self.assertEqual(recovered[0]["kind"], "f4_plus_sat_g")
            self.assertEqual(incomplete, [str(attempt.resolve())])

    def test_recovered_unsat_terminal_is_a_charged_planted_source_contradiction(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tools = root / "tools"
            tools.mkdir()
            magma_path, backend_path, minisat_path = fake_tools(tools)
            magma = external.binary_identity(magma_path, [["--version"]])
            backend = external.binary_identity(backend_path, [["--version"]])
            packet = external.build_dry_run_manifest(
                external.DEFAULT_PROTOCOL,
                external.DEFAULT_PANEL,
                external.DEFAULT_METER,
            )
            task = packet["tasks"][0]
            panel = Path(packet["archive_root"])
            frozen_magma = external.resolve_manifest_path(
                panel, task["files"]["magma_boolean_f4"], "Magma input"
            )
            output = root / "run"
            task_root = output / "tasks" / task["id"]
            task_root.mkdir(parents=True)
            prior_environment = os.environ.copy()
            prior_environment["FAKE_F4_STATUS"] = "UNSAT"
            external.run_metered(
                external.DEFAULT_METER,
                [magma["path"], "-t", "1", "-b", str(frozen_magma)],
                task_root,
                10.0,
                task_root / "f4-attempt-001",
                prior_environment,
            )
            environment = os.environ.copy()
            environment["PATH"] = str(minisat_path.parent) + os.pathsep + environment.get(
                "PATH", ""
            )
            policy = {
                "parallel_workers": 1,
                "magma_threads": 1,
                "gpu_disabled_by_input": True,
                "f4_watchdog_seconds": 10.0,
                "model_watchdog_seconds": 10.0,
                "validation_watchdog_seconds": 10.0,
                "model_cap": 3,
                "unknown_and_timeout_are_inconclusive": True,
                "witness_process_includes_f4_recomputation": True,
                "synthetic_test_mode": True,
                "tool_identity_bound": True,
                "license_entitlement_authenticated": False,
                "archive_contracts_pinned": True,
                "relevant_checkout_clean": packet["relevant_checkout"]["dirty"] is False,
                "relevant_checkout_state": packet["relevant_checkout"],
            }
            result = external.execute_task(
                packet,
                task,
                output,
                external.DEFAULT_METER,
                magma,
                backend,
                environment,
                policy,
            )
            self.assertEqual(result["interrupted_processes"][0]["f4_status"], "unsat")
            self.assertIn(
                "an interrupted clean Magma receipt reported UNSAT for a frozen planted instance",
                result["contradictions"],
            )
            self.assertFalse(
                external.summarize_run(packet, [task], [result], policy)[
                    "no_contradictions"
                ]
            )

    def test_repeated_excluded_model_is_rejected(self) -> None:
        temporary, _, summary, _, _, _, _, result = self.run_fake(
            always_nonlifting=True, duplicate_models=True, model_cap=3
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(result["final_status"], "duplicate_excluded_model")
        self.assertTrue(result["contradictions"])
        self.assertEqual(len(result["model_attempts"]), 2)
        self.assertEqual(summary["charged_totals"]["validator"]["processes"], 1)
        self.assertFalse(summary["no_contradictions"])

    def test_exhaustion_after_nonlifting_model_is_a_planted_source_contradiction(self) -> None:
        temporary, _, summary, _, _, _, _, result = self.run_fake(
            always_nonlifting=True,
            exhaust_after_nonlifting=True,
            model_cap=3,
        )
        self.addCleanup(temporary.cleanup)
        self.assertEqual(
            result["final_status"], "nonlifting_models_exhausted_contradiction"
        )
        self.assertTrue(result["contradictions"])
        self.assertFalse(summary["no_contradictions"])
        self.assertFalse(summary["full_gate_passed"])

    def test_resume_rejects_command_output_metric_assignment_and_validator_tampering(self) -> None:
        temporary, _, _, packet, progress, task, task_path, result = self.run_fake()
        self.addCleanup(temporary.cleanup)
        tools = progress["tools"]
        policy = progress["policy"]

        def assert_tamper(path: Path, replacement: bytes) -> None:
            original = path.read_bytes()
            path.write_bytes(replacement)
            with self.assertRaises(external.ExternalMagmaError):
                external.validate_completed_task(task_path, task, packet, tools, policy)
            path.write_bytes(original)
            external.validate_completed_task(task_path, task, packet, tools, policy)

        f4_stdout = Path(result["f4"]["receipt"]["stdout"]["path"])
        assert_tamper(f4_stdout, f4_stdout.read_bytes() + b"tamper\n")

        f4_metrics = f4_stdout.parent / "metrics.json"
        changed_metrics = json.loads(f4_metrics.read_text())
        changed_metrics["metrics"]["wall_seconds"] += 1
        assert_tamper(f4_metrics, (json.dumps(changed_metrics) + "\n").encode())

        assignment = Path(result["model_attempts"][0]["assignment"]["path"])
        assert_tamper(assignment, b"[true]\n")

        validation_stdout = Path(
            result["model_attempts"][0]["validation"]["receipt"]["stdout"]["path"]
        )
        assert_tamper(validation_stdout, validation_stdout.read_bytes() + b"{}\n")

        original_task = task_path.read_bytes()
        changed_task = json.loads(original_task)
        changed_task["f4"]["receipt"]["process"]["command"][0] = "/wrong/magma"
        task_path.write_text(json.dumps(changed_task))
        with self.assertRaises(external.ExternalMagmaError):
            external.validate_completed_task(task_path, task, packet, tools, policy)
        task_path.write_bytes(original_task)
        external.validate_completed_task(task_path, task, packet, tools, policy)


if __name__ == "__main__":
    unittest.main()
