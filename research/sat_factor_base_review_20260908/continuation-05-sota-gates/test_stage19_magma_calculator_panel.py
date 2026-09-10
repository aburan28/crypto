#!/usr/bin/env python3
"""Adversarial and crash-boundary tests for the Stage 19 calculator panel."""

from __future__ import annotations

from contextlib import ExitStack
from datetime import datetime, timedelta, timezone
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent


def load_module(name: str, filename: str):
    path = HERE / filename
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RENDERER = load_module("stage19_renderer_tests", "render_stage19_magma_calculator_panel.py")
RUNNER = load_module("stage19_runner_tests", "run_stage19_magma_calculator_panel.py")
VERIFIER = load_module("stage19_verifier_tests", "verify_stage19_magma_calculator_panel.py")
CHILD = RUNNER.CHILD
PREPARER = load_module("stage19_manifest_tests", "prepare_stage19_magma_execution_manifest.py")
AMENDMENT = load_module("stage19_amendment_tests", "verify_stage19_amendment01.py")
AMENDMENT02 = load_module("stage19_amendment02_tests", "verify_stage19_amendment02.py")
TLS_PROBE = load_module("stage19_tls_probe_tests", "probe_stage19_magma_tls.py")


def stamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def fresh_probe_fixture() -> dict:
    archived = TLS_PROBE.verify_archived()
    receipt = {
        "schema": TLS_PROBE.FRESH_SCHEMA,
        "started_at": "2026-09-10T14:00:00Z",
        "finished_at": "2026-09-10T14:00:00.125000Z",
        "wall_seconds": 0.125,
        "request": archived["request"],
        "response": {
            **archived["response"],
            "selected_headers_without_date": {},
        },
        "fresh_immediately_before_manifest_generation": True,
        "classification": archived["classification"],
        "claim_boundary": archived["claim_boundary"],
    }
    receipt["receipt_sha256"] = TLS_PROBE.sha256_bytes(
        TLS_PROBE.canonical_bytes(receipt)
    )
    return receipt


class Stage19PanelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cached_plan, cls.cached_inputs = RENDERER.build_plan()

    def fake_binding(self) -> dict:
        body = {
            "schema": "koblitz_magma_calculator_stage19_execution_binding.v1",
            "first_preexecution_commit": "1" * 40,
            "first_preexecution_tree": "2" * 40,
            "execution_commit": "3" * 40,
            "execution_tree": "4" * 40,
            "tree_delta": [{"status": "A", "path": "execution-manifest.json"}],
            "execution_manifest": {
                "path": "execution-manifest.json",
                "bytes": 1,
                "sha256": "5" * 64,
                "git_mode": "100644",
                "git_blob_oid": "6" * 40,
            },
            "relevant_blob_count": 21,
            "relevant_path_list_sha256": "7" * 64,
            "ca_bundle": CHILD.expected_ca_bundle_record(),
            "fresh_preexecution_tls_get_probe": {"receipt_sha256": "8" * 64},
            "checkout_clean_outside_runtime_artifacts": True,
        }
        body["binding_sha256"] = RUNNER.sha256_bytes(RUNNER.canonical_bytes(body))
        return body

    def materialized(self) -> tuple[Path, dict, dict]:
        temporary = tempfile.TemporaryDirectory()
        artifact = Path(temporary.name) / "artifact"
        plan = self.cached_plan
        CHILD.atomic_write(artifact / "plan.json", RENDERER.canonical_bytes(plan))
        for filename, data in self.cached_inputs.items():
            CHILD.atomic_write(artifact / "inputs" / filename, data)
        binding = self.fake_binding()
        prepared = VERIFIER.summary_body(plan, [], "prepared_not_executed", None)
        CHILD.atomic_write(artifact / "prepared-summary.json", VERIFIER.canonical_bytes(prepared))
        self.addCleanup(temporary.cleanup)
        return artifact, plan, binding

    def initialize_run(self, artifact: Path, plan: dict, binding: dict) -> dict:
        state = {
            "schema": RUNNER.RUN_SCHEMA,
            "plan_sha256": RUNNER.sha256_bytes(RUNNER.canonical_bytes(plan)),
            "protocol_sha256": plan["protocol"]["sha256"],
            "task_order": [task["id"] for task in plan["tasks"]],
            "execution_policy": plan["execution_policy"],
            "execution_binding": binding,
            "started_at": stamp(datetime.now(timezone.utc) - timedelta(seconds=10)),
            "status": "running",
            "receipts": {},
        }
        RUNNER.atomic_json(artifact / "run.json", state)
        return state

    def response_xml(self, task: dict, status: str = "SAT", duplicate_identity: bool = False) -> bytes:
        identity = [
            f"KOBLITZ_MAGMA_TASK_ID={task['id']}",
            f"KOBLITZ_MAGMA_SOURCE_SHA256={task['source_instance_sha256']}",
        ]
        if duplicate_identity:
            identity.append(f"KOBLITZ_MAGMA_TASK_ID={task['id']}")
        basis = 1 if status == "UNSAT" else 43
        lines = identity + [
            "KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1",
            "KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse",
            f"KOBLITZ_MAGMA_STATUS={status}",
            "KOBLITZ_MAGMA_F4_DEGREES=[ 2, 3, 3, 3, 3 ]",
            f"KOBLITZ_MAGMA_BASIS_SIZE={basis}",
            "KOBLITZ_MAGMA_CPU_SECONDS=0.1",
            "KOBLITZ_MAGMA_WALL_SECONDS=0.1",
            "",
        ]
        results = "".join(f"<line>{line}</line>" for line in lines)
        return (
            '<?xml version="1.0"?>\n<calculator><headers><max_time>60</max_time>'
            '<max_input>50000</max_input><seed>1</seed><version>2.29-10</version>'
            '<time>0.1</time><memory>32.09MB</memory></headers><results>'
            + results + "</results></calculator>\n"
        ).encode()

    def fake_process(
        self,
        artifact: Path,
        attempt: Path,
        plan: dict,
        task: dict,
        kind: str,
        *,
        timed_out: bool = False,
        orphan: bool = False,
        returncode: int = 0,
    ) -> dict:
        input_bytes = (artifact / task["named_input"]["path"]).read_bytes()
        request_body = RUNNER.request_body(input_bytes, plan["service"]["form_field"])
        if kind == "clean":
            body, http_status, final_url, error, exceeded = (
                self.response_xml(task), 200, plan["service"]["endpoint"], None, False
            )
        elif kind == "wrong_url":
            body, http_status, final_url, error, exceeded = (
                self.response_xml(task), 200, "https://example.invalid/", None, False
            )
        elif kind == "body_limit":
            body, http_status, final_url, error, exceeded = (
                b"x" * 32, 200, plan["service"]["endpoint"], None, True
            )
        else:
            body, http_status, final_url, error, exceeded = (
                None, None, None, "TimeoutError: synthetic", False
            )
        result = {
            "http_status": http_status,
            "final_url": final_url,
            "headers": {"content-type": "text/xml"} if body is not None else {},
            "body": body,
            "body_limit_exceeded": exceeded,
            "transport_error": error,
            "request_body_bytes": len(request_body),
            "request_body_sha256": RUNNER.sha256_bytes(request_body),
        }
        authorization = attempt / "launch-authorization.json"
        consumed = attempt / "launch-authorization.consumed.json"
        if authorization.is_file():
            os.replace(authorization, consumed)
        CHILD.persist_result(attempt, input_bytes, result, RUNNER.now())
        command = VERIFIER.expected_transport_command(artifact, task)
        CHILD.atomic_write(attempt / "transport.stdout", b"{}\n")
        CHILD.atomic_write(attempt / "transport.stderr", b"")
        metrics = {
            "command": command,
            "returncode": returncode,
            "watchdog_seconds": 75.0,
            "timed_out": timed_out,
            "orphan_group_terminated": orphan,
            "metrics": {
                "wall_seconds": 0.1,
                "user_seconds": 0.01,
                "system_seconds": 0.01,
                "total_core_seconds": 0.02,
                "single_core_seconds": 0.02,
                "peak_rss_bytes": 1024,
                "meter": "fresh-process getrusage(RUSAGE_CHILDREN)",
            },
        }
        RUNNER.atomic_json(attempt / "transport-metrics.json", metrics)
        return {
            "metrics": metrics,
            "stdout": RUNNER.regular_record(attempt / "transport.stdout", attempt),
            "stderr": RUNNER.regular_record(attempt / "transport.stderr", attempt),
            "metrics_file": RUNNER.regular_record(attempt / "transport-metrics.json", attempt),
        }

    def assert_recovery_boundary(self, boundary: str, populate) -> None:  # noqa: ANN001
        artifact, plan, binding = self.materialized()
        self.initialize_run(artifact, plan, binding)
        task = plan["tasks"][0]
        attempt = artifact / "attempts" / RUNNER.task_stem(task)
        attempt.mkdir(parents=True)
        populate(artifact, plan, binding, task, attempt)
        with mock.patch.object(RUNNER, "run_metered_transport") as post:
            state = RUNNER.run(artifact, resume=True, execution_binding_override=binding)
        post.assert_not_called()
        self.assertEqual("halted_after_nonclean_receipt", state["status"])
        receipt = RUNNER.read_json(attempt / "receipt.json")
        self.assertEqual(boundary, receipt["outcome"])
        self.assertFalse(receipt["claim_admitted"])
        verified = VERIFIER.verify(artifact, execution_binding_override=binding)
        self.assertEqual("halted_after_nonclean_receipt", verified["artifact_status"])

    def test_exact_inventory_reverifies_prior_stages_and_identity_markers(self):
        plan, inputs = RENDERER.build_plan()
        self.assertEqual(10, len(inputs))
        self.assertEqual(5, len(plan["retained_prior_receipts"]))
        self.assertEqual([15, 16, 17], [row["stage"] for row in plan["retained_prior_expected_summaries"]])
        for task in plan["tasks"]:
            data = inputs[Path(task["named_input"]["path"]).name]
            self.assertIn(task["id"].encode(), data)
            self.assertIn(task["source_instance_sha256"].encode(), data)
            self.assertLessEqual(len(data), 50_000)

    def test_zero_attempt_no_retry_is_not_applicable(self):
        artifact, _, binding = self.materialized()
        summary = VERIFIER.verify(artifact, execution_binding_override=binding)
        self.assertFalse(summary["one_attempt_no_retry_verified"])
        self.assertEqual("not_applicable", summary["one_attempt_no_retry_status"])

    def test_response_swap_and_duplicate_identity_fail(self):
        artifact, plan, _ = self.materialized()
        first, second = plan["tasks"][:2]
        path = artifact / "swapped.xml"
        CHILD.atomic_write(path, self.response_xml(first))
        observed = VERIFIER.classify_response(path, 200, plan["service"]["endpoint"], plan, second)
        self.assertEqual("identity_or_terminal_mismatch", observed["outcome"])
        CHILD.atomic_write(path, self.response_xml(first, duplicate_identity=True))
        observed = VERIFIER.classify_response(path, 200, plan["service"]["endpoint"], plan, first)
        self.assertEqual("identity_or_terminal_mismatch", observed["outcome"])

    def test_final_url_mismatch_is_nonclean(self):
        artifact, plan, _ = self.materialized()
        task = plan["tasks"][0]
        path = artifact / "response.xml"
        CHILD.atomic_write(path, self.response_xml(task))
        observed = VERIFIER.classify_response(path, 200, "https://example.invalid/", plan, task)
        self.assertEqual("final_url_mismatch", observed["outcome"])
        self.assertFalse(observed["clean_terminal"])

    def test_body_reader_is_bounded(self):
        data, exceeded = CHILD.bounded_read(io.BytesIO(b"x" * (CHILD.MAX_RESPONSE_BYTES + 1)))
        self.assertTrue(exceeded)
        self.assertEqual(CHILD.MAX_RESPONSE_BYTES, len(data))

    def test_atomic_write_rejects_symlink_and_hardlink_destination(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "target"
            target.write_bytes(b"old")
            symlink = root / "symlink"
            symlink.symlink_to(target)
            with self.assertRaises(CHILD.ChildError):
                CHILD.atomic_write(symlink, b"new")
            hardlink = root / "hardlink"
            os.link(target, hardlink)
            with self.assertRaises(CHILD.ChildError):
                CHILD.atomic_write(target, b"new")

    def test_initial_execution_requires_committed_manifest_and_clean_checkout(self):
        with self.assertRaises(RUNNER.RunError):
            RUNNER.validate_execution_binding(RUNNER.DEFAULT_ARTIFACT, resume=False)
        with self.assertRaises(PREPARER.PreparationError):
            PREPARER.build_manifest()

    def test_direct_child_arbitrary_attempt_cannot_reach_transport(self):
        with tempfile.TemporaryDirectory() as temporary:
            attempt = Path(temporary) / "attempt"
            attempt.mkdir()
            with mock.patch.object(CHILD, "perform_request") as transport:
                with self.assertRaises(CHILD.ChildError):
                    CHILD.execute(attempt)
            transport.assert_not_called()

    def test_launch_nonce_is_hash_bound_and_not_serialized(self):
        artifact, plan, binding = self.materialized()
        task = plan["tasks"][0]
        attempt, start, nonce = RUNNER.create_attempt(artifact, plan, task, binding)
        expected_hash = RUNNER.sha256_bytes(nonce.encode())
        self.assertEqual(expected_hash, start["launch_nonce_sha256"])
        authorization = RUNNER.read_json(attempt / "launch-authorization.json")
        self.assertEqual(expected_hash, authorization["launch_nonce_sha256"])
        self.assertNotIn(nonce.encode(), (attempt / "attempt-start.json").read_bytes())
        self.assertNotIn(nonce.encode(), (attempt / "launch-authorization.json").read_bytes())

    def test_parent_passes_nonce_only_through_meter_environment(self):
        artifact, plan, binding = self.materialized()
        task = plan["tasks"][0]
        attempt, _, nonce = RUNNER.create_attempt(artifact, plan, task, binding)

        def fake_meter(command, cwd, env, check):  # noqa: ANN001
            self.assertEqual(nonce, env[CHILD.LAUNCH_NONCE_ENV])
            self.assertEqual(str(CHILD.CA_BUNDLE), env[CHILD.CA_BUNDLE_ENV])
            self.assertFalse(any(name in env for name in CHILD.FORBIDDEN_TLS_ENV))
            stdout = Path(command[command.index("--stdout") + 1])
            stderr = Path(command[command.index("--stderr") + 1])
            metrics_path = Path(command[command.index("--metrics") + 1])
            child_command = command[command.index("--") + 1 :]
            CHILD.atomic_write(stdout, b"")
            CHILD.atomic_write(stderr, b"")
            RUNNER.atomic_json(metrics_path, {
                "command": child_command, "returncode": 1, "watchdog_seconds": 75.0,
                "timed_out": False, "orphan_group_terminated": False,
                "metrics": {
                    "wall_seconds": 0.01, "user_seconds": 0.0, "system_seconds": 0.0,
                    "total_core_seconds": 0.0, "single_core_seconds": 0.0,
                    "peak_rss_bytes": 0, "meter": "fresh-process getrusage(RUSAGE_CHILDREN)",
                },
            })
            return SimpleNamespace(returncode=0)

        with mock.patch.object(RUNNER.subprocess, "run", side_effect=fake_meter):
            process = RUNNER.run_metered_transport(attempt, plan, task, nonce)
        self.assertEqual(1, process["metrics"]["returncode"])

    def test_current_child_removes_ps_and_consumes_before_post(self):
        source = (HERE / "post_stage19_magma_calculator_request.py").read_text()
        self.assertNotIn("validate_meter_parent", source)
        self.assertNotIn('["ps"', source)
        self.assertLess(
            source.index("os.replace(authorization_path, consumed_path)"),
            source.index("result = perform_request(input_bytes)"),
        )

    def test_amendment01_authenticates_zero_post_failure(self):
        amendment = AMENDMENT.read_json(AMENDMENT.AMENDMENT)
        original = AMENDMENT.verify_original(amendment)
        self.assertEqual(0, original["post_requests_started"])
        self.assertTrue(original["authorization_present_unconsumed"])
        self.assertEqual(0, original["response_or_envelope_files"])
        self.assertEqual(0.124652, original["charged_transport_process"]["total_core_seconds"])
        self.assertEqual(30244864, original["charged_transport_process"]["peak_rss_bytes"])

    def test_amendment02_preserves_tls_failure_and_prepares_explicit_ca(self):
        summary = AMENDMENT02.summarize()
        self.assertEqual(23, summary["amendment01_tls_failure"]["artifact_files_verified"])
        self.assertEqual("prepared_not_executed", summary["corrected_artifact"]["status"])
        self.assertEqual(0, summary["corrected_artifact"]["attempted_new_requests"])
        tls = summary["explicit_tls_binding"]
        self.assertEqual(CHILD.CA_BUNDLE_SHA256, tls["ca_bundle"]["sha256"])
        self.assertEqual(200, tls["get_probe_http_status"])
        self.assertEqual(TLS_PROBE.EXPECTED_BODY_SHA256, tls["get_probe_body_sha256"])

    def test_tls_environment_overrides_fail_closed(self):
        with mock.patch.dict(os.environ, {CHILD.CA_BUNDLE_ENV: str(CHILD.CA_BUNDLE)}, clear=True):
            CHILD.validate_tls_environment()
            os.environ[CHILD.CA_BUNDLE_ENV] = "/tmp/wrong-ca.pem"
            with self.assertRaises(CHILD.ChildError):
                CHILD.validate_tls_environment()

    def test_static_checks_survive_missing_or_different_live_ca_but_execution_fails(self):
        modules = {
            id(module): module
            for module in (
                CHILD,
                RUNNER.CHILD,
                VERIFIER.CHILD,
                PREPARER.CHILD,
                TLS_PROBE.CHILD,
                AMENDMENT02.CHILD,
                AMENDMENT02.PROBE.CHILD,
                AMENDMENT02.RUNNER.CHILD,
                AMENDMENT02.STAGE19.CHILD,
            )
        }.values()
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            different = temporary_root / "different-ca.pem"
            different.write_bytes(b"not the frozen CA bundle\n")
            for live_path in (temporary_root / "missing-ca.pem", different):
                with self.subTest(live_path=live_path.name), ExitStack() as stack:
                    for module in modules:
                        stack.enter_context(mock.patch.object(module, "CA_BUNDLE", live_path))
                    CHILD.self_test()
                    RUNNER.self_test()
                    VERIFIER.self_test()
                    PREPARER.self_test()
                    TLS_PROBE.self_test()
                    AMENDMENT02.self_test()
                    AMENDMENT02.summarize()
                    with self.assertRaises(CHILD.ChildError):
                        CHILD.ca_bundle_record()
                    with self.assertRaises(PREPARER.CHILD.ChildError):
                        PREPARER.live_external_dependencies()
                    with self.assertRaises(RUNNER.CHILD.ChildError):
                        RUNNER.live_ca_bundle_for_execution(
                            RUNNER.CHILD.expected_ca_bundle_record()
                        )
                    with mock.patch.dict(
                        os.environ,
                        {CHILD.CA_BUNDLE_ENV: str(live_path)},
                        clear=True,
                    ):
                        with self.assertRaises(CHILD.ChildError):
                            CHILD.create_tls_context()
                    with mock.patch.dict(
                        os.environ,
                        {TLS_PROBE.CHILD.CA_BUNDLE_ENV: str(live_path)},
                        clear=True,
                    ):
                        with self.assertRaises(TLS_PROBE.CHILD.ChildError):
                            TLS_PROBE.exact_context()

    def test_manifest_calls_one_fresh_probe_and_rejects_archived_substitution(self):
        fresh = fresh_probe_fixture()
        events = []

        def fake_git(*args):  # noqa: ANN001
            if args == ("status", "--porcelain=v1", "--untracked-files=all"):
                return b""
            if args == ("rev-parse", "HEAD"):
                return ("1" * 40 + "\n").encode()
            if args == ("rev-parse", "HEAD^{tree}"):
                return ("2" * 40 + "\n").encode()
            raise AssertionError(f"unexpected Git call: {args}")

        def fake_blob(_revision, relative):  # noqa: ANN001
            return {
                "path": relative,
                "git_mode": "100644",
                "git_blob_oid": "3" * 40,
                "bytes": 1,
                "sha256": "4" * 64,
            }

        def live_dependencies():
            events.append("live-ca")
            return {"ca_bundle": CHILD.expected_ca_bundle_record()}

        def perform_probe():
            events.append("fresh-probe")
            return fresh

        with (
            mock.patch.object(PREPARER, "git", side_effect=fake_git),
            mock.patch.object(PREPARER, "blob_record", side_effect=fake_blob),
            mock.patch.object(
                PREPARER, "live_external_dependencies", side_effect=live_dependencies
            ),
            mock.patch.object(
                PREPARER.TLS_PROBE, "perform_probe", side_effect=perform_probe
            ) as probe,
        ):
            manifest = PREPARER.build_manifest()
        probe.assert_called_once_with()
        self.assertEqual(["live-ca", "fresh-probe"], events)
        self.assertEqual(fresh, manifest["fresh_preexecution_tls_get_probe"])
        self.assertEqual(
            PREPARER.TLS_PROBE.archived_records(), manifest["archived_tls_get_probe"]
        )

        with (
            mock.patch.object(PREPARER, "git", side_effect=fake_git),
            mock.patch.object(PREPARER, "blob_record", side_effect=fake_blob),
            mock.patch.object(
                PREPARER,
                "live_external_dependencies",
                return_value={"ca_bundle": CHILD.expected_ca_bundle_record()},
            ),
            mock.patch.object(
                PREPARER.TLS_PROBE,
                "perform_probe",
                return_value=PREPARER.TLS_PROBE.verify_archived(),
            ),
        ):
            with self.assertRaises(PREPARER.TLS_PROBE.ProbeError):
                PREPARER.build_manifest()

        with mock.patch.object(
            PREPARER.TLS_PROBE,
            "perform_probe",
            side_effect=AssertionError("self-test attempted network probe"),
        ) as probe:
            PREPARER.self_test()
        probe.assert_not_called()
        with mock.patch.dict(
            os.environ,
            {
                CHILD.CA_BUNDLE_ENV: str(CHILD.CA_BUNDLE),
                CHILD.FORBIDDEN_TLS_ENV[0]: "/tmp/conflicting-ca-dir",
            },
            clear=True,
        ):
            with self.assertRaises(CHILD.ChildError):
                CHILD.validate_tls_environment()

    def test_original_terminal_ledger_cannot_be_selected_for_execution(self):
        original = HERE / "stage-19-magma-calculator-panel-20260909"
        with self.assertRaises(RUNNER.RunError):
            RUNNER.verify_prepared_artifact(original)

    def test_nonclean_halt_is_terminal_and_old_nine_request_p0_is_closed(self):
        artifact, plan, binding = self.materialized()
        with mock.patch.object(
            RUNNER,
            "run_metered_transport",
            side_effect=lambda attempt, current_plan, task, nonce: self.fake_process(
                artifact, attempt, current_plan, task, "transport_error"
            ),
        ) as post:
            state = RUNNER.run(artifact, resume=False, execution_binding_override=binding)
        self.assertEqual(1, post.call_count)
        self.assertEqual("halted_after_nonclean_receipt", state["status"])
        with mock.patch.object(RUNNER, "run_metered_transport") as second_post:
            with self.assertRaises(RUNNER.RunError):
                RUNNER.run(artifact, resume=True, execution_binding_override=binding)
        second_post.assert_not_called()

        state = RUNNER.read_json(artifact / "run.json")
        state["status"] = "running"
        state.pop("finished_at", None)
        state.pop("summary", None)
        RUNNER.atomic_json(artifact / "run.json", state)
        with mock.patch.object(RUNNER, "run_metered_transport") as old_p0_post:
            repaired = RUNNER.run(artifact, resume=True, execution_binding_override=binding)
        old_p0_post.assert_not_called()
        self.assertEqual("halted_after_nonclean_receipt", repaired["status"])

    def test_crash_empty_attempt_directory_is_terminal(self):
        self.assert_recovery_boundary("interrupted_before_attempt_start", lambda *args: None)

    def test_crash_after_start_before_transport_is_terminal(self):
        def populate(artifact, plan, binding, task, attempt):  # noqa: ANN001
            attempt.rmdir()
            RUNNER.create_attempt(artifact, plan, task, binding)
        self.assert_recovery_boundary("interrupted_before_transport_receipt", populate)

    def test_crash_after_response_before_envelope_is_terminal(self):
        def populate(artifact, plan, binding, task, attempt):  # noqa: ANN001
            attempt.rmdir()
            attempt, _, _ = RUNNER.create_attempt(artifact, plan, task, binding)
            CHILD.atomic_write(attempt / "response.xml", self.response_xml(task))
        self.assert_recovery_boundary("interrupted_with_unbound_response", populate)

    def test_crash_after_response_and_envelope_before_receipt_is_terminal(self):
        def populate(artifact, plan, binding, task, attempt):  # noqa: ANN001
            attempt.rmdir()
            attempt, _, _ = RUNNER.create_attempt(artifact, plan, task, binding)
            input_bytes = (artifact / task["named_input"]["path"]).read_bytes()
            body = RUNNER.request_body(input_bytes, "input")
            CHILD.persist_result(attempt, input_bytes, {
                "http_status": 200, "final_url": plan["service"]["endpoint"], "headers": {},
                "body": self.response_xml(task), "body_limit_exceeded": False, "transport_error": None,
                "request_body_bytes": len(body), "request_body_sha256": RUNNER.sha256_bytes(body),
            }, RUNNER.now())
        self.assert_recovery_boundary("interrupted_with_retained_response", populate)

    def test_crash_after_envelope_without_body_is_terminal(self):
        def populate(artifact, plan, binding, task, attempt):  # noqa: ANN001
            attempt.rmdir()
            attempt, _, _ = RUNNER.create_attempt(artifact, plan, task, binding)
            input_bytes = (artifact / task["named_input"]["path"]).read_bytes()
            body = RUNNER.request_body(input_bytes, "input")
            CHILD.persist_result(attempt, input_bytes, {
                "http_status": None, "final_url": None, "headers": {}, "body": None,
                "body_limit_exceeded": False, "transport_error": "synthetic",
                "request_body_bytes": len(body), "request_body_sha256": RUNNER.sha256_bytes(body),
            }, RUNNER.now())
        self.assert_recovery_boundary("interrupted_with_incomplete_transport_envelope", populate)

    def test_crash_after_meter_before_receipt_is_terminal(self):
        def populate(artifact, plan, binding, task, attempt):  # noqa: ANN001
            attempt.rmdir()
            attempt, _, _ = RUNNER.create_attempt(artifact, plan, task, binding)
            RUNNER.atomic_json(attempt / "transport-metrics.json", {"timed_out": True})
        self.assert_recovery_boundary("interrupted_after_transport_process", populate)

    def test_timeout_after_complete_response_binds_all_bytes_nonclaiming(self):
        artifact, plan, binding = self.materialized()
        task = plan["tasks"][0]
        attempt, start, _ = RUNNER.create_attempt(artifact, plan, task, binding)
        process = self.fake_process(
            artifact, attempt, plan, task, "clean", timed_out=True, returncode=-15
        )
        receipt = RUNNER.build_final_receipt(attempt, plan, task, start, process)
        self.assertEqual("transport_timeout", receipt["outcome"])
        self.assertFalse(receipt["claim_admitted"])
        self.assertIsNotNone(receipt["response"])
        self.assertIsNotNone(receipt["transport_envelope"])
        self.assertIn("response.xml", {row["path"] for row in receipt["retained_files"]})
        verified = VERIFIER.verify_attempt(artifact, plan, task, attempt, binding)
        self.assertFalse(verified["clean_f4_terminal"])

    def test_nonzero_after_complete_response_binds_all_bytes_nonclaiming(self):
        artifact, plan, binding = self.materialized()
        task = plan["tasks"][0]
        attempt, start, _ = RUNNER.create_attempt(artifact, plan, task, binding)
        process = self.fake_process(
            artifact, attempt, plan, task, "clean", returncode=2
        )
        receipt = RUNNER.build_final_receipt(attempt, plan, task, start, process)
        self.assertEqual("transport_child_error", receipt["outcome"])
        self.assertFalse(receipt["claim_admitted"])
        self.assertIsNotNone(receipt["response"])
        self.assertIsNotNone(receipt["transport_envelope"])
        VERIFIER.verify_attempt(artifact, plan, task, attempt, binding)

    def test_killed_during_atomic_write_preserves_temp_hash_nonclaiming(self):
        artifact, plan, binding = self.materialized()
        task = plan["tasks"][0]
        attempt, start, _ = RUNNER.create_attempt(artifact, plan, task, binding)
        os.replace(
            attempt / "launch-authorization.json",
            attempt / "launch-authorization.consumed.json",
        )
        CHILD.atomic_write(attempt / "transport.stdout", b"")
        CHILD.atomic_write(attempt / "transport.stderr", b"killed\n")
        command = VERIFIER.expected_transport_command(artifact, task)
        metrics = {
            "command": command,
            "returncode": -9,
            "watchdog_seconds": 75.0,
            "timed_out": True,
            "orphan_group_terminated": False,
            "metrics": {
                "wall_seconds": 75.0, "user_seconds": 0.01, "system_seconds": 0.01,
                "total_core_seconds": 0.02, "single_core_seconds": 0.02,
                "peak_rss_bytes": 1024, "meter": "fresh-process getrusage(RUSAGE_CHILDREN)",
            },
        }
        RUNNER.atomic_json(attempt / "transport-metrics.json", metrics)
        temporary = attempt / ".response.xml.tmp-killed"
        temporary.write_bytes(b"partial")
        process = {
            "metrics": metrics,
            "stdout": RUNNER.regular_record(attempt / "transport.stdout", attempt),
            "stderr": RUNNER.regular_record(attempt / "transport.stderr", attempt),
            "metrics_file": RUNNER.regular_record(attempt / "transport-metrics.json", attempt),
        }
        receipt = RUNNER.build_final_receipt(attempt, plan, task, start, process)
        retained = {row["path"]: row["sha256"] for row in receipt["retained_files"]}
        self.assertEqual(RUNNER.sha256_bytes(b"partial"), retained[temporary.name])
        self.assertFalse(receipt["claim_admitted"])
        VERIFIER.verify_attempt(artifact, plan, task, attempt, binding)

    def test_crash_after_nonclean_receipt_before_run_index_halts_without_post(self):
        artifact, plan, binding = self.materialized()
        self.initialize_run(artifact, plan, binding)
        task = plan["tasks"][0]
        attempt, start, _ = RUNNER.create_attempt(artifact, plan, task, binding)
        process = self.fake_process(artifact, attempt, plan, task, "transport_error")
        receipt = RUNNER.build_final_receipt(attempt, plan, task, start, process)
        self.assertFalse(receipt["clean_terminal"])
        with mock.patch.object(RUNNER, "run_metered_transport") as post:
            state = RUNNER.run(artifact, resume=True, execution_binding_override=binding)
        post.assert_not_called()
        self.assertEqual("halted_after_nonclean_receipt", state["status"])

    def test_crash_after_clean_receipt_before_index_never_retries_first_task(self):
        artifact, plan, binding = self.materialized()
        self.initialize_run(artifact, plan, binding)
        first, second = plan["tasks"][:2]
        past = datetime.now(timezone.utc) - timedelta(seconds=20)
        with mock.patch.object(
            RUNNER,
            "now",
            side_effect=[stamp(past), stamp(past), stamp(past + timedelta(seconds=1))],
        ):
            attempt, start, _ = RUNNER.create_attempt(artifact, plan, first, binding)
            process = self.fake_process(artifact, attempt, plan, first, "clean")
            receipt = RUNNER.build_final_receipt(attempt, plan, first, start, process)
        self.assertTrue(receipt["clean_terminal"])
        seen = []

        def one_later_failure(attempt, current_plan, task, nonce):  # noqa: ANN001
            seen.append(task["id"])
            return self.fake_process(artifact, attempt, current_plan, task, "transport_error")

        with mock.patch.object(RUNNER, "run_metered_transport", side_effect=one_later_failure):
            state = RUNNER.run(artifact, resume=True, execution_binding_override=binding)
        self.assertEqual([second["id"]], seen)
        self.assertEqual("halted_after_nonclean_receipt", state["status"])
        self.assertEqual(2, len(state["receipts"]))

    def test_complete_requires_ten_clean_identity_bound_receipts(self):
        artifact, plan, binding = self.materialized()
        state = self.initialize_run(artifact, plan, binding)
        base = datetime(2026, 9, 10, tzinfo=timezone.utc)
        receipts = []
        for index, task in enumerate(plan["tasks"]):
            start_stamp = stamp(base + timedelta(seconds=3 * index))
            finish_stamp = stamp(base + timedelta(seconds=3 * index + 1))
            with mock.patch.object(
                RUNNER, "now", side_effect=[start_stamp, start_stamp, finish_stamp]
            ):
                attempt, start, _ = RUNNER.create_attempt(artifact, plan, task, binding)
                process = self.fake_process(artifact, attempt, plan, task, "clean")
                receipt = RUNNER.build_final_receipt(attempt, plan, task, start, process)
            receipts.append(receipt)
            state["receipts"][task["id"]] = RUNNER.receipt_index_entry(artifact, attempt, receipt)
        state["status"] = "complete"
        state["finished_at"] = stamp(base + timedelta(seconds=31))
        state["summary"] = RUNNER.summarize_receipts(plan, receipts, "complete")
        RUNNER.atomic_json(artifact / "run.json", state)
        RUNNER.atomic_json(artifact / "summary.json", state["summary"])
        summary = VERIFIER.verify(artifact, execution_binding_override=binding)
        self.assertEqual("ten_clean_requests_complete", summary["artifact_status"])
        self.assertEqual(10, summary["clean_new_f4_terminals"])


if __name__ == "__main__":
    unittest.main()
