#!/usr/bin/env python3
"""Verify Stage 19 Amendment 02, its TLS-failure custody, and explicit CA binding."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
AMENDMENT = HERE / "stage-19-amendment-02-explicit-ca.json"
AMENDMENT01 = HERE / "stage-19-magma-calculator-panel-amendment-01-20260910"
CORRECTED = HERE / "stage-19-magma-calculator-panel-amendment-02-20260910"
SUMMARY = HERE / "stage-19-amendment-02-summary.json"
SCHEMA = "koblitz_magma_calculator_stage19_amendment_verification.v2"


class VerificationError(RuntimeError):
    """Amendment 02 or one of its immutable inputs is inconsistent."""


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise VerificationError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CHILD = load_module("stage19_child_for_amendment02", HERE / "post_stage19_magma_calculator_request.py")
RUNNER = load_module("stage19_runner_for_amendment02", HERE / "run_stage19_magma_calculator_panel.py")
STAGE19 = load_module("stage19_verifier_for_amendment02", HERE / "verify_stage19_magma_calculator_panel.py")
PROBE = load_module("stage19_probe_for_amendment02", HERE / "probe_stage19_magma_tls.py")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise VerificationError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise VerificationError(f"expected JSON object in {path}")
    return value


def file_record(path: Path, relative_to: Path) -> dict:
    try:
        info = path.lstat()
    except FileNotFoundError as error:
        raise VerificationError(f"missing immutable file {path}") from error
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise VerificationError(f"not a single-link regular file: {path}")
    data = path.read_bytes()
    return {
        "path": str(path.relative_to(relative_to)),
        "bytes": len(data),
        "sha256": sha256_bytes(data),
    }


def artifact_inventory(root: Path) -> list[dict]:
    if root.is_symlink() or not root.is_dir():
        raise VerificationError(f"artifact is absent or a symlink: {root}")
    files = [path for path in sorted(root.rglob("*")) if path.is_file() and not path.is_symlink()]
    if any(path.is_symlink() for path in root.rglob("*")):
        raise VerificationError(f"artifact contains a symlink: {root}")
    return [file_record(path, root) for path in files]


def verify_amendment01_failure(amendment: dict) -> dict:
    frozen = amendment["amendment01_tls_failure"]
    inventory = artifact_inventory(AMENDMENT01)
    if (
        len(inventory) != frozen["artifact_file_count"]
        or sha256_bytes(canonical_bytes(inventory)) != frozen["artifact_inventory_sha256"]
    ):
        raise VerificationError("Amendment 01 TLS-failure artifact bytes changed")
    attempt = AMENDMENT01 / "attempts" / "01-seed-2026091301-n31-l5-m3-ggmp-a0-f0"
    envelope = read_json(attempt / "transport-envelope.json")
    receipt = read_json(attempt / "receipt.json")
    run = read_json(AMENDMENT01 / "run.json")
    summary = read_json(AMENDMENT01 / "summary.json")
    if (
        sha256_bytes((AMENDMENT01 / "run.json").read_bytes()) != frozen["run_sha256"]
        or sha256_bytes((AMENDMENT01 / "summary.json").read_bytes()) != frozen["summary_sha256"]
        or sha256_bytes((attempt / "transport-envelope.json").read_bytes())
        != frozen["transport_envelope_sha256"]
        or sha256_bytes((attempt / "request-body.bin").read_bytes())
        != frozen["request_body_sha256"]
    ):
        raise VerificationError("Amendment 01 indexed TLS-failure bytes changed")
    if (
        run.get("status") != frozen["terminal_status"]
        or run.get("summary") != summary
        or summary.get("attempted_tasks") != 1
        or summary.get("clean_f4_terminals") != 0
        or envelope.get("http", {}).get("transport_error") != frozen["transport_error"]
        or envelope.get("http", {}).get("status") is not None
        or envelope.get("response") is not None
        or receipt.get("outcome") != frozen["outcome"]
        or receipt.get("claim_admitted") is not False
        or receipt.get("clean_terminal") is not False
        or receipt.get("response") is not None
        or not (attempt / "launch-authorization.consumed.json").is_file()
        or (attempt / "launch-authorization.json").exists()
    ):
        raise VerificationError("Amendment 01 TLS-failure semantics changed")
    return {
        "artifact_files_verified": len(inventory),
        "artifact_inventory_sha256": sha256_bytes(canonical_bytes(inventory)),
        "attempts": 1,
        "post_attempts_started": 1,
        "http_responses_received": 0,
        "calculator_results_received": 0,
        "claim_admitted": False,
        "terminal_status": run["status"],
        "outcome": receipt["outcome"],
    }


def verify_tls_implementation(amendment: dict) -> dict:
    ca = CHILD.expected_ca_bundle_record()
    if ca != amendment["explicit_ca"]:
        raise VerificationError("frozen CA bundle record differs from Amendment 02")
    probe = PROBE.verify_archived()
    recorded_probe = amendment["tls_get_probe"]
    if (
        probe["response"]["http_status"] != recorded_probe["http_status"]
        or probe["response"]["body_bytes"] != recorded_probe["body_bytes"]
        or probe["response"]["body_sha256"] != recorded_probe["body_sha256"]
        or probe["request"]["method"] != recorded_probe["method"]
        or probe["request"]["no_computation_requested"] is not True
        or sha256_bytes((PROBE.ARTIFACT / "receipt.json").read_bytes())
        != recorded_probe["receipt_sha256"]
    ):
        raise VerificationError("archived no-compute TLS GET probe changed")
    child_source = (HERE / "post_stage19_magma_calculator_request.py").read_text()
    preparer_source = (HERE / "prepare_stage19_magma_execution_manifest.py").read_text()
    runner_source = (HERE / "run_stage19_magma_calculator_panel.py").read_text()
    probe_source = (HERE / "probe_stage19_magma_tls.py").read_text()
    required = {
        "child": (
            "ssl.create_default_context(cafile=str(CA_BUNDLE))",
            "manifest.get(\"external_dependencies\", {}).get(\"ca_bundle\") != ca_bundle_record()",
            "SSL_CERT_FILE does not name the bound CA bundle",
            "conflicting TLS environment overrides are set",
        ),
        "preparer": (
            'return {"ca_bundle": CHILD.ca_bundle_record()}',
            '"external_dependencies": live_external_dependencies()',
        ),
        "runner": (
            "observed = CHILD.ca_bundle_record()",
            "ca_bundle = live_ca_bundle_for_execution(",
        ),
        "probe": (
            "CHILD.ca_bundle_record()",
            "ssl.create_default_context(cafile=str(CHILD.CA_BUNDLE))",
        ),
    }
    sources = {
        "child": child_source,
        "preparer": preparer_source,
        "runner": runner_source,
        "probe": probe_source,
    }
    for label, markers in required.items():
        if any(marker not in sources[label] for marker in markers):
            raise VerificationError(f"{label} source lacks a live explicit-CA check")
    if not (
        child_source.index('manifest.get("external_dependencies", {}).get("ca_bundle") != ca_bundle_record()')
        < child_source.index("os.replace(authorization_path, consumed_path)")
        < child_source.index("result = perform_request(input_bytes)")
    ):
        raise VerificationError("live CA authentication is not before authorization consume and POST")
    old = {name: os.environ.get(name) for name in (CHILD.CA_BUNDLE_ENV, *CHILD.FORBIDDEN_TLS_ENV)}
    try:
        os.environ[CHILD.CA_BUNDLE_ENV] = "/tmp/not-the-bound-ca.pem"
        try:
            CHILD.validate_tls_environment()
        except CHILD.ChildError:
            pass
        else:
            raise VerificationError("mismatched SSL_CERT_FILE was accepted")
        os.environ[CHILD.CA_BUNDLE_ENV] = str(CHILD.CA_BUNDLE)
        os.environ[CHILD.FORBIDDEN_TLS_ENV[0]] = "/tmp/conflicting-ca-dir"
        try:
            CHILD.validate_tls_environment()
        except CHILD.ChildError:
            pass
        else:
            raise VerificationError("conflicting TLS environment override was accepted")
    finally:
        for name, value in old.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
    return {
        "ca_bundle": ca,
        "context": amendment["tls_policy"]["context"],
        "environment_mismatch_rejected": True,
        "get_probe_http_status": probe["response"]["http_status"],
        "get_probe_body_bytes": probe["response"]["body_bytes"],
        "get_probe_body_sha256": probe["response"]["body_sha256"],
        "get_probe_no_computation_requested": True,
    }


def summarize() -> dict:
    amendment = read_json(AMENDMENT)
    if (
        amendment.get("schema") != "koblitz_magma_calculator_stage19_amendment.v1"
        or amendment.get("amendment") != 2
    ):
        raise VerificationError("unexpected Stage 19 Amendment 02 record")
    prior = verify_amendment01_failure(amendment)
    tls = verify_tls_implementation(amendment)
    corrected = STAGE19.verify(CORRECTED)
    if corrected.get("artifact_status") != "prepared_not_executed":
        raise VerificationError("Amendment 02 corrected artifact is not prepared-only")
    expected_bound = RUNNER.expected_execution_bound_paths()
    if not all(
        any(fragment in path for path in expected_bound)
        for fragment in (
            "stage-19-magma-calculator-panel-20260909/",
            "stage-19-magma-calculator-panel-amendment-01-20260910/",
            "stage-19-tls-get-probe-20260910/",
        )
    ):
        raise VerificationError("next execution manifest omits a prior Stage 19 artifact")
    return {
        "schema": SCHEMA,
        "amendment_sha256": sha256_bytes(AMENDMENT.read_bytes()),
        "amendment01_tls_failure": prior,
        "explicit_tls_binding": tls,
        "corrected_artifact": {
            "path": str(CORRECTED.relative_to(REPO)),
            "plan_sha256": corrected["plan_sha256"],
            "status": corrected["artifact_status"],
            "attempted_new_requests": corrected["attempted_new_requests"],
        },
        "next_manifest_bound_path_count": len(expected_bound),
        "disposition": "verified_tls_failure_with_explicit_ca_prepared_successor",
        "claim_boundary": amendment["claim_boundary"],
    }


def self_test() -> dict:
    amendment = read_json(AMENDMENT)
    prior = verify_amendment01_failure(amendment)
    tls = verify_tls_implementation(amendment)
    return {
        "self_test": "pass",
        "amendment": 2,
        "amendment01_artifact_files_verified": prior["artifact_files_verified"],
        "ca_bundle_sha256": tls["ca_bundle"]["sha256"],
        "get_probe_body_sha256": tls["get_probe_body_sha256"],
        "no_post_performed": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    try:
        summary = summarize()
    except (
        VerificationError,
        CHILD.ChildError,
        RUNNER.RunError,
        STAGE19.VerificationError,
        PROBE.ProbeError,
    ) as error:
        parser.error(str(error))
    rendered = canonical_bytes(summary)
    if args.expected is not None and args.expected.read_bytes() != rendered:
        parser.error("recomputed Amendment 02 summary differs from expected")
    if args.output is not None:
        CHILD.atomic_write(args.output.resolve(), rendered)
    print(rendered.decode(), end="")


if __name__ == "__main__":
    main()
