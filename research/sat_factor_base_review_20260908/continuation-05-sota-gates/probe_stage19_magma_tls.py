#!/usr/bin/env python3
"""Verify or freshly repeat the GET-only Stage 19 TLS endpoint probe."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import socket
import ssl
import time
from typing import Any
import urllib.error
import urllib.request


HERE = Path(__file__).resolve().parent
ARTIFACT = HERE / "stage-19-tls-get-probe-20260910"
CHILD_PATH = HERE / "post_stage19_magma_calculator_request.py"
EXPECTED_BODY = b'<?xml version="1.0"?>\n<calculator/>\n'
EXPECTED_BODY_SHA256 = "c96bf54dfc5469d7a6178b115e82337ea5f1552038028cf2fa14a38db1e72161"
FRESH_SCHEMA = "koblitz_magma_calculator_fresh_tls_get_probe.v1"


class ProbeError(RuntimeError):
    """The exact-context GET-only TLS probe failed closed."""


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ProbeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CHILD = load_module("stage19_child_for_tls_probe", CHILD_PATH)


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise ProbeError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise ProbeError(f"expected JSON object in {path}")
    return value


def validate_probe_environment() -> None:
    selected = os.environ.get(CHILD.CA_BUNDLE_ENV)
    if selected not in {None, "", str(CHILD.CA_BUNDLE)}:
        raise ProbeError("SSL_CERT_FILE conflicts with the selected Stage 19 CA bundle")
    mismatches = {
        name: os.environ[name]
        for name in CHILD.FORBIDDEN_TLS_ENV
        if os.environ.get(name)
    }
    if mismatches:
        raise ProbeError(f"conflicting TLS environment overrides: {sorted(mismatches)}")


def exact_context() -> ssl.SSLContext:
    validate_probe_environment()
    CHILD.ca_bundle_record()
    return ssl.create_default_context(cafile=str(CHILD.CA_BUNDLE))


def archived_records() -> dict:
    receipt = verify_archived()
    return {
        "receipt": CHILD.regular_record(ARTIFACT / "receipt.json", HERE),
        "response": CHILD.regular_record(ARTIFACT / "response.xml", HERE),
        "http_status": receipt["response"]["http_status"],
        "no_computation_requested": receipt["request"]["no_computation_requested"],
    }


def verify_archived() -> dict:
    receipt = read_json(ARTIFACT / "receipt.json")
    body = (ARTIFACT / "response.xml").read_bytes()
    if len(body) != 36 or body != EXPECTED_BODY or sha256_bytes(body) != EXPECTED_BODY_SHA256:
        raise ProbeError("archived GET probe response bytes changed")
    expected = {
        "schema": "koblitz_magma_calculator_tls_get_probe.v1",
        "observed_at_date": "2026-09-10",
        "request": {
            "method": "GET",
            "endpoint": CHILD.ENDPOINT,
            "no_calculator_input": True,
            "no_computation_requested": True,
            "redirects_allowed": False,
            "socket_timeout_seconds": CHILD.SOCKET_TIMEOUT_SECONDS,
            "response_byte_limit": CHILD.MAX_RESPONSE_BYTES,
            "ca_bundle": CHILD.expected_ca_bundle_record(),
            "tls_context": "ssl.create_default_context(cafile='/etc/ssl/cert.pem')",
        },
        "response": {
            "http_status": 200,
            "final_url": CHILD.ENDPOINT,
            "body_path": "response.xml",
            "body_bytes": 36,
            "body_sha256": EXPECTED_BODY_SHA256,
            "body_complete": True,
            "date_header_retained": False,
        },
        "classification": "tls_endpoint_reachability_only",
        "claim_boundary": (
            "A GET-only TLS reachability probe using the selected CA bundle. It submitted no "
            "calculator input, requested no Magma computation, and is not solver, scientific, "
            "performance, or POST-delivery evidence."
        ),
    }
    if receipt != expected:
        raise ProbeError("archived GET probe receipt changed")
    return receipt


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):  # noqa: ANN001
        return None


def perform_probe() -> dict:
    archived = verify_archived()
    request = urllib.request.Request(
        CHILD.ENDPOINT,
        method="GET",
        headers={"User-Agent": CHILD.USER_AGENT},
    )
    opener = urllib.request.build_opener(
        NoRedirect(), urllib.request.HTTPSHandler(context=exact_context())
    )
    started_at = now()
    started = time.monotonic()
    try:
        with opener.open(request, timeout=CHILD.SOCKET_TIMEOUT_SECONDS) as response:
            body = response.read(CHILD.MAX_RESPONSE_BYTES + 1)
            status = int(response.status)
            final_url = response.geturl()
            headers = {
                name: response.headers.get(name)
                for name in ("content-type", "content-length", "server")
                if response.headers.get(name) is not None
            }
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, socket.timeout, OSError) as error:
        raise ProbeError(f"fresh pre-execution TLS GET failed: {type(error).__name__}: {error}") from error
    wall = time.monotonic() - started
    if (
        status != 200
        or final_url != CHILD.ENDPOINT
        or body != EXPECTED_BODY
        or len(body) != 36
        or sha256_bytes(body) != EXPECTED_BODY_SHA256
    ):
        raise ProbeError("fresh TLS GET response differs from the frozen reachability contract")
    result = {
        "schema": FRESH_SCHEMA,
        "started_at": started_at,
        "finished_at": now(),
        "wall_seconds": wall,
        "request": archived["request"],
        "response": {
            **archived["response"],
            "selected_headers_without_date": headers,
        },
        "fresh_immediately_before_manifest_generation": True,
        "classification": archived["classification"],
        "claim_boundary": archived["claim_boundary"],
    }
    result["receipt_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def self_test() -> dict:
    receipt = verify_archived()
    return {
        "self_test": "pass",
        "method": receipt["request"]["method"],
        "body_bytes": receipt["response"]["body_bytes"],
        "body_sha256": receipt["response"]["body_sha256"],
        "ca_bundle_sha256": receipt["request"]["ca_bundle"]["sha256"],
        "no_computation_requested": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    if not args.probe:
        parser.error("a fresh network GET requires --probe")
    try:
        result = perform_probe()
    except ProbeError as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
