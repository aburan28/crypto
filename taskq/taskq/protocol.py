"""The wire protocol: task specs in, task results out.

Both documents are JSON and both are versioned by a `schema` field. A spec is
normalised (defaults filled in) before it is hashed, so two submissions that
mean the same thing carry the same `spec_sha256`, and every result binds the
hash of the exact spec it executed.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
import secrets
import time
from importlib import resources
from pathlib import Path
from typing import Any

import jsonschema

SPEC_SCHEMA_ID = "taskq.task-spec/v1"
RESULT_SCHEMA_ID = "taskq.task-result/v1"

# Terminal statuses, and which of them mean "the command ran to exit".
COMPLETED = ("succeeded", "failed")
NOT_COMPLETED = ("timeout", "cancelled", "infra_error")
TERMINAL_STATES = COMPLETED + NOT_COMPLETED + ("dead",)


class SpecError(ValueError):
    """The spec is not a valid taskq.task-spec/v1 document."""


def _load_schema(name: str) -> dict[str, Any]:
    text = resources.files("taskq.schemas").joinpath(name).read_text()
    return json.loads(text)


SPEC_SCHEMA = _load_schema("task-spec.v1.json")
RESULT_SCHEMA = _load_schema("task-result.v1.json")

_DEFAULTS: dict[str, Any] = {
    "benchmark": {"warmups": 1, "repetitions": 5},
    "limits": {"timeout_seconds": 3600, "setup_timeout_seconds": 3600,
               "memory_mb": None, "max_log_bytes": 16 * 1024 * 1024},
    "placement": {"require_labels": {}, "cpus": None},
    "retry": {"max_attempts": 3},
    "labels": {},
}
_VERIFY_DEFAULTS = {"timeout_seconds": 600, "certificate_file": "certificate.json"}


def normalize_spec(spec: dict[str, Any]) -> dict[str, Any]:
    """Validate a spec and return a copy with every default made explicit."""
    try:
        jsonschema.validate(spec, SPEC_SCHEMA)
    except jsonschema.ValidationError as err:
        where = "/".join(str(p) for p in err.absolute_path) or "<root>"
        raise SpecError(f"{where}: {err.message}") from None
    out = copy.deepcopy(spec)
    for key, default in _DEFAULTS.items():
        if key == "benchmark" and out["kind"] != "benchmark":
            out.pop("benchmark", None)
            continue
        merged = copy.deepcopy(default)
        merged.update(out.get(key) or {})
        out[key] = merged
    if out.get("verify"):
        out["verify"] = {**_VERIFY_DEFAULTS, **out["verify"]}
    cmd = out["command"]
    cmd.setdefault("cwd", ".")
    cmd.setdefault("env", {})
    cmd.setdefault("setup", [])
    cwd = os.path.normpath(cmd["cwd"])
    if os.path.isabs(cwd) or cwd == ".." or cwd.startswith("../"):
        raise SpecError(f"command/cwd: {cmd['cwd']!r} escapes the checkout")
    for sp in out["source"].get("sparse_paths") or []:
        if ".." in Path(sp).parts:
            raise SpecError(f"source/sparse_paths: {sp!r} escapes the checkout")
    return out


def canonical_json(doc: Any) -> bytes:
    return json.dumps(doc, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def spec_sha256(normalized: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_json(normalized)).hexdigest()


def new_task_id() -> str:
    """`T-<ms timestamp base36>-<8 hex>`: sortable by creation, never reused."""
    ms = int(time.time() * 1000)
    digits = "0123456789abcdefghijklmnopqrstuvwxyz"
    b36 = ""
    while ms:
        ms, r = divmod(ms, 36)
        b36 = digits[r] + b36
    return f"T-{b36}-{secrets.token_hex(4)}"


def validate_result(result: dict[str, Any]) -> None:
    jsonschema.validate(result, RESULT_SCHEMA)


def outcome_class(status: str) -> str:
    return "completed" if status in COMPLETED else "not_completed"
