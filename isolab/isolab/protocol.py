"""The wire protocol: job specs in, results out.

Both documents are JSON and both are versioned by a ``schema`` field. A spec
is normalised (every default written out, policy defaults applied) before it
is hashed, so two submissions that mean the same thing carry the same
``spec_sha256``, and every result binds the hash of the exact spec it ran.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
import posixpath
import secrets
import time
from importlib import resources
from typing import Any

import jsonschema

SPEC_SCHEMA_ID = "isolab.job/v1"
RESULT_SCHEMA_ID = "isolab.result/v1"

COMPLETED = ("succeeded", "failed")
NOT_COMPLETED = ("timeout", "cancelled", "infra_error")
TERMINAL_STATES = COMPLETED + NOT_COMPLETED + ("dead",)
ACTIVE_STATES = ("queued", "claimed", "running")

MiB = 1024 * 1024

DEFAULT_PERF_EVENTS = [
    "task-clock", "context-switches", "cpu-migrations", "page-faults",
    "cycles", "instructions", "branches", "branch-misses", "cache-misses",
]

#: Policy presets. An explicit value in the spec always wins over these.
POLICY_DEFAULTS: dict[str, dict[str, Any]] = {
    "strict": {
        "settle_s": 5.0, "max_psi_some_avg10": 0.5, "max_other_cpu": 0.02,
        "max_job_cpu_steal_pct": 0.0, "governor": "performance", "turbo": "off",
        "aslr": "any", "min_isolation_tier": "A",
        "require": ["perf_counters", "cpu_partition", "numa_bind", "evicted",
                    "irq_moved", "no_steal", "bare_metal"],
    },
    "standard": {
        "settle_s": 3.0, "max_psi_some_avg10": 1.0, "max_other_cpu": 0.05,
        "max_job_cpu_steal_pct": 0.5, "governor": "any", "turbo": "any",
        "aslr": "any", "min_isolation_tier": "C", "require": [],
    },
    "best_effort": {
        "settle_s": 1.0, "max_psi_some_avg10": 1e9, "max_other_cpu": 1e9,
        "max_job_cpu_steal_pct": 100.0, "governor": "any", "turbo": "any",
        "aslr": "any", "min_isolation_tier": "D", "require": [],
    },
}

TIER_ORDER = {"A": 4, "B": 3, "C": 2, "D": 1, "none": 0}

_DEFAULTS: dict[str, Any] = {
    "runtime": {"backend": "auto", "oci_runtime": "auto", "image": None,
                "network": "none", "user": None, "extra_args": []},
    "measure": {"warmups": 0, "repeats": 1, "cooldown_s": 0.0, "drop_caches": False,
                "sample_period_s": 1.0, "perf_events": DEFAULT_PERF_EVENTS},
    "resources": {"cpus": 1, "memory_mb": None, "numa_node": "single", "smt": "isolate",
                  "gpus": 0, "pids": 4096, "scratch_mb": 0},
    "placement": {"worker": None, "labels": {}, "arch": None, "cpu_flags": [],
                  "min_kernel": None, "require_images": True},
    "limits": {"timeout_s": 3600.0, "build_timeout_s": 3600.0, "max_log_bytes": 16 * MiB,
               "max_artifact_bytes": 1024 * MiB, "max_artifact_files": 10000},
    "retry": {"max_attempts": 2},
    "labels": {},
}
_VERIFY_DEFAULTS = {"certificate_file": "certificate.json", "timeout_s": 600.0}


class SpecError(ValueError):
    """The spec is not a valid isolab.job/v1 document."""


def _load_schema(name: str) -> dict[str, Any]:
    return json.loads(resources.files("isolab.schemas").joinpath(name).read_text())


SPEC_SCHEMA = _load_schema("job-spec.v1.json")
RESULT_SCHEMA = _load_schema("result.v1.json")


def _safe_relpath(path: str, what: str) -> str:
    norm = posixpath.normpath(path)
    if posixpath.isabs(norm) or norm == ".." or norm.startswith("../") or norm == ".":
        if norm != "." or what != "cwd":
            raise SpecError(f"{what}: {path!r} must be a relative path inside the working directory")
    return norm


def normalize_spec(spec: dict[str, Any], allow_local_files: bool = False) -> dict[str, Any]:
    """Validate a spec and return a copy with every default made explicit.

    ``local_file`` inputs are only legal before submission; the client turns
    them into ``sha256`` blobs (see :func:`isolab.blobs.resolve_local_inputs`).
    """
    try:
        jsonschema.validate(spec, SPEC_SCHEMA)
    except jsonschema.ValidationError as err:
        where = "/".join(str(p) for p in err.absolute_path) or "<root>"
        raise SpecError(f"{where}: {err.message}") from None
    out = copy.deepcopy(spec)
    out.setdefault("pool", "default")
    out.setdefault("name", None)
    out.setdefault("inputs", [])
    out.setdefault("build", [])
    for key, default in _DEFAULTS.items():
        merged = copy.deepcopy(default)
        merged.update(out.get(key) or {})
        out[key] = merged
    fid = dict(out.get("fidelity") or {})
    policy = fid.get("policy", "standard")
    merged_fid = {"policy": policy, **copy.deepcopy(POLICY_DEFAULTS[policy])}
    for k, v in fid.items():
        if k == "require":
            merged_fid["require"] = sorted(set(merged_fid["require"]) | set(v))
        else:
            merged_fid[k] = v
    out["fidelity"] = merged_fid
    if out.get("verify"):
        out["verify"] = {**_VERIFY_DEFAULTS, **out["verify"]}
    else:
        out["verify"] = None
    out.setdefault("idempotency_key", None)

    cmd = out["command"]
    cmd.setdefault("cwd", ".")
    cmd.setdefault("env", {})
    cmd.setdefault("stdin", None)
    cmd["cwd"] = _safe_relpath(cmd["cwd"], "cwd")
    for i, step in enumerate(out["build"]):
        step.setdefault("cwd", ".")
        step.setdefault("env", {})
        step.setdefault("timeout_s", None)
        step["cwd"] = _safe_relpath(step["cwd"], "cwd")

    seen: set[str] = set()
    for i, inp in enumerate(out["inputs"]):
        inp["path"] = _safe_relpath(inp["path"], f"inputs/{i}/path")
        if inp["path"] in seen:
            raise SpecError(f"inputs/{i}/path: {inp['path']!r} given twice")
        seen.add(inp["path"])
        kinds = [k for k in ("content", "local_file", "sha256", "git") if k in inp]
        if len(kinds) != 1:
            raise SpecError(f"inputs/{i}: give exactly one of content, local_file, sha256, git "
                            f"(got {kinds or 'none'})")
        kind = kinds[0]
        if kind == "local_file" and not allow_local_files:
            raise SpecError(f"inputs/{i}: local_file must be resolved to a sha256 blob before submission")
        if kind == "content":
            inp.setdefault("encoding", "utf8")
        if kind == "sha256" and "bytes" not in inp:
            raise SpecError(f"inputs/{i}: a sha256 input needs its byte count")
        if kind == "git":
            for sp in inp["git"].get("sparse_paths") or []:
                _safe_relpath(sp, f"inputs/{i}/git/sparse_paths")
            if not inp["git"].get("sparse_paths"):
                inp["git"]["sparse_paths"] = []
        inp.setdefault("mode", "0644" if kind != "git" else None)
        if inp.get("mode") and not inp["mode"].startswith("0"):
            inp["mode"] = "0" + inp["mode"]
    if out["runtime"]["backend"] == "direct" and out["runtime"]["image"]:
        raise SpecError("runtime/image: the direct backend takes no image")
    if out["runtime"]["oci_runtime"] == "runsc" and out["runtime"]["backend"] == "direct":
        raise SpecError("runtime/oci_runtime: runsc needs a container backend")
    res = out["resources"]
    if res["memory_mb"] is not None and res["memory_mb"] < 16:
        raise SpecError("resources/memory_mb: at least 16")
    return out


def canonical_json(doc: Any) -> bytes:
    return json.dumps(doc, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def spec_sha256(normalized: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_json(normalized)).hexdigest()


def new_job_id() -> str:
    """``J-<ms timestamp base36>-<6 hex>``: sortable by creation, never reused."""
    ms = int(time.time() * 1000)
    digits = "0123456789abcdefghijklmnopqrstuvwxyz"
    b36 = ""
    while ms:
        ms, r = divmod(ms, 36)
        b36 = digits[r] + b36
    return f"J-{b36}-{secrets.token_hex(3)}"


def validate_result(result: dict[str, Any]) -> None:
    jsonschema.validate(result, RESULT_SCHEMA)


def outcome_class(status: str) -> str:
    return "completed" if status in COMPLETED else "not_completed"


def tier_at_least(have: str, want: str | None) -> bool:
    return want is None or TIER_ORDER.get(have, 0) >= TIER_ORDER.get(want, 0)


def summarize_spec(spec: dict[str, Any]) -> dict[str, Any]:
    """A one-line view of a normalised spec for listings."""
    return {
        "name": spec.get("name"), "pool": spec["pool"],
        "image": spec["runtime"]["image"], "backend": spec["runtime"]["backend"],
        "argv": spec["command"]["argv"], "cpus": spec["resources"]["cpus"],
        "memory_mb": spec["resources"]["memory_mb"], "gpus": spec["resources"]["gpus"],
        "repeats": spec["measure"]["repeats"], "policy": spec["fidelity"]["policy"],
        "worker": spec["placement"]["worker"],
    }


def env_contract(job_id: str, attempt: int, repeat: int | None, warmup: bool | None,
                 output_dir: str, work_dir: str, scratch_dir: str) -> dict[str, str]:
    env = {"ISOLAB_JOB_ID": job_id, "ISOLAB_ATTEMPT": str(attempt),
           "ISOLAB_OUTPUT_DIR": output_dir, "ISOLAB_WORK": work_dir,
           "ISOLAB_SCRATCH": scratch_dir}
    if repeat is not None:
        env["ISOLAB_REPEAT"] = str(repeat)
        env["ISOLAB_WARMUP"] = "1" if warmup else "0"
    return env


def host_env_passthrough() -> dict[str, str]:
    keep = ("PATH", "HOME", "LANG", "LC_ALL", "TZ", "USER", "CARGO_HOME", "RUSTUP_HOME",
            "SAGE_ROOT", "PYTHONPATH", "VIRTUAL_ENV", "CUDA_HOME", "LD_LIBRARY_PATH")
    return {k: os.environ[k] for k in keep if k in os.environ}
