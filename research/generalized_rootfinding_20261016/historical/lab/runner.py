"""Bounded, persistent experiment orchestration; scientific release is external.

No shell commands, background scheduling, automatic retries, or implicit release
claims. Every attempt freezes identities before starting one cold process.
"""
from __future__ import annotations

import argparse
import contextlib
import copy
import datetime as dt
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile
import time


SCHEMA = "rootfinding_research_program_v1"
RESULT_SCHEMA = "scientific_result_v1"
OUTCOMES = {"passed", "refuted", "inconclusive", "blocked", "error"}
ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}\Z")


class LabError(ValueError):
    pass


class LabBusy(LabError):
    pass


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(value).hexdigest()


def strict_json(value):
    def pairs(items):
        result = {}
        for key, item in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = item
        return result
    def constant(value):
        raise ValueError(f"nonfinite JSON number: {value}")
    return json.loads(value, object_pairs_hook=pairs, parse_constant=constant)


def run_cold_process(command, *, timeout, workdir, memory_mib, cpu, cancelled=lambda: False):
    """One cancellable cold process, with the frozen harness's resource limits.

    wait4 measures this child, not cumulative children. Both interruption and
    timeout terminate and reap its process group. No benchmarking is parallel.
    """
    from performance_gate import _limits
    started = time.monotonic()
    started_unix = time.time()
    with tempfile.TemporaryFile() as stdout, tempfile.TemporaryFile() as stderr:
        try:
            child = subprocess.Popen(command, cwd=workdir, stdout=stdout, stderr=stderr,
                                     start_new_session=True, preexec_fn=_limits(memory_mib, cpu))
        except (OSError, subprocess.SubprocessError) as error:
            return {"status": "launch_error", "command": command, "reason": str(error),
                    "exit_code": None, "stdout": "", "stderr": "", "wall_seconds": time.monotonic() - started}
        outcome = None
        reaped = False
        try:
            while True:
                pid, status, usage = os.wait4(child.pid, os.WNOHANG)
                if pid:
                    reaped = True
                    break
                if cancelled() or time.monotonic() - started >= timeout:
                    outcome = "interrupted" if cancelled() else "timeout"
                    try:
                        os.killpg(child.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    deadline = time.monotonic() + .20
                    while time.monotonic() < deadline:
                        pid, status, usage = os.wait4(child.pid, os.WNOHANG)
                        if pid:
                            reaped = True
                            break
                        time.sleep(.005)
                    if not reaped:
                        try:
                            os.killpg(child.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        pid, status, usage = os.wait4(child.pid, 0)
                        reaped = True
                    break
                time.sleep(.005)
            child.returncode = os.waitstatus_to_exitcode(status)
        finally:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            if not reaped:
                try:
                    _, status, usage = os.wait4(child.pid, 0)
                    child.returncode = os.waitstatus_to_exitcode(status)
                except ChildProcessError:
                    pass
        stdout.seek(0)
        stderr.seek(0)
        out = stdout.read().decode("utf-8", errors="replace")
        err = stderr.read().decode("utf-8", errors="replace")
    return {"status": outcome or ("ok" if child.returncode == 0 else "process_error"),
            "command": command, "exit_code": child.returncode, "stdout": out, "stderr": err,
            "pid": child.pid, "process_group": child.pid, "wall_seconds": time.monotonic() - started,
            "started_unix": started_unix, "ended_unix": time.time(),
            "stdout_sha256": digest(out.encode()), "user_seconds": usage.ru_utime,
            "system_seconds": usage.ru_stime, "peak_process_rss_kib": usage.ru_maxrss,
            "measurement_scope": "development diagnostic; isolation must be established separately",
            "memory_limit_scope": "per-process address space inherited by descendants"}


def timestamp():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def positive(value, name, *, integer=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise LabError(f"{name} must be finite and positive")
    if integer and type(value) is not int:
        raise LabError(f"{name} must be an integer")
    return value


def contained(root, value):
    path = (root / value).resolve()
    if not path.is_relative_to(root):
        raise LabError(f"path escapes project root: {value}")
    return path


def normalize_manifest(path):
    path = Path(path).resolve()
    manifest = strict_json(path.read_bytes())
    if manifest.get("schema") != SCHEMA or not ID.fullmatch(str(manifest.get("program_id", ""))):
        raise LabError("invalid research program schema or id")
    root = (path.parent / manifest.get("project_root", "..")).resolve()
    if not root.is_dir():
        raise LabError("project_root does not exist")
    ledger = contained(root, manifest.get("ledger_dir", "lab/ledger"))
    items = manifest.get("experiments")
    if not isinstance(items, list):
        raise LabError("experiments must be a list")
    seen = set()
    for item in items:
        eid = item.get("id", "")
        if not isinstance(eid, str) or not ID.fullmatch(eid) or eid in seen:
            raise LabError("experiment ids must be unique safe identifiers")
        seen.add(eid)
        if not isinstance(item.get("hypothesis"), str) or not item["hypothesis"].strip():
            raise LabError(f"{eid}: hypothesis required")
        command = item.get("command")
        if (not isinstance(command, list) or len(command) < 3 or
                any(not isinstance(arg, str) or "\x00" in arg for arg in command) or
                command[1] != "-m" or
                not re.fullmatch(r"(?:lab|research)\.[A-Za-z_][A-Za-z_0-9.]*", command[2]) or
                (command[0] not in ("python", "python3", sys.executable) and
                 Path(command[0]).resolve() != Path(sys.executable).resolve())):
            raise LabError(f"{eid}: command must invoke a project Python module as an argument array")
        for name in ("source_paths", "input_paths", "artifact_paths", "depends_on"):
            values = item.setdefault(name, [])
            if not isinstance(values, list) or any(not isinstance(v, str) for v in values) or len(set(values)) != len(values):
                raise LabError(f"{eid}: {name} must be unique strings")
        if not item["source_paths"]:
            raise LabError(f"{eid}: explicit source_paths required")
        worker_path = root.joinpath(*command[2].split(".")).with_suffix(".py")
        if not worker_path.is_file() or worker_path not in {contained(root, value) for value in item["source_paths"]}:
            raise LabError(f"{eid}: source_paths must include the invoked worker module")
        for value in item["source_paths"] + item["input_paths"]:
            if not contained(root, value).is_file():
                raise LabError(f"{eid}: source/input file missing: {value}")
        for value in item["artifact_paths"]:
            artifact = contained(root, value)
            if artifact == ledger or artifact.is_relative_to(ledger):
                raise LabError(f"{eid}: worker artifact paths cannot point inside the ledger")
        positive(item.setdefault("timeout_seconds", 30), f"{eid} timeout")
        positive(item.setdefault("memory_limit_mb", 4096), f"{eid} memory", integer=True)
        cpu = item.setdefault("cpu", None)
        if cpu is not None and (type(cpu) is not int or cpu < 0):
            raise LabError(f"{eid}: cpu must be a nonnegative integer or null")
        if "retry_of" in item and (not isinstance(item.get("retry_reason"), str) or not item["retry_reason"].strip()):
            raise LabError(f"{eid}: an explicit retry needs retry_reason")
    for item in items:
        if item["id"] in item["depends_on"] or any(v not in seen for v in item["depends_on"]):
            raise LabError(f"{item['id']}: unknown or self dependency")
        if item.get("retry_of") is not None and item["retry_of"] not in seen:
            raise LabError(f"{item['id']}: unknown retry_of")
    def visit(eid, ancestry):
        if eid in ancestry:
            raise LabError("cyclic dependencies")
        entry = next(v for v in items if v["id"] == eid)
        for dep in entry["depends_on"]:
            visit(dep, ancestry | {eid})
    for eid in seen:
        visit(eid, set())
    return {**manifest, "project_root": str(root), "ledger_dir": str(ledger),
            "manifest_path": str(path), "manifest_sha256": digest(path.read_bytes())}


def atomic_json(path, value, *, exclusive=False):
    path = Path(path)
    data = json.dumps(value, sort_keys=True, indent=2, allow_nan=False).encode() + b"\n"
    if exclusive:
        with path.open("xb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        return
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with temporary.open("xb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


class Laboratory:
    def __init__(self, ledger, process_runner=None):
        self.ledger = Path(ledger).resolve()
        self.process_runner = process_runner
        self.state = None
        self.sequence = 0
        self.previous_hash = "0" * 64
        self.interrupted = False

    @contextlib.contextmanager
    def locked(self):
        self.ledger.mkdir(parents=True, exist_ok=True)
        with (self.ledger / ".lock").open("a+") as handle:
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise LabBusy("another runner owns this ledger; no concurrent experiments started") from error
            try:
                yield self
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)

    def journal(self):
        path = self.ledger / "journal.jsonl"
        records = []
        previous = "0" * 64
        if path.exists():
            for number, line in enumerate(path.read_bytes().splitlines(), 1):
                try:
                    record = strict_json(line)
                    checksum = record.pop("record_sha256")
                except (ValueError, KeyError) as error:
                    raise LabError(f"journal record {number} is incomplete or invalid; preserve evidence and repair manually") from error
                if record.get("sequence") != number or record.get("previous_sha256") != previous or checksum != digest(canonical(record)):
                    raise LabError(f"journal record {number} fails integrity")
                previous = checksum
                records.append({**record, "record_sha256": checksum})
        return records

    def record(self, event, details=None):
        payload = {"sequence": self.sequence + 1, "previous_sha256": self.previous_hash,
                   "time": timestamp(), "event": event, "details": details or {},
                   "state": self.state}
        checksum = digest(canonical(payload))
        data = canonical({**payload, "record_sha256": checksum}) + b"\n"
        descriptor = os.open(self.ledger / "journal.jsonl", os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        try:
            view = memoryview(data)
            while view:
                written = os.write(descriptor, view)
                view = view[written:]
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        self.sequence += 1
        self.previous_hash = checksum
        atomic_json(self.ledger / "state.json", self.state)

    def initialize(self, manifest_path):
        if (self.ledger / "journal.jsonl").exists():
            raise LabError("ledger already initialized; use resume")
        manifest = normalize_manifest(manifest_path)
        if Path(manifest["ledger_dir"]) != self.ledger:
            raise LabError("actual ledger must match the manifest ledger_dir")
        self.ledger.mkdir(parents=True, exist_ok=True)
        self.state = {"schema": "research_ledger_state_v1", "program": manifest,
                      "manual_paused": False, "program_status": "active", "experiments": {},
                      "release_gate": {"approved": False, "evaluation": "external_gate_required"}}
        for spec in manifest["experiments"]:
            self.state["experiments"][spec["id"]] = self.registered_item(spec)
        self.refresh()
        self.record("program_initialized", {"manifest_sha256": manifest["manifest_sha256"]})
        return self.status()

    def load(self, *, recover=True):
        records = self.journal()
        if not records:
            raise LabError("ledger is not initialized")
        self.state = copy.deepcopy(records[-1]["state"])
        self.sequence = records[-1]["sequence"]
        self.previous_hash = records[-1]["record_sha256"]
        if recover:
            recovered = []
            for eid, item in self.state["experiments"].items():
                orphan = self.ledger / "attempts" / eid / "0001"
                if item["status"] == "running" or (item["status"] == "pending" and orphan.exists()):
                    item["status"] = "inconclusive"
                    item["result"] = {"status": "inconclusive", "summary": "runner interrupted; no automatic retry",
                                      "evidence": {"preserved_attempt_directory": str(orphan)},
                                      "next_actions": ["Inspect preserved artifacts and register a new retry experiment id."]}
                    recovered.append(eid)
            if recovered:
                self.refresh()
                self.record("interrupted_attempts_recovered", {"experiments": recovered})
            else:
                atomic_json(self.ledger / "state.json", self.state)
        return self

    def ready_ids(self):
        items = self.state["experiments"]
        return [eid for eid, item in items.items() if item["status"] == "pending" and
                all(items[dep]["status"] == "passed" for dep in item["spec"]["depends_on"])]

    def refresh(self):
        self.state["program_status"] = ("manual_paused" if self.state["manual_paused"] else
                                         "active" if self.ready_ids() else "awaiting_design")

    def status(self):
        self.refresh()
        ready = self.ready_ids()
        next_actions = []
        for eid, item in self.state["experiments"].items():
            if item["status"] == "pending" and eid not in ready:
                unmet = [dep for dep in item["spec"]["depends_on"] if self.state["experiments"][dep]["status"] != "passed"]
                next_actions.append(f"{eid}: redesign or satisfy dependencies {', '.join(unmet)}")
            elif item["status"] in OUTCOMES:
                next_actions.extend(item.get("result", {}).get("next_actions", []))
        if not ready:
            next_actions.append("Register a new hypothesis with a distinct experiment id; existing outcomes remain preserved.")
        return {"program_id": self.state["program"]["program_id"], "program_status": self.state["program_status"],
                "ready": ready, "experiments": {eid: item["status"] for eid, item in self.state["experiments"].items()},
                "next_actions": list(dict.fromkeys(next_actions)), "journal_records": self.sequence,
                "release_gate": self.state["release_gate"]}

    def merge(self, manifest_path):
        new = normalize_manifest(manifest_path)
        old = self.state["program"]
        if (new["program_id"] != old["program_id"] or new["project_root"] != old["project_root"] or
                Path(new["ledger_dir"]) != self.ledger):
            raise LabError("resume cannot change program id, project root, or ledger_dir")
        added = []
        for spec in new["experiments"]:
            item = self.state["experiments"].get(spec["id"])
            if item and canonical(spec) != canonical(item["spec"]):
                raise LabError(f"{spec['id']}: registered experiment is immutable; use a new id")
            if not item:
                self.state["experiments"][spec["id"]] = self.registered_item(spec)
                added.append(spec["id"])
        self.state["program"] = new
        self.refresh()
        self.record("experiments_registered", {"added": added, "manifest_sha256": new["manifest_sha256"]})
        return added

    def registered_item(self, spec):
        item = {"spec": spec, "status": "pending", "attempts": [], "registered_at": timestamp()}
        try:
            identity, fingerprint = self.provenance(spec)
            item.update(registered_identity=identity, registered_fingerprint=fingerprint)
        except (LabError, OSError) as error:
            item["registration_error"] = str(error)
        return item

    def provenance(self, spec):
        root = Path(self.state["program"]["project_root"])
        hashes = {}
        for value in spec["source_paths"] + spec["input_paths"]:
            local = contained(root, value)
            hashes[value] = digest(local.read_bytes())
        expected = spec.get("expected_sha256", {})
        if not isinstance(expected, dict) or any(hashes.get(path) != sha for path, sha in expected.items()):
            raise LabError("expected source/input identity changed")
        runner_sha = digest(Path(__file__).read_bytes())
        process_path = root / "performance_gate.py"
        process_sha = digest(process_path.read_bytes()) if process_path.is_file() else None
        identity = {"command": spec["command"], "hypothesis": spec["hypothesis"],
                    "files_sha256": hashes, "source_paths": spec["source_paths"], "input_paths": spec["input_paths"],
                    "timeout_seconds": spec["timeout_seconds"], "memory_limit_mb": spec["memory_limit_mb"],
                    "cpu": spec["cpu"], "runner_sha256": runner_sha, "process_runner_sha256": process_sha}
        # An experiment-id or fresh output pathname is bookkeeping, not a new
        # scientific experiment. Suppress renamed copies unless retry is explicit.
        artifact_args = {value for path in spec["artifact_paths"]
                         for value in (path, str(contained(root, path)))}
        fingerprint_identity = copy.deepcopy(identity)
        fingerprint_identity.pop("hypothesis")
        fingerprint_identity["command"] = ["<experiment_id>" if index > 0 and spec["command"][index - 1] in ("--experiment-id", "--experiment_id") else
                                            "<artifact>" if arg in artifact_args else arg
                                            for index, arg in enumerate(spec["command"])]
        return identity, digest(canonical(fingerprint_identity))

    def validate(self):
        records = self.journal()
        problems = []
        for eid, item in self.state["experiments"].items():
            for attempt in item["attempts"]:
                prereg_path = self.ledger / attempt["preregistration"]
                if not prereg_path.is_file() or digest(prereg_path.read_bytes()) != attempt["preregistration_sha256"]:
                    problems.append(f"{eid}: preregistration missing or changed")
                result_path = self.ledger / attempt["directory"] / "process.json"
                if "process_sha256" in attempt and (not result_path.is_file() or digest(result_path.read_bytes()) != attempt["process_sha256"]):
                    problems.append(f"{eid}: process evidence missing or changed")
                for filename, expected in attempt.get("evidence_sha256", {}).items():
                    local = self.ledger / attempt["directory"] / filename
                    if not local.is_file() or digest(local.read_bytes()) != expected:
                        problems.append(f"{eid}: preserved evidence missing or changed: {filename}")
            try:
                current, _ = self.provenance(item["spec"])
                if current != item.get("registered_identity"):
                    problems.append(f"{eid}: source/input/runner identity differs from immutable registration")
                for attempt in item["attempts"]:
                    prereg_path = self.ledger / attempt["preregistration"]
                    if prereg_path.is_file() and strict_json(prereg_path.read_bytes())["identity"] != current:
                        problems.append(f"{eid}: current identity differs from immutable attempt preregistration")
            except (ValueError, KeyError, OSError) as error:
                problems.append(f"{eid}: {error}")
        return {"integrity": "passed" if not problems else "failed", "journal_records": len(records),
                "problems": problems, "release_gate": self.state["release_gate"]}

    def _outcome(self, eid, process):
        if self.interrupted or process.get("status") == "timeout":
            return {"status": "inconclusive", "summary": "interrupted or exhausted process budget; no scientific conclusion",
                    "next_actions": [f"Inspect {eid} evidence; register an explicit new experiment for further work."]}
        if process.get("status") != "ok" or process.get("exit_code") != 0:
            return {"status": "error", "summary": "worker process did not complete successfully", "next_actions": [f"Repair {eid} worker and register a new experiment id."]}
        try:
            value = strict_json(process["stdout"])
            if (not isinstance(value, dict) or value.get("schema") != RESULT_SCHEMA or value.get("experiment_id") != eid or
                    value.get("status") not in OUTCOMES or not isinstance(value.get("summary"), str) or not value["summary"].strip() or
                    not isinstance(value.get("evidence"), dict) or not isinstance(value.get("next_actions"), list) or
                    any(not isinstance(v, str) for v in value["next_actions"])):
                raise ValueError("worker result contract mismatch")
            allowed = {"schema", "experiment_id", "status", "summary", "evidence", "next_actions"}
            if set(value) != allowed:
                raise ValueError("worker result keys must match the scientific result contract exactly")
        except (ValueError, KeyError, TypeError) as error:
            return {"status": "error", "summary": f"invalid scientific result: {error}",
                    "next_actions": [f"Repair {eid} result contract and register a new experiment id."]}
        return value

    def run_one(self, eid, timeout, *, deadline=None):
        item = self.state["experiments"][eid]
        spec = item["spec"]
        try:
            identity, fingerprint = self.provenance(spec)
            if identity != item.get("registered_identity"):
                raise LabError("source/input/runner identity differs from immutable registration")
            root = Path(self.state["program"]["project_root"])
            occupied = [value for value in spec["artifact_paths"] if contained(root, value).exists()]
            if occupied:
                raise LabError(f"worker output already exists; refusing overwrite: {', '.join(occupied)}")
        except (OSError, LabError) as error:
            item["status"] = "blocked"
            item["result"] = {"status": "blocked", "summary": str(error), "next_actions": [f"Repair {eid} provenance in a new registered experiment."]}
            self.refresh()
            self.record("provenance_blocked", {"experiment_id": eid, "reason": str(error)})
            return
        previous = [(old_id, attempt) for old_id, old in self.state["experiments"].items()
                    for attempt in old["attempts"] if attempt["fingerprint"] == fingerprint]
        if previous and spec.get("retry_of") not in {old_id for old_id, _ in previous}:
            item["status"] = "blocked"
            item["result"] = {"status": "blocked", "summary": "duplicate experiment fingerprint; automatic no-op suppressed",
                              "next_actions": [f"Change {eid} hypothesis or implementation; explicit retry_of and retry_reason are required to repeat evidence."]}
            self.refresh()
            self.record("duplicate_suppressed", {"experiment_id": eid, "fingerprint": fingerprint})
            return
        directory = self.ledger / "attempts" / eid / "0001"
        try:
            directory.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            item["status"] = "inconclusive"
            item["result"] = {"status": "inconclusive", "summary": "orphan attempt directory preserved; refusing overwrite",
                              "next_actions": [f"Inspect {eid} orphan attempt and register an explicit new experiment."]}
            self.refresh()
            self.record("orphan_attempt_preserved", {"experiment_id": eid, "directory": str(directory)})
            return
        prereg = {"schema": "experiment_preregistration_v1", "experiment_id": eid, "time": timestamp(),
                  "spec": spec, "identity": identity, "fingerprint": fingerprint,
                  "effective_timeout_seconds": timeout, "measurement_scope": "development diagnostic; external isolation and release gate required"}
        prereg_path = directory / "preregistration.json"
        atomic_json(prereg_path, prereg, exclusive=True)
        attempt = {"directory": str(directory.relative_to(self.ledger)),
                   "preregistration": str(prereg_path.relative_to(self.ledger)),
                   "preregistration_sha256": digest(prereg_path.read_bytes()), "fingerprint": fingerprint}
        item["attempts"].append(attempt)
        item["status"] = "running"
        self.record("experiment_started", {"experiment_id": eid, "fingerprint": fingerprint})
        root = Path(self.state["program"]["project_root"])
        artifacts_before = {value: None for value in spec["artifact_paths"]}
        runner = self.process_runner
        if runner is None:
            if str(root) not in sys.path:
                sys.path.insert(0, str(root))
            runner = lambda command, **limits: run_cold_process(command, **limits, cancelled=lambda: self.interrupted)
        command = [sys.executable, *spec["command"][1:]]
        launch_timeout = min(timeout, deadline - time.monotonic()) if deadline is not None else timeout
        try:
            if launch_timeout <= .05:
                process = {"status": "timeout", "command": command, "stdout": "", "stderr": "",
                           "summary": "cycle wall budget exhausted during preregistration; worker was not launched",
                           "effective_timeout_seconds": 0}
            else:
                process = runner(command, timeout=launch_timeout, workdir=str(root), memory_mib=spec["memory_limit_mb"], cpu=spec["cpu"])
                process["effective_timeout_seconds"] = launch_timeout
        except KeyboardInterrupt:
            self.interrupted = True
            process = {"status": "interrupted", "command": command, "stdout": "", "stderr": "runner interrupted"}
        except Exception as error:
            process = {"status": "runner_error", "command": command, "stdout": "", "stderr": f"{type(error).__name__}: {error}"}
        atomic_json(directory / "process.json", process, exclusive=True)
        (directory / "stdout.txt").write_text(process.get("stdout", ""))
        (directory / "stderr.txt").write_text(process.get("stderr", ""))
        attempt["process_sha256"] = digest((directory / "process.json").read_bytes())
        artifact_records = []
        artifact_errors = []
        for index, value in enumerate(spec["artifact_paths"]):
            try:
                local = contained(root, value)
                if local.is_file():
                    data = local.read_bytes()
                    saved = directory / f"artifact-{index:03d}-{local.name}"
                    saved.write_bytes(data)
                    artifact_records.append({"path": value, "sha256": digest(data), "before_sha256": artifacts_before[value],
                                             "snapshot": saved.name, "changed": digest(data) != artifacts_before[value]})
                else:
                    raise LabError(f"declared artifact missing: {value}")
            except (OSError, LabError) as error:
                artifact_records.append({"path": value, "capture_error": str(error)})
                artifact_errors.append(str(error))
        atomic_json(directory / "artifacts.json", artifact_records, exclusive=True)
        result = self._outcome(eid, process)
        if artifact_errors and result["status"] in ("passed", "refuted"):
            result = {"status": "error", "summary": "artifact capture failed: " + "; ".join(artifact_errors),
                      "next_actions": [f"Repair {eid} artifact collection and register a new experiment."]}
        elif artifact_errors:
            result.setdefault("evidence", {})["artifact_capture_errors"] = artifact_errors
        try:
            after_identity, after_fingerprint = self.provenance(spec)
            if after_fingerprint != fingerprint:
                result = {"status": "inconclusive", "summary": "source or input identity changed during execution",
                          "next_actions": [f"Freeze {eid} data and implementation and register a new experiment."]}
        except (OSError, LabError):
            result = {"status": "inconclusive", "summary": "source or input identity unavailable after execution", "next_actions": []}
        item["status"] = result["status"]
        item["result"] = result
        atomic_json(directory / "scientific_result.json", result, exclusive=True)
        attempt["evidence_sha256"] = {local.name: digest(local.read_bytes()) for local in directory.iterdir() if local.is_file()}
        self.refresh()
        self.record("experiment_finished", {"experiment_id": eid, "status": item["status"]})

    def run_ready(self, max_experiments, wall_budget):
        positive(max_experiments, "max_experiments", integer=True)
        positive(wall_budget, "wall_budget")
        start = time.monotonic()
        executed = []
        if self.state["manual_paused"]:
            return {**self.status(), "executed": executed, "budget_exhausted": False}
        while len(executed) < max_experiments and not self.interrupted:
            remaining = wall_budget - (time.monotonic() - start)
            ready = self.ready_ids()
            if remaining <= .05 or not ready:
                break
            eid = ready[0]
            self.run_one(eid, min(remaining, self.state["experiments"][eid]["spec"]["timeout_seconds"]), deadline=start + wall_budget)
            executed.append(eid)
        self.refresh()
        exhausted = time.monotonic() - start >= wall_budget - .05
        self.record("bounded_cycle_finished", {"executed": executed, "budget_exhausted": exhausted, "interrupted": self.interrupted})
        return {**self.status(), "executed": executed, "budget_exhausted": exhausted}

    def daemon(self, cycles, interval, max_experiments, wall_budget):
        positive(cycles, "cycles", integer=True)
        if isinstance(interval, bool) or not isinstance(interval, (int, float)) or not math.isfinite(interval) or not 0 <= interval <= 60:
            raise LabError("interval must be finite, between zero and 60 seconds")
        positive(wall_budget, "wall_budget")
        positive(max_experiments, "max_experiments", integer=True)
        start = time.monotonic()
        results = []
        for _ in range(cycles):
            remaining = wall_budget - (time.monotonic() - start)
            if remaining <= .05 or self.interrupted:
                break
            results.append(self.run_ready(max_experiments, remaining))
            if not self.ready_ids() or self.state["manual_paused"]:
                break
            delay = min(interval, max(0, wall_budget - (time.monotonic() - start)))
            if delay:
                wake_at = time.monotonic() + delay
                while not self.interrupted and time.monotonic() < wake_at:
                    time.sleep(max(0, min(.1, wake_at - time.monotonic())))
        return {**self.status(), "completed_cycles": len(results), "cycle_results": results,
                "budget_exhausted": time.monotonic() - start >= wall_budget - .05}


def cli(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    init = sub.add_parser("init")
    init.add_argument("--manifest", required=True)
    init.add_argument("--ledger")
    for action in ("status", "validate", "resume", "pause", "run-ready", "daemon"):
        child = sub.add_parser(action)
        child.add_argument("--ledger", default="lab/ledger")
        if action == "resume":
            child.add_argument("--manifest")
        if action == "pause":
            child.add_argument("--reason", default="manual pause")
        if action in ("run-ready", "daemon"):
            child.add_argument("--max-experiments", type=int, default=1)
            child.add_argument("--wall-budget", type=float, required=True)
        if action == "daemon":
            child.add_argument("--cycles", type=int, required=True)
            child.add_argument("--interval", type=float, default=5)
    args = parser.parse_args(argv)
    try:
        ledger = args.ledger
        if args.action == "init" and ledger is None:
            ledger = normalize_manifest(args.manifest)["ledger_dir"]
        lab = Laboratory(ledger)
        old_handlers = {}
        if args.action in ("run-ready", "daemon"):
            def interrupt(signum, frame):
                lab.interrupted = True
            for signum in (signal.SIGINT, signal.SIGTERM):
                old_handlers[signum] = signal.signal(signum, interrupt)
        try:
            with lab.locked():
                if args.action == "init":
                    output = lab.initialize(args.manifest)
                else:
                    lab.load()
                    if args.action == "status":
                        output = lab.status()
                    elif args.action == "validate":
                        output = lab.validate()
                    elif args.action == "resume":
                        if args.manifest:
                            lab.merge(args.manifest)
                        lab.state["manual_paused"] = False
                        lab.refresh()
                        lab.record("program_resumed")
                        output = lab.status()
                    elif args.action == "pause":
                        lab.state["manual_paused"] = True
                        lab.refresh()
                        lab.record("program_paused", {"reason": args.reason})
                        output = lab.status()
                    elif args.action == "run-ready":
                        output = lab.run_ready(args.max_experiments, args.wall_budget)
                    else:
                        output = lab.daemon(args.cycles, args.interval, args.max_experiments, args.wall_budget)
        finally:
            for signum, handler in old_handlers.items():
                signal.signal(signum, handler)
        print(json.dumps(output, sort_keys=True, indent=2))
        return 0 if output.get("integrity") != "failed" else 3
    except (ValueError, OSError) as error:
        print(json.dumps({"status": "runner_error", "reason": str(error), "goal_complete": False}))
        return 2


if __name__ == "__main__":
    raise SystemExit(cli())
