"""Re-derive a session from its files and refuse anything that does not match.

A committed session is evidence only if a reader can recompute it.  For every
session directory this checks, without running anything:

* session.json, capsule.json, every record and every frozen comparison
  validates against its schema in docs/ic/measurement/schema/;
* the capsule, records.jsonl and every spec copy hash to what session.json says;
* the capsule's env_class is the env_class of its stable facts and hashes to
  env_class_id, and every record names that class and capsule;
* every spec copy loads under the current standard and hashes to the spec id
  and workload ids its arms and records carry;
* records follow the plan in order, carry the session id, and hash to their
  own record_id; every stdout/stderr file under exec/ hashes to the record;
* the isolation level of every record is what the gate computes from the
  record's own raw observations, the session and the capsule;
* every vocabulary key (unit, window, phase) is in the registry;
* every frozen comparison is exactly what ``compare`` computes now from the
  same records and request.

Anything that fails is a problem to fix by producing a new session, never by
editing the old one (sessions are immutable evidence).
"""
from __future__ import annotations

import json
import os
from typing import Any

from . import CAPSULE_SCHEMA, COMPARISON_SCHEMA, RECORD_SCHEMA, SESSION_SCHEMA
from .canonical import env_class_id, record_id, sha256_file
from .environment import env_class
from .gates import evaluate, level_at_least
from .registry import Registry
from .schemacheck import Validator
from .spec import REPO, SpecError, load as load_spec

SCHEMA_DIR = os.path.join(REPO, "docs", "ic", "measurement", "schema")
SCHEMA_FILES = {SESSION_SCHEMA: "session.v1.json", RECORD_SCHEMA: "record.v1.json",
                CAPSULE_SCHEMA: "capsule.v1.json", COMPARISON_SCHEMA: "comparison.v1.json"}


def validators() -> dict[str, Validator]:
    out = {}
    for sid, name in SCHEMA_FILES.items():
        with open(os.path.join(SCHEMA_DIR, name)) as fh:
            out[sid] = Validator(json.load(fh))
    return out


def _load_json(path: str) -> Any:
    with open(path) as fh:
        return json.load(fh)


def _normalise(obj: Any) -> Any:
    """JSON round trip, so tuples and lists compare equal as they would on disk."""
    return json.loads(json.dumps(obj, sort_keys=True))


def find_sessions(path: str) -> list[str]:
    """A session directory, or a directory whose children are session directories."""
    if os.path.isfile(os.path.join(path, "session.json")):
        return [path]
    if not os.path.isdir(path):
        return []
    return sorted(os.path.join(path, d) for d in os.listdir(path)
                  if os.path.isfile(os.path.join(path, d, "session.json")))


def audit_session(path: str, reg: Registry | None = None, vals: dict[str, Validator] | None = None) -> list[str]:
    reg = reg or Registry.load()
    vals = vals or validators()
    problems: list[str] = []

    def bad(msg: str) -> None:
        problems.append(f"{os.path.basename(os.path.normpath(path))}: {msg}")

    for need in ("session.json", "capsule.json", "records.jsonl", "specs"):
        if not os.path.exists(os.path.join(path, need)):
            bad(f"missing {need}")
    if problems:
        return problems
    session = _load_json(os.path.join(path, "session.json"))
    for e in vals[SESSION_SCHEMA].errors(session):
        bad(f"session.json {e}")
    if problems:
        return problems

    # ---- capsule ---------------------------------------------------------------
    cap_path = os.path.join(path, "capsule.json")
    cap = _load_json(cap_path)
    for e in vals[CAPSULE_SCHEMA].errors(cap):
        bad(f"capsule.json {e}")
    if sha256_file(cap_path) != session["capsule_sha256"]:
        bad("capsule.json does not hash to session.capsule_sha256")
    if not problems:
        derived = _normalise(env_class(cap["stable"]))
        if derived != cap["env_class"]:
            bad("capsule env_class is not the env_class of its stable facts")
        if env_class_id(cap["env_class"]) != cap["env_class_id"]:
            bad("capsule env_class_id is not the hash of its env_class")
        if cap["env_class_id"] != session["env_class_id"]:
            bad("session env_class_id differs from the capsule's")

    # ---- specs -----------------------------------------------------------------
    loaded: dict[str, dict[str, Any]] = {}
    on_disk = sorted(os.listdir(os.path.join(path, "specs")))
    if on_disk != sorted(session["spec_files"]):
        bad(f"specs/ holds {on_disk}, session lists {sorted(session['spec_files'])}")
    for name, digest in session["spec_files"].items():
        sp = os.path.join(path, "specs", name)
        if sha256_file(sp) != digest:
            bad(f"specs/{name} does not hash to session.spec_files")
            continue
        try:
            loaded[name] = load_spec(sp, reg)
        except SpecError as exc:
            bad(f"specs/{name} no longer validates: {exc.problems}")
    arms = {a["index"]: a for a in session["arms"]}
    if sorted(arms) != list(range(len(session["arms"]))):
        bad("arm indices are not 0..n-1")
    for a in session["arms"]:
        sp = loaded.get(os.path.basename(a["spec_path"]))
        if sp is None:
            bad(f"arm {a['index']}: no spec copy for {a['spec_path']}")
            continue
        if sp["spec_id"] != a["spec_id"]:
            bad(f"arm {a['index']}: spec copy hashes to {sp['spec_id']}, arm says {a['spec_id']}")
        if a["workload_id"] not in [w["workload_id"] for w in sp["workloads"]]:
            bad(f"arm {a['index']}: workload {a['workload_id']} is not one of its spec's workloads")

    # ---- records ---------------------------------------------------------------
    rec_path = os.path.join(path, "records.jsonl")
    if sha256_file(rec_path) != session["records_sha256"]:
        bad("records.jsonl does not hash to session.records_sha256")
    with open(rec_path) as fh:
        records = [json.loads(line) for line in fh if line.strip()]
    if len(records) != len(session["plan"]):
        bad(f"{len(records)} records for a plan of {len(session['plan'])} executions")
    units = {u["id"] for u in reg.doc["units"]["list"]}
    windows = {w["name"] for w in reg.doc["windows"]["list"]}
    phases = set(reg.phases)
    for n, rec in enumerate(records):
        where = f"record {n}"
        errs = vals[RECORD_SCHEMA].errors(rec)
        for e in errs:
            bad(f"{where} {e}")
        if errs:
            continue
        body = {k: v for k, v in rec.items() if k != "record_id"}
        if record_id(body) != rec["record_id"]:
            bad(f"{where}: content does not hash to its record_id {rec['record_id']}")
        if rec["session_id"] != session["session_id"]:
            bad(f"{where}: session_id {rec['session_id']} is not this session's")
        if rec["execution_index"] != n:
            bad(f"{where}: execution_index {rec['execution_index']} out of order")
        if n < len(session["plan"]):
            step = session["plan"][n]
            if (rec["arm"], rec["round"], rec["warmup"]) != (step["arm"], step["round"], step["warmup"]):
                bad(f"{where}: (arm, round, warmup) differs from plan step {step}")
        arm = arms.get(rec["arm"])
        if arm is None:
            bad(f"{where}: arm {rec['arm']} is not in the session")
            continue
        if (rec["spec_id"], rec["workload_id"]) != (arm["spec_id"], arm["workload_id"]):
            bad(f"{where}: spec/workload differ from arm {rec['arm']}")
        if rec["environment"] != {"env_class_id": session["env_class_id"], "capsule_sha256": session["capsule_sha256"]}:
            bad(f"{where}: environment does not name this session's capsule")
        if not set(rec["execution"]["pinned_cpus"]) <= set(session["reservation"]["cpus"]):
            bad(f"{where}: pinned outside the session's reservation")
        exec_dir = os.path.join(path, "exec", f"{n:04d}")
        for stream in ("stdout", "stderr"):
            out = rec["execution"]["outputs"][stream]
            fp = os.path.join(exec_dir, out["path"])
            if os.path.exists(fp):
                if sha256_file(fp) != out["sha256"]:
                    bad(f"{where}: exec/{n:04d}/{out['path']} does not hash to the record")
            elif out["bytes"]:
                bad(f"{where}: exec/{n:04d}/{out['path']} is missing ({out['bytes']} bytes recorded)")
        sp = loaded.get(os.path.basename(arm["spec_path"]))
        if sp is not None:
            check_gate(rec, session, cap, sp["spec"], bad, where)
        if rec["unit"] not in units:
            bad(f"{where}: unit {rec['unit']} is not in the registry")
        if rec["window"] not in windows:
            bad(f"{where}: window {rec['window']} is not in the registry")
        for k in (rec.get("units") or {}):
            if k not in units:
                bad(f"{where}: units.{k} is not in the registry")
        for k, w in (rec.get("windows") or {}).items():
            if k not in windows:
                bad(f"{where}: windows.{k} is not in the registry")
            if w.get("ops_unit") is not None and w["ops_unit"] not in units:
                bad(f"{where}: windows.{k}.ops_unit {w['ops_unit']} is not in the registry")
        for k, ph in (rec.get("phases") or {}).items():
            if k not in phases:
                bad(f"{where}: phases.{k} is not in the registry")
            if ph and ph.get("ops_unit") is not None and ph["ops_unit"] not in units:
                bad(f"{where}: phases.{k}.ops_unit {ph['ops_unit']} is not in the registry")

    # ---- comparisons -----------------------------------------------------------
    cmp_dir = os.path.join(path, "comparisons")
    if os.path.isdir(cmp_dir):
        from .compare import compare
        for name in sorted(os.listdir(cmp_dir)):
            if not name.endswith(".json"):
                bad(f"comparisons/{name}: not a .json file")
                continue
            frozen = _load_json(os.path.join(cmp_dir, name))
            errs = vals[COMPARISON_SCHEMA].errors(frozen)
            for e in errs:
                bad(f"comparisons/{name} {e}")
            if errs:
                continue
            req = frozen["request"]
            again = _normalise(compare(path, req["a"], req["b"], req["declared"]))
            if again != frozen:
                diff = sorted(k for k in set(again) | set(frozen) if again.get(k) != frozen.get(k))
                bad(f"comparisons/{name}: recomputing from the records gives a different result in {diff}")
    return problems


def check_gate(rec: dict[str, Any], session: dict[str, Any], cap: dict[str, Any], spec: dict[str, Any],
               bad, where: str) -> None:
    """Recompute the isolation gate from the record's raw observations."""
    gate = evaluate(rec["execution"], session, cap["stable"], spec["execution"]["threads"],
                    spec["measurement"].get("thresholds"))
    required = spec["measurement"]["isolation_required"]
    want = _normalise({**gate, "required": required, "wall_admissible": level_at_least(gate["earned_level"], required)})
    if want != rec["isolation"]:
        keys = sorted(k for k in want if want.get(k) != rec["isolation"].get(k))
        bad(f"{where}: isolation differs from the gate recomputed from its observations in {keys}")


def audit(paths: list[str]) -> tuple[int, list[str]]:
    reg, vals = Registry.load(), validators()
    sessions: list[str] = []
    for p in paths:
        sessions += find_sessions(p)
    problems: list[str] = []
    for s in sessions:
        problems += audit_session(s, reg, vals)
    return len(sessions), problems
