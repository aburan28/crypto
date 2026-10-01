"""Load, validate and identify an ICMS specification.

A spec is validated in three layers, and any failure refuses the spec:

1. the JSON Schema in docs/ic/measurement/schema/spec.v1.json (structure);
2. the registry (vocabulary: families, units, windows, references, the solver
   options that must be pinned);
3. policy rules that span sections (one target in the primary window, CPUs
   at least threads, a wall-clock budget flagged, ...).

Defaults are written out before hashing, the way taskq normalises its specs,
so two files that differ only by an omitted default have the same spec_id.
"""
from __future__ import annotations

import copy
import json
import os
from typing import Any

from . import SPEC_SCHEMA
from .canonical import spec_id as _spec_id, workload_id as _workload_id
from .registry import Registry
from .schemacheck import Validator

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SCHEMA_DIR = os.path.join(REPO, "docs", "ic", "measurement", "schema")

NON_IDENTITY = ("label", "description", "role")
NON_IDENTITY_NESTED = (("execution", "notes"),)

DEFAULTS = {
    "factor_base": {"quotient": [], "large_primes": {"mode": "none"}, "params": {}},
    "decomposition": {"encoding": "none", "splitting": "none", "params": {}},
    "descent": {"method": "embedded"},
    "measurement": {"warmup": 1, "interleave": "abab", "isolation_required": "L2", "sample_interval_ms": 50,
                    "thresholds": {}},
    "execution": {"env": {}, "args": {}},
    "extensions": {},
}


class SpecError(ValueError):
    def __init__(self, problems: list[str]):
        super().__init__("; ".join(problems))
        self.problems = problems


def load_document(path: str) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    if path.endswith((".yaml", ".yml")):
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise SpecError([f"{path}: reading YAML needs PyYAML; install it or use JSON"]) from exc
        doc = yaml.safe_load(text)
    else:
        doc = json.loads(text)
    if not isinstance(doc, dict):
        raise SpecError([f"{path}: top level must be a mapping"])
    return doc


def _merge_defaults(spec: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(spec)
    for section, defaults in DEFAULTS.items():
        if section not in out:
            if section in ("descent", "extensions"):
                out[section] = copy.deepcopy(defaults)
            continue
        if isinstance(defaults, dict):
            for k, v in defaults.items():
                out[section].setdefault(k, copy.deepcopy(v))
    out["execution"].setdefault("cpus", out["execution"]["threads"])
    return out


def _exact(obj: Any) -> Any:
    """Floats become their shortest round-trip decimal string, so identity never hashes a float."""
    if isinstance(obj, float):
        return repr(obj)
    if isinstance(obj, dict):
        return {k: _exact(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_exact(v) for v in obj]
    return obj


def identity_view(spec: dict[str, Any]) -> dict[str, Any]:
    view = {k: _exact(v) for k, v in spec.items() if k not in NON_IDENTITY}
    for section, key in NON_IDENTITY_NESTED:
        if section in view and key in view[section]:
            view[section] = {k: v for k, v in view[section].items() if k != key}
    return view


def _policy(spec: dict[str, Any], reg: Registry) -> list[str]:
    p: list[str] = []
    fb, dec, meas, ex = spec["factor_base"], spec["decomposition"], spec["measurement"], spec["execution"]
    fam = reg.family(fb["family"])
    if fam is None:
        p.append(f"factor_base.family {fb['family']!r} is not in the registry")
    else:
        missing = [k for k in fam.get("params", []) if k not in fb.get("params", {})]
        if missing:
            p.append(f"factor_base.params is missing {missing} for family {fb['family']}")
        for k, allowed in (fam.get("param_values") or {}).items():
            v = (fb.get("params") or {}).get(k)
            if v is not None and v not in allowed:
                p.append(f"factor_base.params.{k} {v!r} is not one of {allowed}")
    if not reg.has_unit(spec["accounting"]["unit"]):
        p.append(f"accounting.unit {spec['accounting']['unit']!r} is not in the registry")
    window = reg.window(meas["window"])
    if window is None:
        p.append(f"measurement.window {meas['window']!r} is not in the registry")
    elif window["name"] == "online_one_target" and spec["instance"]["workload"]["targets"] != 1:
        p.append("the online_one_target window solves exactly one target (AGENTS.md 'Primary comparison uses one target')")
    ref = spec.get("reference") or {}
    if ref.get("rho") and not reg.has_rho(ref["rho"]):
        p.append(f"reference.rho {ref['rho']!r} is not in the registry")
    if ref.get("floor") and not reg.has_floor(ref["floor"]):
        p.append(f"reference.floor {ref['floor']!r} is not in the registry")
    solver = dec.get("solver")
    if solver:
        need = reg.solver_options(solver["name"])
        have = solver.get("options") or {}
        missing = [o for o in need if o not in have]
        if missing:
            p.append(f"decomposition.solver.options must pin {missing} for {solver['name']} "
                     "(registry: solver_options_that_change_search)")
    if ex["cpus"] < ex["threads"]:
        p.append(f"execution.cpus ({ex['cpus']}) is fewer than execution.threads ({ex['threads']})")
    if spec["relations"]["stop"] == "count" and "count" not in spec["relations"]:
        p.append("relations.stop = count needs relations.count")
    return p


def validate(spec: dict[str, Any], reg: Registry | None = None) -> dict[str, Any]:
    """Return the normalised spec, or raise SpecError with every problem found."""
    reg = reg or Registry.load()
    with open(os.path.join(SCHEMA_DIR, "spec.v1.json"), encoding="utf-8") as fh:
        schema = json.load(fh)
    problems = Validator(schema).errors(spec)
    if problems:
        raise SpecError(problems)
    norm = _merge_defaults(spec)
    problems = _policy(norm, reg)
    if problems:
        raise SpecError(problems)
    return norm


def spec_id(norm: dict[str, Any]) -> str:
    return _spec_id(identity_view(norm))


def workloads(norm: dict[str, Any]) -> list[dict[str, Any]]:
    """One workload record per seed, in the cryptanalysis AGENTS.md sense:
    curve/subgroup, targets, input law, seed and cold target count."""
    inst = norm["instance"]
    out = []
    for seed in inst["workload"]["seeds"]:
        rec = {
            "curve": inst["curve"],
            "targets": inst["workload"]["targets"],
            "law": inst["workload"]["law"],
            "seed": seed,
            "cold_target_count": inst["workload"]["targets"],
            "file": inst["workload"].get("file"),
        }
        out.append({"workload_id": _workload_id(rec), "record": rec})
    return out


def load(path: str, reg: Registry | None = None) -> dict[str, Any]:
    raw = load_document(path)
    if raw.get("icms") != SPEC_SCHEMA:
        raise SpecError([f"{path}: icms must be {SPEC_SCHEMA!r}"])
    norm = validate(raw, reg)
    return {"path": os.path.relpath(path, REPO) if os.path.abspath(path).startswith(REPO) else path,
            "spec": norm, "spec_id": spec_id(norm), "workloads": workloads(norm)}
