"""Adapters: compile a spec into a producer's command, and parse what it reports.

An adapter does three things and nothing else:

* ``check(spec)`` returns every reason it cannot run the spec faithfully
  (a family it has no spelling for, a window its producer cannot measure, a
  unit it does not emit).  A non-empty list refuses the spec before any run.
* ``command(spec, workload, ctx)`` returns the argv, working directory and
  declared environment for one execution.  No shell is involved.
* ``parse(spec, workload, ctx, stdout_path)`` turns the producer's report into
  the record's ``outcome``, ``metrics``, ``phases`` and ``windows``.  Anything
  the producer did not report is ``None``; a declared-versus-observed
  mismatch (for example a quotient the spec declares but the report does not
  show) is listed in ``consistency`` and refuses the record's comparisons.
"""
from __future__ import annotations

from typing import Any, Protocol


class Context(dict):
    """repo_root, run_dir, binary paths: whatever the session resolved."""


class Adapter(Protocol):
    name: str

    def check(self, spec: dict[str, Any]) -> list[str]: ...

    def command(self, spec: dict[str, Any], workload: dict[str, Any], ctx: Context) -> dict[str, Any]: ...

    def parse(self, spec: dict[str, Any], workload: dict[str, Any], ctx: Context, stdout_path: str) -> dict[str, Any]: ...

    def provenance(self, spec: dict[str, Any], ctx: Context) -> dict[str, Any]: ...


def get(name: str) -> "Adapter":
    if name == "crypto.ic_bench":
        from .crypto_ic_bench import CryptoIcBench
        return CryptoIcBench()
    if name == "command":
        from .command import Command
        return Command()
    if name == "cryptanalysis.ic_bench":
        from .cryptanalysis_ic_bench import CryptanalysisIcBench
        return CryptanalysisIcBench()
    if name == "autoresearcher.index_calculus":
        from .autoresearcher_ic import AutoresearcherIc
        return AutoresearcherIc()
    raise KeyError(f"no adapter {name!r}")


def none_if_missing(d: dict[str, Any] | None, *path: str) -> Any:
    cur: Any = d
    for p in path:
        if not isinstance(cur, dict) or p not in cur:
            return None
        cur = cur[p]
    return cur


# The record's metric sections and the descriptor keys each may carry besides
# the registry's size fields (factor_base) and PDP fields (system).  Anything
# else a producer reports is kept, under ``native`` (``extra`` for the
# solver), so the standard's field names never silently change meaning.
_SECTIONS = ("instance", "factor_base", "decomposition", "system", "solver", "linear_algebra")
_DESCRIPTORS = {
    "factor_base": ("family", "producer_name", "points_per_column", "native"),
    "decomposition": ("oracle", "arity", "targets_tried", "relations_found", "hit_rate", "native"),
    "system": ("note", "native"),
    "solver": ("name", "calls", "ops", "op_unit", "priced_by", "wall_ns", "budget_exceeded", "solving_degree_mean",
               "solving_degree_max", "peak_bytes", "conflicts", "decisions", "propagations", "restarts",
               "matrix_rows_max", "matrix_cols_max", "extra"),
    "linear_algebra": ("method", "rows", "rank", "dependent", "work", "work_unit", "native"),
}
_SOLVER_REQUIRED = ("name", "calls", "ops", "op_unit", "priced_by", "wall_ns")


def normalise(parsed: dict[str, Any], spec: dict[str, Any], registry: Any) -> dict[str, Any]:
    """Put an adapter's parse result into the record's shape.

    Every registry field of a reported section is present (``None`` when the
    producer did not report it); unknown keys move under ``native``; a unit
    that does not say it is deterministic is not, and says why.
    """
    out = {k: parsed.get(k) for k in ("outcome", "units", "reference", "metrics", "phases", "windows",
                                      "consistency", "producer")}
    out["consistency"] = list(out["consistency"] or [])
    if not isinstance(out["outcome"], dict):
        out["outcome"] = {"status": "error", "verified": False, "reason": "the producer reported no outcome"}
    out["outcome"].setdefault("verified", False)
    units = out["units"]
    if isinstance(units, dict):
        for uid, u in units.items():
            if not isinstance(u, dict):
                units[uid] = u = {"total": None}
            u.setdefault("total", None)
            u.setdefault("deterministic", False)
            u.setdefault("host_dependent_because", [])
            if not u["deterministic"] and not u["host_dependent_because"]:
                u["host_dependent_because"] = ["the producer did not state that this unit is deterministic"]
    metrics = out["metrics"]
    if isinstance(metrics, dict):
        fields = {
            "factor_base": list(registry.doc["factor_base_size_fields"]["fields"]),
            "system": list(registry.doc["pdp_metrics"]["fields"]),
        }
        # Search-effort fields of registry.pdp_metrics describe the solver's run,
        # not the system it was handed: the record keeps them under solver.
        solver_registry = {"conflicts", "decisions", "propagations", "restarts", "matrix_rows_max", "matrix_cols_max"}
        system = metrics.get("system")
        if isinstance(system, dict) and (solver_registry | {"solving_degree"}) & set(system):
            metrics = dict(metrics)
            system = dict(system)
            solver = dict(metrics.get("solver") or {})
            for k in sorted(solver_registry & set(system)):
                solver.setdefault(k, system.pop(k))
            if "solving_degree" in system:
                solver.setdefault("solving_degree_max", system.pop("solving_degree"))
            metrics["system"], metrics["solver"] = system, solver
        moved_to = {"solver": "extra"}
        norm: dict[str, Any] = {}
        for sec in _SECTIONS:
            val = metrics.get(sec)
            if not isinstance(val, dict) or sec == "instance":
                norm[sec] = val if isinstance(val, dict) else None
                continue
            allowed = set(fields.get(sec, [])) | set(_DESCRIPTORS[sec])
            if sec == "system":
                allowed -= solver_registry | {"solving_degree"}
            bucket = moved_to.get(sec, "native")
            kept = {k: v for k, v in val.items() if k in allowed}
            extra = {k: v for k, v in val.items() if k not in allowed}
            if extra:
                kept[bucket] = {**(kept.get(bucket) or {}), **extra}
            for f in fields.get(sec, []):
                if f in allowed:
                    kept.setdefault(f, None)
            if sec == "system" and kept.get("max_degree") is None and kept.get("degrees"):
                kept["max_degree"] = max(kept["degrees"])
            if sec == "factor_base":
                kept.setdefault("family", spec["factor_base"]["family"])
            if sec == "solver":
                for f in _SOLVER_REQUIRED:
                    kept.setdefault(f, None)
            norm[sec] = kept
        unknown = sorted(set(metrics) - set(_SECTIONS))
        if unknown:
            norm.setdefault("instance", None)
            norm["instance"] = {**(norm["instance"] or {}), "native_sections": {k: metrics[k] for k in unknown}}
        out["metrics"] = norm
    return out
