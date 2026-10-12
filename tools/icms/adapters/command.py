"""The generic adapter: any producer that writes an ICMS metrics file.

The spec gives the argv in ``execution.args.argv``.  ``{seed}``, ``{targets}``
and ``{metrics_path}`` in any argument are replaced per execution.  The
command must write a JSON object to ``$ICMS_METRICS_PATH`` with any of these
keys, each in the record's own shape:

    outcome   {"status": "complete|insufficient_relations|budget|timeout|oom|error", "verified": bool, ...}
    units     {"<unit id>": {"total": ..., "S": ..., "deterministic": bool, "host_dependent_because": [...]}}
    reference {...}, metrics {...}, phases {...}, windows {...}, consistency [...]

A missing file or a missing key is recorded as unknown, never as zero.  This
is how a new producer joins the standard without an adapter of its own: the
cryptanalysis and autoresearcher adapters are thin layers over it.
"""
from __future__ import annotations

import json
import os
from typing import Any

from . import Context


class Command:
    name = "command"

    def check(self, spec: dict[str, Any]) -> list[str]:
        argv = (spec["execution"].get("args") or {}).get("argv")
        if not argv or not isinstance(argv, list) or not all(isinstance(a, str) for a in argv):
            return ["the command adapter needs execution.args.argv, a list of strings"]
        return []

    def metrics_path(self, ctx: Context, workload: dict[str, Any]) -> str:
        return os.path.join(ctx["exec_dir"], "metrics.json")

    def command(self, spec: dict[str, Any], workload: dict[str, Any], ctx: Context) -> dict[str, Any]:
        mpath = self.metrics_path(ctx, workload)
        subst = {"seed": str(workload["record"]["seed"]), "targets": str(workload["record"]["targets"]),
                 "metrics_path": mpath, "repo_root": ctx["repo_root"]}
        argv = [a.format(**subst) for a in spec["execution"]["args"]["argv"]]
        cwd = spec["execution"]["args"].get("cwd", ".")
        cwd = cwd if os.path.isabs(cwd) else os.path.join(ctx["repo_root"], cwd)
        env = dict(spec["execution"].get("env") or {})
        env["ICMS_METRICS_PATH"] = mpath
        return {"argv": argv, "cwd": cwd, "env": env}

    def provenance(self, spec: dict[str, Any], ctx: Context) -> dict[str, Any]:
        return {"binary": spec["execution"].get("binary")}

    def parse(self, spec: dict[str, Any], workload: dict[str, Any], ctx: Context, stdout_path: str) -> dict[str, Any]:
        mpath = self.metrics_path(ctx, workload)
        if not os.path.exists(mpath):
            return {"outcome": {"status": "error", "verified": False, "reason": "no metrics file written"}}
        try:
            with open(mpath, encoding="utf-8") as fh:
                doc = json.load(fh)
        except ValueError as exc:
            return {"outcome": {"status": "error", "verified": False, "reason": f"metrics file unparseable: {exc}"}}
        out = {k: doc.get(k) for k in ("outcome", "units", "reference", "metrics", "phases", "windows", "consistency", "producer")}
        if not isinstance(out.get("outcome"), dict):
            out["outcome"] = {"status": "error", "verified": False, "reason": "metrics file has no outcome"}
        return out
