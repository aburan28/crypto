"""Compact views of worker records and results, shared by the MCP server and the CLI (no mcp dependency)."""
from __future__ import annotations

from typing import Any


def worker_line(w: dict[str, Any]) -> dict[str, Any]:
    cap = w.get("capacity") or {}
    return {"id": w["id"], "pools": w.get("pools"), "labels": w.get("labels"), "summary": w.get("summary"),
            "lab_cpus": w.get("lab_cpus"), "max_cpus_single_node": {k: v.get("single_node_max") for k, v in cap.items()},
            "max_tier": w.get("max_tier"), "backends": w.get("backends"), "oci_runtimes": w.get("oci_runtimes"),
            "images": sorted({n for img in w.get("images") or [] for n in (img.get("names") or [])})[:20],
            "gpus": [g.get("name") for g in w.get("gpus") or []], "busy": w.get("busy"),
            "heartbeat_age_s": w.get("heartbeat_age_s")}


def summary_view(res: dict[str, Any]) -> dict[str, Any]:
    pl = res.get("placement") or {}
    iso = pl.get("isolation") or {}
    return {
        "job_id": res["job_id"], "status": res["status"], "outcome_class": res["outcome_class"], "error": res.get("error"),
        "worker": res["worker"], "attempt": res["attempt"],
        "fidelity": {k: res["fidelity"].get(k) for k in ("policy", "tier", "grade", "contended", "violations")},
        "placement": {"cpus": pl.get("cpus"), "idle_siblings": pl.get("idle_siblings"), "nodes": pl.get("nodes"),
                      "partition": iso.get("partition_state"), "threads_moved": iso.get("threads_moved"),
                      "irqs_moved": iso.get("irqs_moved"), "gpus": pl.get("gpus")},
        "runtime": {k: (res.get("runtime") or {}).get(k) for k in ("backend", "oci_runtime", "image", "image_digest")},
        "summary": res.get("summary"), "verification": res.get("verification"),
        "timing": res.get("timing"), "notes": res.get("notes"),
        "runs": [{"index": r["index"], "warmup": r["warmup"], "exit_code": r["exit_code"], "wall_s": r["wall_s"],
                  "instructions": (r.get("counters") or {}).get("instructions"), "cycles": (r.get("counters") or {}).get("cycles"),
                  "contended": r.get("contended"), "metrics": r.get("metrics")} for r in res.get("runs") or []],
        "artifacts": [a["path"] for a in res.get("artifacts") or []],
        "logs": res.get("logs"),
    }


def result_view(res: dict[str, Any], section: str) -> dict[str, Any]:
    """One section of a result, shared by the `isolab_result` MCP tool and `isolab result --section`."""
    if section == "full":
        return res
    if section == "summary":
        return summary_view(res)
    if section == "fidelity":
        return {"fidelity": res["fidelity"], "placement": res.get("placement"),
                "per_run_checks": [{"index": r["index"], "contended": r.get("contended"), "checks": r.get("checks"),
                                    "conditions": r.get("conditions")} for r in res.get("runs") or []]}
    if section in res:
        return {section: res[section]}
    raise ValueError(f"unknown section {section!r}")
