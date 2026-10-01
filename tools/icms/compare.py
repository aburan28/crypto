"""Paired comparison of two arms of a session.

Two figures, with different admission rules (AGENTS.md sections 2, 6, 8, 10):

* **Operations** in the spec's unit, over the spec's window.  Admitted when
  both arms measured the same problem (instance, window, unit, reference) and
  differ only in declared variables.  A deterministic unit gives an exact
  ratio; a host-dependent one (a solver priced from wall time, a measured
  calibration ratio) is admitted only on the same env_class and is labelled.
* **Wall time** of the window: additionally both arms at or above their
  declared isolation level, the same env_class_id and pinned-CPU count, and
  paired by round.  Reported as the median per-round ratio with a percentile
  bootstrap 95% interval, beside the A/A spread when the session has one.

A comparison is never computed across sessions with different env classes,
and never across units: a refusal names the reason instead.
"""
from __future__ import annotations

import json
import os
import random
import statistics
from typing import Any

from . import COMPARISON_SCHEMA
from .spec import SpecError, identity_view, load as load_spec

# Identity paths that define the measured problem; never a declared variable.
FORBIDDEN_PREFIXES = ("instance.", "measurement.window", "accounting.", "reference.", "relations.stop",
                      "relations.count", "execution.threads", "execution.cpus")


def load_session(path: str) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    """The session, its records, and each spec copy loaded and normalised the
    way its spec id was computed (None for a copy that no longer loads)."""
    with open(os.path.join(path, "session.json")) as fh:
        session = json.load(fh)
    with open(os.path.join(path, "records.jsonl")) as fh:
        records = [json.loads(line) for line in fh if line.strip()]
    specs: dict[str, Any] = {}
    for name in session["spec_files"]:
        try:
            specs[name] = load_spec(os.path.join(path, "specs", name))
        except (SpecError, OSError, ValueError):
            specs[name] = None
    return session, records, specs


def _flatten(d: Any, prefix: str = "") -> dict[str, Any]:
    """Leaf paths of a nested mapping.  An empty mapping is a leaf, and a key
    that itself contains a dot is bracketed, so two different documents never
    flatten to the same map."""
    out: dict[str, Any] = {}
    if isinstance(d, dict) and d:
        for k, v in d.items():
            part = f'["{k}"]' if "." in k else k
            out.update(_flatten(v, (f"{prefix}{part}" if part.startswith("[") else f"{prefix}.{part}") if prefix else part))
    else:
        out[prefix] = d
    return out


def spec_differences(a: dict[str, Any], b: dict[str, Any], declared: list[str]) -> dict[str, list[dict[str, Any]]]:
    fa, fb = _flatten(identity_view(a)), _flatten(identity_view(b))
    out: dict[str, list[dict[str, Any]]] = {"declared": [], "confound": [], "forbidden": []}
    for k in sorted(set(fa) | set(fb)):
        if fa.get(k) == fb.get(k):
            continue
        item = {"field": k, "a": fa.get(k), "b": fb.get(k)}
        if k.startswith(FORBIDDEN_PREFIXES) or k == "instance.workload.seeds":
            out["forbidden"].append(item)
        elif any(k == d or k.startswith(d + ".") for d in declared):
            out["declared"].append(item)
        else:
            out["confound"].append(item)
    return out


def bootstrap_median_ratio(pairs: list[tuple[float, float]], iters: int = 20000, seed: int = 20261001) -> dict[str, Any]:
    ratios = [b / a for a, b in pairs if a > 0]
    if not ratios:
        return {"n_pairs": 0, "median_ratio": None, "ci95": None}
    if len(ratios) < 5:
        return {"n_pairs": len(ratios), "median_ratio": statistics.median(ratios), "ci95": None,
                "note": "fewer than 5 pairs: AGENTS.md section 10 asks for at least five ABAB rounds; no interval"}
    rng = random.Random(seed)
    meds = sorted(statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(iters))
    lo, hi = meds[int(0.025 * iters)], meds[int(0.975 * iters) - 1]
    return {"n_pairs": len(ratios), "median_ratio": statistics.median(ratios), "min_ratio": min(ratios),
            "max_ratio": max(ratios), "ci95": [lo, hi], "excludes_1": hi < 1 or lo > 1,
            "method": f"percentile bootstrap of the median per-round ratio, {iters} resamples, seed {seed}"}


def _arm_records(records: list[dict[str, Any]], arm: int) -> list[dict[str, Any]]:
    return [r for r in records if r["arm"] == arm and not r["warmup"]]


def _median(xs: list[float]) -> float | None:
    xs = [x for x in xs if x is not None]
    return statistics.median(xs) if xs else None


def structural(recs: list[dict[str, Any]]) -> dict[str, Any]:
    """Median of every numeric structural metric across an arm's runs."""
    acc: dict[str, list[float]] = {}
    for r in recs:
        for section in ("factor_base", "system", "solver", "decomposition", "linear_algebra"):
            for k, v in _flatten((r.get("metrics") or {}).get(section) or {}, section).items():
                if isinstance(v, (int, float)) and not isinstance(v, bool):
                    acc.setdefault(k, []).append(float(v))
    return {k: {"median": statistics.median(v), "min": min(v), "max": max(v), "n": len(v)} for k, v in sorted(acc.items())}


def compare(session_dir: str, arm_a: int, arm_b: int, declared: list[str]) -> dict[str, Any]:
    session, records, specs = load_session(session_dir)
    arms = {a["index"]: a for a in session["arms"]}
    A, B = arms[arm_a], arms[arm_b]
    la, lb = specs.get(os.path.basename(A["spec_path"])), specs.get(os.path.basename(B["spec_path"]))
    refusals = []
    for arm, loaded in ((A, la), (B, lb)):
        if loaded is None or loaded["spec_id"] != arm["spec_id"]:
            refusals.append({"reason": f"arm {arm['index']}'s spec copy does not load to its spec id {arm['spec_id']}"})
    spec_a = (la or {}).get("spec") or {}
    spec_b = (lb or {}).get("spec") or {}
    diffs = spec_differences(spec_a, spec_b, declared) if la and lb else {"declared": [], "confound": [], "forbidden": []}
    if la and lb and A["spec_id"] != B["spec_id"] and not any(diffs.values()):
        refusals.append({"reason": "the spec ids differ but no field difference was found; refusing rather than guessing"})
    ra, rb = _arm_records(records, arm_a), _arm_records(records, arm_b)
    if diffs["forbidden"]:
        refusals.append({"reason": "the arms measure different problems", "fields": [d["field"] for d in diffs["forbidden"]]})
    if diffs["confound"]:
        refusals.append({"reason": "undeclared differences (confounds); declare them or remove them",
                         "fields": [d["field"] for d in diffs["confound"]]})
    if A["workload_id"] != B["workload_id"]:
        refusals.append({"reason": "different workloads", "a": A["workload_id"], "b": B["workload_id"]})
    inconsistent = [r["record_id"] for r in ra + rb if r.get("consistency")]
    if inconsistent:
        refusals.append({"reason": "records where the producer's report contradicts the spec", "record_ids": inconsistent})
    out: dict[str, Any] = {
        "schema": COMPARISON_SCHEMA, "session_id": session["session_id"], "env_class_id": session["env_class_id"],
        "records_sha256": session["records_sha256"],
        "request": {"a": arm_a, "b": arm_b, "declared": list(declared)},
        "a": {"arm": arm_a, "spec_id": A["spec_id"], "label": A.get("label")},
        "b": {"arm": arm_b, "spec_id": B["spec_id"], "label": B.get("label")},
        "aa": A["spec_id"] == B["spec_id"],
        "declared_variables": diffs["declared"], "confounds": diffs["confound"], "forbidden": diffs["forbidden"],
        "refusals": refusals,
        "structure": {"a": structural(ra), "b": structural(rb)},
        "outcomes": {"a": [r["outcome"]["status"] for r in ra], "b": [r["outcome"]["status"] for r in rb]},
    }
    if refusals:
        out["ops"] = out["wall"] = {"admitted": False}
        return out

    unit = spec_a["accounting"]["unit"]
    window = spec_a["measurement"]["window"]

    def ops_of(r: dict[str, Any]) -> float | None:
        # The window's own operation count when it is in the spec's unit; else the
        # unit's total, which every adapter reports over the cold window only.
        w = (r.get("windows") or {}).get(window) or {}
        if w.get("ops") is not None and w.get("ops_unit") == unit:
            return float(w["ops"])
        if window != "cold_end_to_end":
            return None
        u = (r.get("units") or {}).get(unit) or {}
        return None if u.get("total") is None else float(u["total"])

    def det(r: dict[str, Any]) -> bool:
        return bool(((r.get("units") or {}).get(unit) or {}).get("deterministic"))

    def unpriced(r: dict[str, Any]) -> list[str]:
        return list(((r.get("units") or {}).get(unit) or {}).get("unpriced_counters") or [])

    ok_a = [r for r in ra if r["outcome"]["status"] == "complete" and r["outcome"].get("verified")]
    ok_b = [r for r in rb if r["outcome"]["status"] == "complete" and r["outcome"].get("verified")]
    va, vb = [ops_of(r) for r in ok_a], [ops_of(r) for r in ok_b]
    deterministic = all(det(r) for r in ok_a + ok_b)
    zero_priced = sorted({c for r in ok_a + ok_b for c in unpriced(r)})
    # Every measured run must complete and verify: a median over the survivors
    # of a budget or timeout is biased toward the cheap runs.
    all_done = bool(ra and rb) and len(ok_a) == len(ra) and len(ok_b) == len(rb)
    ops: dict[str, Any] = {"unit": unit, "window": window, "admitted": all_done and None not in va + vb,
                           "deterministic": deterministic, "lower_bound": bool(zero_priced),
                           "completed": {"a": f"{len(ok_a)}/{len(ra)}", "b": f"{len(ok_b)}/{len(rb)}"}}
    if zero_priced:
        ops["unpriced_counters"] = zero_priced
    if ops["admitted"]:
        ma, mb = statistics.median(va), statistics.median(vb)
        ops.update({"a_median": ma, "b_median": mb, "ratio_b_over_a": mb / ma if ma else None,
                    "a_spread": [min(va), max(va)], "b_spread": [min(vb), max(vb)],
                    "speedup_a_over_b": ma / mb if mb else None})
        if not deterministic:
            ops["label"] = ("host-dependent: some cost in this unit was priced from wall time on this host "
                            "(see units.<unit>.host_dependent_because); valid on this env_class only")
    else:
        ops["reason"] = ("not every measured run completed and verified, or the unit/window was not reported; "
                         "the total is unknown, not zero, and a median of the survivors would be biased")
    out["ops"] = ops

    wall_refusals = []
    low = [r["record_id"] for r in ra + rb if not r["isolation"]["wall_admissible"]]
    if low:
        levels = sorted({r["isolation"]["earned_level"] for r in ra + rb})
        wall_refusals.append({"reason": f"runs below their declared isolation level (earned {levels})",
                              "count": len(low), "blocking": sorted({b for r in ra + rb for b in r["isolation"]["blocking_next_level"]})})
    conf = {((r.get("windows") or {}).get(window) or {}).get("conformance") for r in ra + rb}
    if window == "online_one_target" and conf != {"exact"}:
        wall_refusals.append({"reason": "online window derived, not timed on the five exclusive online clocks; "
                                        "not admissible as a primary wall-clock speedup", "conformance": sorted(map(str, conf))})

    def wall_of(r: dict[str, Any]) -> float | None:
        w = (r.get("windows") or {}).get(window) or {}
        return float(w["wall_ns"]) if w.get("wall_ns") is not None else None

    by_round_a = {r["round"]: r for r in ok_a}
    by_round_b = {r["round"]: r for r in ok_b}
    pairs = [(wall_of(by_round_a[k]), wall_of(by_round_b[k])) for k in sorted(set(by_round_a) & set(by_round_b))]
    pairs = [(x, y) for x, y in pairs if x is not None and y is not None]
    wall = {"window": window, "pairs": len(pairs), "estimate": bootstrap_median_ratio(pairs),
            "refusals": wall_refusals, "admitted": not wall_refusals and len(pairs) >= 5,
            "contended_runs": sum(1 for r in ra + rb if r["isolation"]["earned_level"] in ("L0", "L1"))}
    if not wall["admitted"] and not wall_refusals:
        wall["refusals"].append({"reason": "fewer than five paired rounds (AGENTS.md section 10)"})
    out["wall"] = wall
    # The A/A spread of this session, if it has an arm repeated.
    seen: dict[str, int] = {}
    for a in session["arms"]:
        if a["spec_id"] in seen and a["workload_id"] == arms[seen[a["spec_id"]]]["workload_id"]:
            x, y = _arm_records(records, seen[a["spec_id"]]), _arm_records(records, a["index"])
            bx = {r["round"]: r for r in x if r["outcome"]["status"] == "complete"}
            by = {r["round"]: r for r in y if r["outcome"]["status"] == "complete"}
            aa_pairs = [(wall_of(bx[k]), wall_of(by[k])) for k in sorted(set(bx) & set(by))]
            aa_pairs = [(p, q) for p, q in aa_pairs if p is not None and q is not None]
            out.setdefault("aa_noise", []).append({"arms": [seen[a["spec_id"]], a["index"]], "spec_id": a["spec_id"],
                                                   "estimate": bootstrap_median_ratio(aa_pairs)})
        else:
            seen.setdefault(a["spec_id"], a["index"])
    return out
