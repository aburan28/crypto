#!/usr/bin/env python3
"""Summarize completed, hash-verified runs without changing their artifacts."""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent

def read(path):
    return json.loads(path.read_text())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def validate_run_names(names):
    # Repeating a passing primary is not independent confirmation.
    assert names in (["run_01"], ["run_01", "run_02"]), "Require one primary and, optionally, its distinct confirmation."
    if len(names) == 1:
        return
    primary, confirmation = [HERE / name for name in names]
    first = read(primary / "protocol.json")
    second = read(confirmation / "protocol.json")
    assert second == read(HERE / "confirmation_protocol.json")
    assert second["primary_manifest_sha256"] == sha(primary / "manifest.json")
    assert second["primary_protocol_sha256"] == sha(primary / "protocol.json")
    assert second["confirmation_plan_sha256"] == sha(HERE / "confirmation_plan.json")
    assert read(primary / "metadata.json")["source_hashes"] == read(confirmation / "metadata.json")["source_hashes"]
    for key in ["variables", "families", "variants", "repetitions", "limits", "new_candidates", "reference_arms", "matched_controls", "final_gate", "checksum_contract"]:
        assert first[key] == second[key], key
    assert second["regression_seeds"] == first["holdout_seeds"]
    seen = set(first["discovery_seeds"] + first["regression_seeds"] + first["holdout_seeds"])
    assert not seen.intersection(second["holdout_seeds"])

def summarize(name):
    root = HERE / name
    manifest = read(root / "manifest.json")
    for path, digest in manifest["files"].items():
        assert sha(root / path) == digest, (name, path)
    metadata, protocol, result = [read(root / file) for file in ["metadata.json", "protocol.json", "results.json"]]
    assert metadata["complete"]
    groups = defaultdict(list)
    pairs = defaultdict(dict)
    fixtures = {}
    for n in protocol["variables"]:
        for line in (root / f"raw-n{n}.jsonl").read_text().splitlines():
            row = json.loads(line)
            if row["type"] == "fixture":
                fixtures[row["cell"]] = row
                continue
            fixture = fixtures[row["cell"]]
            groups[fixture["split"], n, fixture["family"], row["variant"]].append(row)
            pairs[row["cell"]][row["rep"], row["variant"]] = row
    costs = []
    for n in protocol["variables"]:
        for arm in protocol["variants"]:
            family_costs = {}
            for family in protocol["families"]:
                values = groups["holdout", n, family, arm]
                complete = all(v["verified"] and v["outcome"] != "UNKNOWN" for v in values)
                family_costs[family] = {
                    "completion_ns": statistics.median(v["total_ns"] for v in values) if complete else None,
                    "observed_ns": statistics.median(v["total_ns"] for v in values),
                    "nodes": statistics.median(v["nodes"] for v in values),
                    "enumeration_points": statistics.median(v["enumeration_points"] for v in values),
                }
            costs.append({"n": n, "variant": arm, "families": family_costs})
    for cost in costs:
        reference = next(c for c in costs if c["n"] == cost["n"] and c["variant"] == "packed_untraced")
        before = reference["families"]["planted"]["completion_ns"]
        after = cost["families"]["planted"]["completion_ns"]
        cost["display_reference"] = "packed_untraced"
        cost["display_planted_ratio_of_medians"] = before / after if before is not None and after is not None else None
    case_ratios = []
    for cell, samples in sorted(pairs.items()):
        fixture = fixtures[cell]
        for candidate in protocol["new_candidates"]:
            controls = [a for a in protocol["reference_arms"] if a != candidate]
            ratios = []
            for rep in range(protocol["repetitions"]):
                values = [samples[rep, a] for a in controls + [candidate]]
                if all(v["verified"] and v["outcome"] != "UNKNOWN" for v in values):
                    ratios.append(min(samples[rep, a]["total_ns"] for a in controls) / samples[rep, candidate]["total_ns"])
            case_ratios.append({
                "cell": cell, "n": fixture["n"], "family": fixture["family"], "split": fixture["split"],
                "seed": fixture["seed"], "candidate": candidate,
                "paired_ratio_median": statistics.median(ratios) if len(ratios) == protocol["repetitions"] else None,
            })
    return {
        "name": name, "manifest_sha256": sha(root / "manifest.json"), "results_sha256": sha(root / "results.json"),
        "cells": result["cells"], "observations": result["samples"],
        "all_complete_verified": result["all_results_completed_and_verified"],
        "gates": result["new_gates"], "comparisons": result["new_details"],
        "costs": costs, "case_ratios": case_ratios,
        "outcomes_per_cell": dict(Counter(c["arms"]["gray_simd"]["outcome"] for c in result["cells_detail"])),
        "peak_worker_rss_bytes": max(c["peak_rss_bytes"] for c in result["cells_detail"]),
        "campaign_seconds": metadata["campaign_seconds"],
        "fixture_fingerprints": sorted(hashlib.sha256(json.dumps([f["n"], f["polys"]], separators=(",", ":")).encode()).hexdigest() for f in fixtures.values()),
    }

def main():
    names = sys.argv[1:] or ["run_01", "run_02"]
    validate_run_names(names)
    runs = [summarize(name) for name in names]
    candidates = read(HERE / names[0] / "protocol.json")["new_candidates"]
    accepted = {c: all(r["all_complete_verified"] and r["gates"][c]["retained_frontier"] == "PASS" for r in runs) for c in candidates}
    value = {
        "schema_version": 1, "classification": "engineering", "runs": runs,
        "combined_gates": {
            c: ("PASS" if len(runs) >= 2 else "PENDING_CONFIRMATION") if okay else "REJECTED"
            for c, okay in accepted.items()
        },
        "run_cells": sum(r["cells"] for r in runs), "observations": sum(r["observations"] for r in runs),
        "distinct_systems": len(set(f for r in runs for f in r["fixture_fingerprints"])),
        "production_solver_cost": None, "full_ic_cost": None, "rho_ratio": None, "calibrated_operation_ratio": None,
        "scope": "Complete bounded generated Boolean solves. Paired timing intervals concern fixed fixtures; no production, generic-exponent or cryptanalytic crossover claim.",
        "confirmation_guard_sha256": sha(HERE / "CONFIRMATION_GUARD.json"),
    }
    (HERE / "SUMMARY.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    ledger = {
        "schema_version": 1, "scope": value["scope"], "classification": "engineering",
        "source_lineage_sha256": sha(HERE / "SOURCE_LINEAGE.json"),
        "summary_sha256": sha(HERE / "SUMMARY.json"), "summarizer_sha256": sha(Path(__file__)),
        "artifacts": [{"path": r["name"], "manifest_sha256": r["manifest_sha256"], "kind": "complete-solve comparison"} for r in runs],
        "combined_gates": value["combined_gates"],
        "confirmation_guard_sha256": value["confirmation_guard_sha256"],
        "production_solver_cost": None, "full_ic_cost": None, "rho_ratio": None,
    }
    for path in sorted(HERE.glob("*resource_probe_*/manifest.json")):
        ledger["artifacts"].append({"path": path.parent.name, "manifest_sha256": sha(path), "kind": "discovery only"})
    (HERE / "RUN_LEDGER.json").write_text(json.dumps(ledger, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: value[k] for k in ["combined_gates", "run_cells", "observations", "distinct_systems"]}))

if __name__ == "__main__":
    main()
