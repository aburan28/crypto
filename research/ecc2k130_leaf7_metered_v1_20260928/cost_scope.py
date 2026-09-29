#!/usr/bin/env python3
"""Audit a future cost receipt without constructing or timing an isogeny.

The producer's field counters still need independent code review; a JSON flag
cannot prove instrumentation. This checker prevents incomplete scopes and
unbound structural outputs from being promoted by arithmetic alone.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_digest(path: Path) -> str:
    return digest(path.read_bytes())


def canonical(data: object) -> bytes:
    return (json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n").encode()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parent_output_digest(spec: dict, producer: dict) -> str:
    labels = spec["point_panel"]["labels"]
    cases = producer["cases"]
    require(len(cases) == len(labels), "parent certificate panel length differs")
    rows = []
    for label, case in zip(labels, cases):
        require(case["label"] == label, f"parent case label differs: {label}")
        direct = case["direct_leaf1"]
        bridge = case["sage_oriented_leaf1"]
        require(direct == bridge, f"parent bridge/direct point mismatch: {label}")
        require(isinstance(direct, list) and len(direct) == 2
                and all(isinstance(x, str) and x.startswith("0x") for x in direct),
                f"parent lacks a canonical full affine point: {label}")
        rows.append([label, direct])
    return digest(canonical(rows))


def phase_total(phases: dict, required: list[str], square_weight: float) -> float:
    require(set(phases) == set(required), "missing or extra metered phase")
    total = 0.0
    for name in required:
        row = phases[name]
        require(set(row) == {"mul", "sqr", "inv_calls"},
                f"wrong operation-counter fields: {name}")
        for key in ("mul", "sqr", "inv_calls"):
            require(type(row[key]) is int and row[key] >= 0,
                    f"nonnegative integer count required: {name}.{key}")
        total += row["mul"] + square_weight * row["sqr"]
    return total


def audit(
    spec: dict,
    receipt: dict,
    producer: dict,
    producer_sha: str,
    extension: dict,
    extension_sha: str,
    bridge_freeze_sha: str,
    own_freeze_sha: str,
) -> dict:
    require(spec["release_main_head"] is not None,
            "metered gate remains held; refreeze before a measured receipt")
    require(receipt["schema"] == "ecc2k130-leaf7-metered-receipt-v1",
            "wrong cost receipt schema")
    require(receipt["freeze_sha256"] == own_freeze_sha, "cost freeze hash differs")
    require(receipt["release_main_head"] == spec["release_main_head"],
            "cost release head differs")
    require(receipt["parent_producer_sha256"] == producer_sha,
            "wrong parent producer")
    require(receipt["parent_extension_replay_sha256"] == extension_sha,
            "wrong parent extension replay")
    require(producer["status"] == spec["prerequisite"]["producer_status"],
            "parent producer is not complete")
    require(producer["freeze_sha256"] == bridge_freeze_sha,
            "parent bridge freeze differs")
    require(extension["status"] == spec["prerequisite"]["extension_replay_status"]
            and extension["structural_gate_status"]
            == spec["prerequisite"]["structural_gate_status"],
            "parent structural gate is not PASS")
    require(extension["freeze_sha256"] == bridge_freeze_sha,
            "parent extension freeze differs")
    require(extension["source_receipt_sha256"] == producer_sha,
            "parent extension does not bind producer")
    expected = parent_output_digest(spec, producer)
    outputs = receipt["outputs"]
    require(outputs["case_labels"] == spec["point_panel"]["labels"],
            "cost output labels differ")
    require(outputs["direct_sha256"] == expected
            and outputs["bridge_sha256"] == expected,
            "cost output does not match independently checked parent points")

    if receipt["status"] == "STOP":
        return {"schema": "ecc2k130-leaf7-metered-audit-v1", "status": "STOP",
                "promotion": False, "incremental_ratio": None,
                "cold_ratio": None, "reason": "producer stopped or exceeded a cap"}
    require(receipt["status"] == "MEASURED", "invalid cost receipt status")
    coverage = receipt["coverage"]
    required_flags = (
        "all_phase_counts_exclusive",
        "torsion_discovery_metered",
        "sage_7_division_and_kernel_discovery_metered",
        "all_extension_arithmetic_reduced_to_Fq",
        "inversion_internals_in_base_counts",
    )
    require(isinstance(coverage.get("unmetered_steps"), list)
            and all(isinstance(step, str) for step in coverage["unmetered_steps"]),
            "unmetered_steps must be a list of phase names")
    incomplete = [flag for flag in required_flags if coverage.get(flag) is not True]
    incomplete += coverage["unmetered_steps"]
    if incomplete:
        return {"schema": "ecc2k130-leaf7-metered-audit-v1", "status": "HOLD",
                "promotion": False, "incremental_ratio": None,
                "cold_ratio": None, "reason": "unmetered scope: " + ", ".join(incomplete)}

    calibration = receipt["calibration"]
    frozen = spec["cost"]
    require(calibration["label"] == frozen["square_weight_label"],
            "calibration label differs")
    require(calibration["batches"] == frozen["calibration_batches"]
            and calibration["operands_per_batch"]
            == frozen["calibration_operands_per_batch"],
            "calibration sample size differs")
    require(calibration["same_host"] is True, "calibration used another host")
    weight = calibration["square_per_mul"]
    require(type(weight) in (float, int) and math.isfinite(weight) and weight > 0,
            "invalid square weight")
    mads = [calibration["multiply_relative_mad"], calibration["square_relative_mad"]]
    require(all(type(mad) in (float, int) and math.isfinite(mad) and mad >= 0
                for mad in mads), "invalid relative MAD")
    if any(mad > frozen["max_calibration_relative_mad"] for mad in mads):
        return {"schema": "ecc2k130-leaf7-metered-audit-v1", "status": "HOLD",
                "promotion": False, "incremental_ratio": None,
                "cold_ratio": None, "reason": "calibration MAD exceeds frozen limit"}

    resources = receipt["resources"]
    for block in ("common", "direct", "bridge"):
        row = resources[block]
        require(row["status"] in ("COMPLETE", "STOP"), "invalid child status")
        require(type(row["wall_seconds"]) in (float, int)
                and math.isfinite(row["wall_seconds"]) and row["wall_seconds"] >= 0,
                "invalid child wall time")
        require(type(row["peak_rss_bytes"]) is int and row["peak_rss_bytes"] >= 0,
                "invalid child peak RSS")
        if (row["status"] == "STOP"
                or row["wall_seconds"] > spec["caps"]["child_wall_seconds"]
                or row["peak_rss_bytes"] > spec["caps"]["child_peak_rss_bytes"]):
            return {"schema": "ecc2k130-leaf7-metered-audit-v1", "status": "STOP",
                    "promotion": False, "incremental_ratio": None,
                    "cold_ratio": None, "reason": f"{block} child stopped or exceeded cap"}

    counts = receipt["counts"]
    require(set(counts) == {"common", "direct", "bridge"}, "wrong cost blocks")
    common = phase_total(counts["common"], frozen["common_phases"], weight)
    direct = phase_total(counts["direct"], frozen["direct_phases"], weight)
    bridge = phase_total(counts["bridge"], frozen["bridge_phases"], weight)
    require(direct > 0 and common + direct > 0, "zero direct reference cost")
    incremental = bridge / direct
    cold = (common + bridge) / (common + direct)
    promotion = (
        incremental <= frozen["incremental_bridge_over_direct_at_most"]
        and cold < frozen["cold_bridge_over_direct_below"]
    )
    return {
        "schema": "ecc2k130-leaf7-metered-audit-v1",
        "status": "PASS",
        "promotion": promotion,
        "classification": ("cold_route_engineering" if promotion
                           else "stage_diagnostic_or_negative"),
        "unit": frozen["unit"],
        "common_mul_equivalent": common,
        "direct_incremental_mul_equivalent": direct,
        "bridge_incremental_mul_equivalent": bridge,
        "direct_cold_mul_equivalent": common + direct,
        "bridge_cold_mul_equivalent": common + bridge,
        "incremental_ratio": incremental,
        "cold_ratio": cold,
        "pdp_or_dlp_claim": False,
    }


def self_test(spec: dict) -> None:
    demo = copy.deepcopy(spec)
    demo["release_main_head"] = "a" * 40
    labels = demo["point_panel"]["labels"]
    cases = [{"label": label, "direct_leaf1": ["0x1", "0x2"],
              "sage_oriented_leaf1": ["0x1", "0x2"]} for label in labels]
    producer = {"status": "PRODUCER_PASS", "freeze_sha256": "f" * 64, "cases": cases}
    producer_sha = digest(canonical(producer))
    extension = {"status": "PASS", "structural_gate_status": "PASS",
                 "freeze_sha256": "f" * 64,
                 "source_receipt_sha256": producer_sha}
    extension_sha = digest(canonical(extension))
    output = parent_output_digest(demo, producer)
    frozen = demo["cost"]
    counts = {
        block: {name: {"mul": value, "sqr": 0, "inv_calls": 0}
                for name in frozen[block + "_phases"]}
        for block, value in (("common", 20), ("direct", 20), ("bridge", 10))
    }
    receipt = {
        "schema": "ecc2k130-leaf7-metered-receipt-v1",
        "status": "MEASURED",
        "freeze_sha256": "e" * 64,
        "release_main_head": demo["release_main_head"],
        "parent_producer_sha256": producer_sha,
        "parent_extension_replay_sha256": extension_sha,
        "outputs": {"case_labels": labels, "direct_sha256": output,
                    "bridge_sha256": output},
        "coverage": {
            "all_phase_counts_exclusive": True,
            "torsion_discovery_metered": True,
            "sage_7_division_and_kernel_discovery_metered": True,
            "all_extension_arithmetic_reduced_to_Fq": True,
            "inversion_internals_in_base_counts": True,
            "unmetered_steps": [],
        },
        "calibration": {
            "label": frozen["square_weight_label"],
            "batches": frozen["calibration_batches"],
            "operands_per_batch": frozen["calibration_operands_per_batch"],
            "same_host": True,
            "square_per_mul": 0.8,
            "multiply_relative_mad": 0.01,
            "square_relative_mad": 0.01,
        },
        "resources": {block: {"status": "COMPLETE", "wall_seconds": 1.0,
                              "peak_rss_bytes": 1000}
                      for block in ("common", "direct", "bridge")},
        "counts": counts,
    }
    args = (demo, receipt, producer, producer_sha, extension,
            extension_sha, "f" * 64, "e" * 64)
    assert audit(*args)["promotion"] is True
    bad = copy.deepcopy(receipt)
    del bad["counts"]["common"]["torsion_discovery"]
    try:
        audit(demo, bad, producer, producer_sha, extension, extension_sha,
              "f" * 64, "e" * 64)
    except ValueError:
        pass
    else:
        raise AssertionError("missing cold phase was accepted")
    bad = copy.deepcopy(receipt)
    bad["coverage"]["all_extension_arithmetic_reduced_to_Fq"] = False
    assert audit(demo, bad, producer, producer_sha, extension, extension_sha,
                 "f" * 64, "e" * 64)["cold_ratio"] is None
    bad = copy.deepcopy(receipt)
    bad["outputs"]["bridge_sha256"] = "0" * 64
    try:
        audit(demo, bad, producer, producer_sha, extension, extension_sha,
              "f" * 64, "e" * 64)
    except ValueError:
        pass
    else:
        raise AssertionError("mismatched point outputs were accepted")
    bad = copy.deepcopy(receipt)
    bad["calibration"]["square_relative_mad"] = 0.06
    assert audit(demo, bad, producer, producer_sha, extension, extension_sha,
                 "f" * 64, "e" * 64)["incremental_ratio"] is None
    bad = copy.deepcopy(receipt)
    bad["resources"]["bridge"]["peak_rss_bytes"] = 2147483649
    assert audit(demo, bad, producer, producer_sha, extension, extension_sha,
                 "f" * 64, "e" * 64)["status"] == "STOP"
    print("Cost-scope self-test PASS; no map or measurement was run.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--bridge-producer", type=Path)
    parser.add_argument("--bridge-extension-replay", type=Path)
    args = parser.parse_args()
    spec = json.loads((HERE / "FROZEN.json").read_text())
    if args.self_test:
        require(args.receipt is None, "self-test cannot audit a measured receipt")
        self_test(spec)
        return
    require(all([args.receipt, args.bridge_producer, args.bridge_extension_replay]),
            "receipt, bridge producer, and bridge extension replay are all required")
    parent_freeze = REPO / spec["prerequisite"]["structural_directory"] / "FROZEN.json"
    result = audit(
        spec,
        json.loads(args.receipt.read_text()),
        json.loads(args.bridge_producer.read_text()),
        file_digest(args.bridge_producer),
        json.loads(args.bridge_extension_replay.read_text()),
        file_digest(args.bridge_extension_replay),
        file_digest(parent_freeze),
        file_digest(HERE / "FROZEN.json"),
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
