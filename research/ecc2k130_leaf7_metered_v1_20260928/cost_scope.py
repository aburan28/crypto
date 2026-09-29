#!/usr/bin/env python3
"""Audit an archived route-cost receipt; never construct or time an isogeny.

This checks exact scope and provenance. A declared operation counter is not
proof of instrumentation: the future metering producer needs separate review.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
HEX64 = re.compile(r"[0-9a-f]{64}\Z")
HEX40 = re.compile(r"[0-9a-f]{40}\Z")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_digest(path: Path) -> str:
    return digest(path.read_bytes())


def canonical(data: object) -> bytes:
    return (json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n").encode()


def held_result(status: str, reason: str) -> dict:
    return {
        "schema": "ecc2k130-leaf7-metered-audit-v1",
        "status": status,
        "promotion": False,
        "incremental_ratio": None,
        "cold_ratio": None,
        "reason": reason,
    }


def full_point(value: object, label: str) -> list[str]:
    require(isinstance(value, list) and len(value) == 2,
            f"{label}: full affine point (x,y) required")
    q = 1 << 131
    for component in value:
        require(isinstance(component, str) and re.fullmatch(r"0x[0-9a-f]+", component),
                f"{label}: noncanonical coordinate text")
        require(0 <= int(component, 16) < q and hex(int(component, 16)) == component,
                f"{label}: coordinate outside canonical Fq range")
    return value


def parent_output_rows(spec: dict, producer: dict) -> list[list[object]]:
    labels = spec["point_panel"]["labels"]
    cases = producer["cases"]
    require(isinstance(cases, list) and len(cases) == len(labels),
            "parent certificate panel length differs")
    rows = []
    for label, case in zip(labels, cases):
        require(case["label"] == label, f"parent case label differs: {label}")
        direct = full_point(case["direct_leaf1"], f"parent {label} direct")
        bridge = full_point(case["sage_oriented_leaf1"], f"parent {label} bridge")
        require(direct == bridge, f"parent bridge/direct point mismatch: {label}")
        rows.append([label, direct])
    return rows


def archived_output_rows(spec: dict, value: object, arm: str) -> list[list[object]]:
    labels = spec["point_panel"]["labels"]
    require(isinstance(value, list) and len(value) == len(labels),
            f"{arm}: archived full-point panel missing or wrong length")
    rows = []
    for label, saved in zip(labels, value):
        require(isinstance(saved, dict) and set(saved) == {"label", "point"},
                f"{arm}: labelled full affine row required")
        require(saved["label"] == label, f"{arm}: label/order mismatch at {label}")
        rows.append([label, full_point(saved["point"], f"{arm} {label}")])
    return rows


def phase_total(phases: dict, required: list[str], square_weight: float) -> float:
    require(isinstance(phases, dict) and set(phases) == set(required),
            "missing or extra metered phase")
    total = 0.0
    for name in required:
        row = phases[name]
        require(isinstance(row, dict)
                and set(row) == {"mul", "sqr", "inv_calls"},
                f"wrong operation-counter fields: {name}")
        for key in ("mul", "sqr", "inv_calls"):
            require(type(row[key]) is int and row[key] >= 0,
                    f"nonnegative integer count required: {name}.{key}")
        total += row["mul"] + square_weight * row["sqr"]
    return total


def raw_calibration(spec: dict, calibration: dict, host_sha: str) -> tuple[float, float, float]:
    frozen = spec["cost"]
    require(calibration["label"] == frozen["square_weight_label"],
            "calibration label differs")
    require(calibration["host_sha256"] == host_sha, "calibration host differs")
    raw = calibration["raw_batches"]
    require(isinstance(raw, list)
            and len(raw) == frozen["calibration_batches"] == 11,
            "eleven raw calibration batches required")
    n = frozen["calibration_operands_per_batch"]
    mul_samples = []
    sq_samples = []
    for i, row in enumerate(raw):
        require(isinstance(row, dict)
                and set(row) == {
                    "batch", "operand_pair_sha256", "multiply_operations",
                    "square_operations", "multiply_seconds", "square_seconds",
                    "host_sha256",
                }, f"raw calibration batch {i} has wrong fields")
        require(row["batch"] == i
                and row["operand_pair_sha256"]
                == frozen["calibration_operand_pair_sha256"][i],
                f"raw calibration batch {i} uses wrong operands")
        require(row["multiply_operations"] == n
                and row["square_operations"] == n,
                f"raw calibration batch {i} has wrong operation count")
        require(row["host_sha256"] == host_sha,
                f"raw calibration batch {i} used another host")
        for name, samples in (("multiply_seconds", mul_samples),
                              ("square_seconds", sq_samples)):
            value = row[name]
            require(type(value) in (float, int)
                    and math.isfinite(value) and value > 0,
                    f"raw calibration batch {i} has invalid {name}")
            samples.append(value / n)
    mul_median = statistics.median(mul_samples)
    sq_median = statistics.median(sq_samples)
    mul_mad = statistics.median(abs(x - mul_median) for x in mul_samples) / mul_median
    sq_mad = statistics.median(abs(x - sq_median) for x in sq_samples) / sq_median
    weight = sq_median / mul_median
    for name, value in (
        ("square_per_mul", weight),
        ("multiply_relative_mad", mul_mad),
        ("square_relative_mad", sq_mad),
    ):
        reported = calibration[name]
        require(type(reported) in (float, int)
                and math.isfinite(reported)
                and math.isclose(reported, value, rel_tol=1e-9, abs_tol=1e-12),
                f"reported calibration {name} differs from raw batches")
    return weight, mul_mad, sq_mad


def audit(
    spec: dict,
    receipt: dict,
    own_freeze_sha: str,
    *,
    producer: dict | None = None,
    producer_sha: str | None = None,
    fq: dict | None = None,
    fq_sha: str | None = None,
    extension: dict | None = None,
    extension_sha: str | None = None,
    bridge_freeze_sha: str | None = None,
    metering_producer_sha: str | None = None,
) -> dict:
    require(receipt.get("schema") == "ecc2k130-leaf7-metered-receipt-v1",
            "wrong cost receipt schema")
    require(receipt.get("freeze_sha256") == own_freeze_sha,
            "cost freeze hash differs")
    require(receipt.get("release_main_head") == spec["release_main_head"],
            "cost release head differs")
    status = receipt.get("status")
    if status == "STOP":
        # A cap, crash, or setup failure need not have any target outputs.
        # Preserve it without requiring a structural producer or cost arrays.
        return held_result("STOP", str(receipt.get("error", "producer STOP")))
    require(status == "MEASURED", "invalid cost receipt status")
    require(spec["release_main_head"] is not None,
            "metered gate remains held; refreeze before a measured receipt")
    metering = spec["metering_producer"]
    require(metering["status"] == "REVIEWED"
            and isinstance(metering["sha256"], str)
            and HEX64.fullmatch(metering["sha256"]) is not None
            and isinstance(metering["commit"], str)
            and HEX40.fullmatch(metering["commit"]) is not None,
            "reviewed metering producer and frozen source commit required")
    require(metering_producer_sha == metering["sha256"]
            and receipt.get("metering_producer_sha256") == metering["sha256"]
            and receipt.get("metering_producer_commit") == metering["commit"],
            "metering producer source differs")
    require(receipt.get("input_sha256") == spec["input_sha256"],
            "cost receipt input/source hashes differ")
    require(all(x is not None for x in (
        producer, producer_sha, fq, fq_sha, extension, extension_sha,
        bridge_freeze_sha,
    )), "producer, Fq replay, and extension replay are required")
    require(receipt.get("parent_producer_sha256") == producer_sha,
            "wrong parent producer")
    require(receipt.get("parent_fq_replay_sha256") == fq_sha,
            "wrong parent Fq replay")
    require(receipt.get("parent_extension_replay_sha256") == extension_sha,
            "wrong parent extension replay")
    prerequisite = spec["prerequisite"]
    require(producer["status"] == prerequisite["producer_status"]
            and producer["phase"] == "complete",
            "parent producer is not complete")
    require(producer["freeze_sha256"] == bridge_freeze_sha,
            "parent bridge freeze differs")
    require(fq["status"] == prerequisite["fq_replay_status"]
            and fq["freeze_sha256"] == bridge_freeze_sha
            and fq["source_receipt_sha256"] == producer_sha
            and fq["independent_Fq_case_labels"] == spec["point_panel"]["labels"],
            "parent Fq replay is not bound and PASS on all twelve cases")
    require(extension["status"] == prerequisite["extension_replay_status"]
            and extension["structural_gate_status"]
            == prerequisite["structural_gate_status"]
            and extension["freeze_sha256"] == bridge_freeze_sha
            and extension["source_receipt_sha256"] == producer_sha
            and extension["fq_replay_sha256"] == fq_sha,
            "parent extension replay is not bound and PASS")

    expected_rows = parent_output_rows(spec, producer)
    expected_sha = digest(canonical(expected_rows))
    outputs = receipt["outputs"]
    require(isinstance(outputs, dict)
            and set(outputs) == {
                "direct", "bridge", "direct_sha256", "bridge_sha256",
            }, "archived labelled direct and bridge point arrays required")
    for arm in ("direct", "bridge"):
        rows = archived_output_rows(spec, outputs[arm], arm)
        actual_sha = digest(canonical(rows))
        require(actual_sha == outputs[arm + "_sha256"],
                f"{arm}: claimed output hash differs from archived full points")
        require(actual_sha == expected_sha and rows == expected_rows,
                f"{arm}: archived full points differ from parent certificate")

    host = receipt["host_manifest"]
    required_host = spec["cost"]["host_manifest_fields"]
    require(isinstance(host, dict) and set(host) == set(required_host),
            "host manifest fields differ")
    for name in required_host:
        if name == "field_primitive_sha256":
            continue
        require(isinstance(host[name], str) and host[name],
                f"host manifest {name} missing")
    require(host["git_head"] == metering["commit"],
            "host checkout differs from metering producer commit")
    primitives = {
        path: spec["input_sha256"][path]
        for path in spec["cost"]["field_primitive_paths"]
    }
    require(host["field_primitive_sha256"] == primitives,
            "host field primitive hashes differ")
    host_sha = digest(canonical(host))
    require(receipt["host_sha256"] == host_sha,
            "host manifest digest differs")

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
        return held_result("HOLD", "unmetered scope: " + ", ".join(incomplete))

    calibration = receipt["calibration"]
    weight, mul_mad, sq_mad = raw_calibration(spec, calibration, host_sha)
    if max(mul_mad, sq_mad) > spec["cost"]["max_calibration_relative_mad"]:
        return held_result("HOLD", "raw calibration MAD exceeds frozen limit")

    resources = receipt["resources"]
    require(isinstance(resources, dict)
            and set(resources) == {"common", "direct", "bridge"},
            "common/direct/bridge resource records required")
    for block in ("common", "direct", "bridge"):
        row = resources[block]
        require(row["status"] in ("COMPLETE", "STOP"), "invalid child status")
        require(row["host_sha256"] == host_sha, "child used another host")
        for clock_name in ("wall_seconds", "cpu_seconds"):
            require(type(row[clock_name]) in (float, int)
                    and math.isfinite(row[clock_name]) and row[clock_name] >= 0,
                    f"invalid child {clock_name}")
        require(type(row["peak_rss_bytes"]) is int and row["peak_rss_bytes"] >= 0,
                "invalid child peak RSS")
        if (row["status"] == "STOP"
                or row["wall_seconds"] > spec["caps"]["child_wall_seconds"]
                or row["peak_rss_bytes"] > spec["caps"]["child_peak_rss_bytes"]):
            return held_result("STOP", f"{block} child stopped or exceeded cap")

    counts = receipt["counts"]
    require(isinstance(counts, dict)
            and set(counts) == {"common", "direct", "bridge"},
            "wrong cost blocks")
    frozen = spec["cost"]
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
        "host_sha256": host_sha,
        "parent_fq_replay_sha256": fq_sha,
        "parent_extension_replay_sha256": extension_sha,
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
    # Synthetic data exercises gate logic only; it is not a measured result.
    stop = {
        "schema": "ecc2k130-leaf7-metered-receipt-v1",
        "freeze_sha256": "e" * 64,
        "release_main_head": spec["release_main_head"],
        "status": "STOP",
        "error": "synthetic cap",
    }
    assert audit(spec, stop, "e" * 64)["status"] == "STOP"
    assert audit(spec, stop, "e" * 64)["cold_ratio"] is None

    demo = copy.deepcopy(spec)
    demo["release_main_head"] = "a" * 40
    demo["metering_producer"].update(
        status="REVIEWED", sha256="d" * 64, commit="c" * 40)
    labels = demo["point_panel"]["labels"]
    point = ["0x1", "0x2"]
    cases = [{"label": label, "direct_leaf1": point,
              "sage_oriented_leaf1": point} for label in labels]
    producer = {"status": "PRODUCER_PASS", "phase": "complete",
                "freeze_sha256": "f" * 64, "cases": cases}
    producer_sha = digest(canonical(producer))
    fq = {"status": "FQ_REPLAY_PASS", "freeze_sha256": "f" * 64,
          "source_receipt_sha256": producer_sha,
          "independent_Fq_case_labels": labels}
    fq_sha = digest(canonical(fq))
    extension = {"status": "PASS", "structural_gate_status": "PASS",
                 "freeze_sha256": "f" * 64,
                 "source_receipt_sha256": producer_sha,
                 "fq_replay_sha256": fq_sha}
    extension_sha = digest(canonical(extension))
    rows = [[label, point] for label in labels]
    output_sha = digest(canonical(rows))
    frozen = demo["cost"]
    primitives = {p: demo["input_sha256"][p] for p in frozen["field_primitive_paths"]}
    host = {"os": "synthetic", "arch": "synthetic", "cpu_model": "synthetic",
            "python_version": "synthetic", "sage_version": "synthetic",
            "git_head": "c" * 40, "field_primitive_sha256": primitives}
    host_sha = digest(canonical(host))
    counts = {
        block: {name: {"mul": value, "sqr": 0, "inv_calls": 0}
                for name in frozen[block + "_phases"]}
        for block, value in (("common", 20), ("direct", 20), ("bridge", 10))
    }
    raw = [
        {"batch": i, "operand_pair_sha256": pair_sha,
         "multiply_operations": 10000, "square_operations": 10000,
         "multiply_seconds": 1.0, "square_seconds": 0.8,
         "host_sha256": host_sha}
        for i, pair_sha in enumerate(frozen["calibration_operand_pair_sha256"])
    ]
    receipt = {
        "schema": "ecc2k130-leaf7-metered-receipt-v1",
        "status": "MEASURED",
        "freeze_sha256": "e" * 64,
        "release_main_head": demo["release_main_head"],
        "metering_producer_sha256": "d" * 64,
        "metering_producer_commit": "c" * 40,
        "input_sha256": demo["input_sha256"],
        "parent_producer_sha256": producer_sha,
        "parent_fq_replay_sha256": fq_sha,
        "parent_extension_replay_sha256": extension_sha,
        "outputs": {
            "direct": [{"label": label, "point": point} for label in labels],
            "bridge": [{"label": label, "point": point} for label in labels],
            "direct_sha256": output_sha,
            "bridge_sha256": output_sha,
        },
        "host_manifest": host,
        "host_sha256": host_sha,
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
            "host_sha256": host_sha,
            "raw_batches": raw,
            "square_per_mul": 0.8,
            "multiply_relative_mad": 0.0,
            "square_relative_mad": 0.0,
        },
        "resources": {
            block: {"status": "COMPLETE", "wall_seconds": 1.0,
                    "cpu_seconds": 0.8, "peak_rss_bytes": 1000,
                    "host_sha256": host_sha}
            for block in ("common", "direct", "bridge")
        },
        "counts": counts,
    }
    parents = {
        "producer": producer, "producer_sha": producer_sha,
        "fq": fq, "fq_sha": fq_sha,
        "extension": extension, "extension_sha": extension_sha,
        "bridge_freeze_sha": "f" * 64,
        "metering_producer_sha": "d" * 64,
    }

    def check(row: dict) -> dict:
        return audit(demo, row, "e" * 64, **parents)

    assert check(receipt)["promotion"] is True
    bad = copy.deepcopy(receipt)
    del bad["outputs"]
    try:
        check(bad)
    except KeyError:
        pass
    else:
        raise AssertionError("MEASURED without archived point arrays was accepted")
    bad = copy.deepcopy(receipt)
    bad["outputs"]["bridge"][0]["point"] = ["0x3", "0x2"]
    try:
        check(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("tampered full point with stale digest was accepted")
    bad = copy.deepcopy(receipt)
    bad["outputs"]["bridge"][0]["point"] = ["0x03", "0x2"]
    try:
        check(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("noncanonical full point was accepted")
    bad_parents = dict(parents)
    bad_parents["fq"] = dict(fq, status="STOP")
    try:
        audit(demo, receipt, "e" * 64, **bad_parents)
    except ValueError:
        pass
    else:
        raise AssertionError("failed Fq replay was accepted")
    bad = copy.deepcopy(receipt)
    bad["parent_fq_replay_sha256"] = "0" * 64
    try:
        check(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("wrong Fq replay digest was accepted")
    bad = copy.deepcopy(receipt)
    del bad["counts"]["common"]["torsion_discovery"]
    try:
        check(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("missing cold phase was accepted")
    bad = copy.deepcopy(receipt)
    del bad["counts"]["bridge"]["first_leaf_panel_transport"]
    try:
        check(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("unpriced phi_C panel transport was accepted")
    bad = copy.deepcopy(receipt)
    bad["coverage"]["all_extension_arithmetic_reduced_to_Fq"] = False
    assert check(bad)["cold_ratio"] is None
    bad = copy.deepcopy(receipt)
    for i in range(6):
        bad["calibration"]["raw_batches"][i]["square_seconds"] = 1.2
    try:
        check(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("forged calibration summary was accepted")
    bad["calibration"]["square_per_mul"] = 1.2
    bad["calibration"]["square_relative_mad"] = 0.0
    assert math.isclose(raw_calibration(demo, bad["calibration"], host_sha)[0], 1.2)
    bad = copy.deepcopy(receipt)
    for i, batch in enumerate(bad["calibration"]["raw_batches"]):
        batch["multiply_seconds"] = 0.8 + 0.05 * i
    bad["calibration"]["square_per_mul"] = 0.8 / 1.05
    bad["calibration"]["multiply_relative_mad"] = 0.15 / 1.05
    assert check(bad)["status"] == "HOLD"
    bad = copy.deepcopy(receipt)
    bad["resources"]["bridge"]["peak_rss_bytes"] = 2147483649
    assert check(bad)["status"] == "STOP"
    bad = copy.deepcopy(receipt)
    bad["resources"]["bridge"]["host_sha256"] = "0" * 64
    try:
        check(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("cross-host bridge resource was accepted")
    print("Cost-scope self-test PASS; no map or measurement was run.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--bridge-producer", type=Path)
    parser.add_argument("--bridge-fq-replay", type=Path)
    parser.add_argument("--bridge-extension-replay", type=Path)
    args = parser.parse_args()
    spec = json.loads((HERE / "FROZEN.json").read_text())
    if args.self_test:
        require(args.receipt is None, "self-test cannot audit a measured receipt")
        self_test(spec)
        return
    require(args.receipt is not None, "cost receipt required")
    receipt = json.loads(args.receipt.read_text())
    own_freeze = file_digest(HERE / "FROZEN.json")
    if receipt.get("status") == "STOP":
        print(json.dumps(audit(spec, receipt, own_freeze), sort_keys=True))
        return
    require(all([args.bridge_producer, args.bridge_fq_replay,
                 args.bridge_extension_replay]),
            "MEASURED requires producer, Fq replay, and extension replay")
    for relative, expected in spec["input_sha256"].items():
        require(file_digest(REPO / relative) == expected,
                f"source hash changed: {relative}")
    parent_freeze = REPO / spec["prerequisite"]["structural_directory"] / "FROZEN.json"
    metering_path = REPO / spec["metering_producer"]["path"]
    result = audit(
        spec, receipt, own_freeze,
        producer=json.loads(args.bridge_producer.read_text()),
        producer_sha=file_digest(args.bridge_producer),
        fq=json.loads(args.bridge_fq_replay.read_text()),
        fq_sha=file_digest(args.bridge_fq_replay),
        extension=json.loads(args.bridge_extension_replay.read_text()),
        extension_sha=file_digest(args.bridge_extension_replay),
        bridge_freeze_sha=file_digest(parent_freeze),
        metering_producer_sha=file_digest(metering_path),
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
