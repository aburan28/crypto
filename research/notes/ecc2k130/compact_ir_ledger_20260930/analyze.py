#!/usr/bin/env python3
"""Check archived Callgrind evidence and derive same-Q instruction rows."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import tarfile

from run_panel import CONFIG, sha


def archived_report(archive: Path, cell_id: str, entry: dict) -> dict | None:
    raw = archive / entry["raw_path"]
    assert entry["raw_path"] == f"raw/{cell_id}.tar.gz"
    assert raw.stat().st_size == entry["raw_bytes"] and sha(raw) == entry["raw_sha256"]
    hashes = entry["member_sha256"]
    if entry["run_json_sha256"] is not None:
        assert hashes["run.json"] == entry["run_json_sha256"]
    if entry["receipt_sha256"] is not None:
        assert hashes["receipt.json"] == entry["receipt_sha256"]
    report = None
    seen: set[str] = set()
    with tarfile.open(raw, "r:gz") as tar:
        for member in tar:
            assert member.isfile() and member.name.startswith(f"{cell_id}/")
            relative = member.name[len(cell_id) + 1:]
            assert relative in hashes and relative not in seen
            seen.add(relative)
            stream = tar.extractfile(member)
            assert stream is not None
            digest = hashlib.sha256()
            chunks = [] if relative == "run.json" else None
            while chunk := stream.read(1 << 20):
                digest.update(chunk)
                if chunks is not None:
                    chunks.append(chunk)
            assert digest.hexdigest() == hashes[relative], (cell_id, relative)
            if chunks is not None:
                report = json.loads(b"".join(chunks))
    assert seen == set(hashes)
    return report


def analyze(archive: Path) -> dict:
    manifest = json.loads((archive / "MANIFEST.json").read_text())
    assert manifest["schema"] == "ecc2k130-compact-ir-hosted-archive-v1"
    assert manifest["config_sha256"] == sha(CONFIG)
    if manifest["second_host_path"] is not None:
        host_path = archive / manifest["second_host_path"]
        assert sha(host_path) == manifest["second_host_sha256"]
    config = json.loads(CONFIG.read_text())
    assert set(manifest["cases"]) == {cell["id"] for cell in config["cells"]}
    result = {"schema": "ecc2k130-compact-ir-analysis-v1",
              "source_head": manifest["source_head"],
              "github_run_url": manifest["github_run_url"],
              "manifest_sha256": sha(archive / "MANIFEST.json"),
              "config_sha256": sha(CONFIG), "cells": {}}
    for cell in config["cells"]:
        cell_id = cell["id"]
        entry = manifest["cases"][cell_id]
        status = entry["status"]
        if status == "MISSING":
            result["cells"][cell_id] = {"status": "MISSING", "rows": []}
            continue
        report = archived_report(archive, cell_id, entry)
        if report is not None:
            assert report["host"]["git_head"] == manifest["source_head"]
            assert report["config_sha256"] == sha(CONFIG)
            assert report["cell"] == cell
        if status != "PASS":
            result["cells"][cell_id] = {
                "status": status, "rows": [],
                "partial_arms": ([] if report is None else [
                    {"arm": arm["arm"], "exit_code": arm["exit_code"],
                     "stopped_for": arm["stopped_for"], "Ir": arm["Ir"]}
                    for arm in report["runs"]])}
            continue
        receipt_path = archive / entry["receipt_path"]
        assert entry["receipt_path"] == f"receipts/{cell_id}.json"
        assert sha(receipt_path) == entry["receipt_sha256"]
        receipt = json.loads(receipt_path.read_text())
        assert receipt["status"] == "PASS" and receipt["cell"] == cell_id
        assert receipt["n"] == cell["n"] and receipt["L"] == cell["L"]
        assert receipt["k"] == cell["k"]
        assert report is not None and report["status"] == "PASS"
        second_host_replayed = "second_host_replay_path" in entry
        if second_host_replayed:
            replay_path = archive / entry["second_host_replay_path"]
            assert sha(replay_path) == entry["second_host_replay_sha256"]
        assert {arm["arm"] for arm in report["runs"]} == set(report["sequence"])
        runs = {arm["arm"]: arm for arm in report["runs"]}
        rho_ir = runs["rho"]["Ir"]
        assert isinstance(rho_ir, int) and rho_ir > 0
        sqrt_r = math.sqrt(receipt["subgroup_order"])
        rows = []
        for arm in cell["arms"]:
            raw = runs[arm]
            detail = receipt["details"][arm]
            ir = raw["Ir"]
            assert isinstance(ir, int) and ir > 0 and ir == detail["Ir"]
            assert raw["exit_code"] == 0 and raw["stopped_for"] is None
            s_ir = ir / (cell["L"] * sqrt_r)
            assert math.isclose(s_ir, receipt["S_Ir"][arm], rel_tol=1e-12)
            ratio = ir / rho_ir
            if arm != "rho":
                assert math.isclose(ratio,
                                    receipt["instruction_ratio_IC_to_rho"][arm],
                                    rel_tol=1e-12)
            rows.append({"arm": arm, "n": cell["n"], "L": cell["L"],
                         "K": cell["k"] if arm != "rho" else None,
                         "verified_logs": detail["targets_verified"],
                         "complete_Ir": ir, "S_Ir": s_ir,
                         "same_Q_Ir_ratio_to_rho": ratio,
                         "callgrind_peak_rss_bytes": raw["observed_peak_rss_bytes"],
                         "attempt_floor_ratio": detail.get("attempt_floor_ratio")})
        result["cells"][cell_id] = {
            "status": "PASS", "rows": rows,
            "second_host_replayed": second_host_replayed,
            "host_machine": report["host"]["machine"],
            "valgrind": report["host"]["valgrind"],
            "deterministic_control": report.get("deterministic_control")}
    result["all_hosted_cells_pass"] = all(
        cell["status"] == "PASS" for cell in result["cells"].values())
    result["all_cells_second_host_replayed"] = all(
        cell.get("second_host_replayed", False) for cell in result["cells"].values())
    result["method_crossover"] = None
    result["method_crossover_note"] = (
        "This instruction panel has no disjoint confirmation CPU gate or n131 transfer; "
        "Ir is not group-addition equivalent.")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite an analysis result"
    result = analyze(args.archive)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value["status"] for key, value in
                      result["cells"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
