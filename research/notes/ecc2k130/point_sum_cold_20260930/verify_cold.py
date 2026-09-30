#!/usr/bin/env python3
"""Independent full-rank, target and rho replay of one cold point-query cell."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930"))
from run_panel import sha  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field, verify as verify_rank  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_point_panel_20260929"))
from verify_panel import check_target  # noqa: E402
from run_cold import ARM_ORDER, LIMIT_BYTES, LIMIT_SECONDS, checked_spec, schedule  # noqa: E402
from verify_frozen import verify as verify_frozen  # noqa: E402

T_CRIT_95 = {5: 2.7764451051977987, 20: 2.093024054408263}


def rows(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def paired(values: list[float]) -> dict:
    assert len(values) in T_CRIT_95 and all(value > 0 for value in values)
    logs = [math.log(value) for value in values]
    half = T_CRIT_95[len(values)] * statistics.stdev(logs) / math.sqrt(len(values))
    middle = statistics.mean(logs)
    return {"median": statistics.median(values),
            "interval_95pct": [math.exp(middle - half), math.exp(middle + half)]}


def checked_files(run_dir: Path, item: dict) -> dict[str, Path]:
    assert all(".fixture.jsonl" not in value for value in
               item["command"] + list(item["environment"].values()))
    result = {}
    for name, metadata in item["files"].items():
        path = run_dir / metadata["name"]
        assert path.name == metadata["name"]
        assert path.stat().st_size == metadata["bytes"]
        assert sha(path) == metadata["sha256"]
        result[name] = path
    assert {"stdout", "stderr"} <= result.keys()
    return result


def verify(cell: str, run_dir: Path, relocated: bool = False) -> dict:
    input_receipt = verify_frozen()
    assert input_receipt["status"] == "PASS"
    spec = checked_spec(cell)
    report = json.loads((run_dir / "cold_run.json").read_text())
    assert report["schema"] == "ecc2k130-point-sum-cold-run-v1"
    assert report["cell"] == cell and report["spec"] == spec
    assert report["frozen_sha256"] == input_receipt["frozen_sha256"]
    assert report["limits"] == {"arm_timeout_seconds": LIMIT_SECONDS,
                                "arm_rss_bytes": LIMIT_BYTES,
                                "arm_address_space_bytes": LIMIT_BYTES}
    assert report["source"]["commit"] == "d8d67d3b683334db105195b0fc74d8e660dc4b5e"
    assert report["source"]["cargo_lock_sha256"] == sha(HERE / "CARGO_LOCK_SNAPSHOT.txt")
    assert report["source"]["compact_source_sha256"] == input_receipt_source("compact")
    assert report["source"]["rho_source_sha256"] == input_receipt_source("rho")
    assert len(report["source"]["compact_binary_sha256"]) == 64
    assert len(report["source"]["rho_binary_sha256"]) == 64
    assert "release:" in report["source"]["rustc_version_verbose"]
    assert report["source"]["cargo_version"].startswith("cargo 1.")
    assert report["run_child_source_sha256"] == sha(
        ROOT / "research/notes/ecc2k130/compact_ir_ledger_20260930/run_panel.py")
    assert report["runner_sha256"] == sha(HERE / "run_cold.py")
    mode = report["mode"]
    assert mode in ("smoke", "measure")
    expected = schedule(spec["blocks"], mode)
    assert report["plan"] == [{"block": b, "arm": arm} for b, arm in expected]
    completed = [(item["block"], item["arm"]) for item in report["runs"]]
    assert completed == expected[:len(completed)]
    for item in report["runs"]:
        checked_files(run_dir, item)
    if report["status"] != ("PASS" if mode == "measure" else "SMOKE_PASS"):
        assert report["status"] == "FAIL"
        return {"status": "PRODUCER_FAILURE", "cell": cell,
                "completed_children": len(completed), "timing_eligible": False,
                "cold_run_sha256": sha(run_dir / "cold_run.json")}
    assert completed == expected
    assert report["host"]["machine"] == ("x86_64" if mode == "measure"
                                               else report["host"]["machine"])
    if mode == "measure":
        assert report["host"]["reserved_cpu"] is not None
    labels = rows(HERE / spec["fixture_file"])
    curve = Curve(Field(spec["n"], spec["field_modulus_low_terms"]), 0)
    generator = tuple(spec["generator"])
    checked_labels: set[tuple[str, int]] = set()
    by_block: dict[int, dict[str, dict]] = {}
    checks = []
    for item in report["runs"]:
        paths = checked_files(run_dir, item)
        assert item["exit_code"] == 0 and item["stopped_for"] is None
        assert item["observed_peak_rss_bytes"] <= LIMIT_BYTES
        assert item["child_max_rss_kib_linux"] is None or (
            item["child_max_rss_kib_linux"] * 1024 <= LIMIT_BYTES)
        assert item["environment"]["RAYON_NUM_THREADS"] == "1"
        assert item["environment"]["LC_ALL"] == "C"
        if mode == "measure":
            assert item["command"][:3] == ["taskset", "-c",
                                            str(report["host"]["reserved_cpu"])]
        arm = item["arm"]
        block = item["block"]
        cpu = item["child_user_cpu_seconds"] + item["child_system_cpu_seconds"]
        assert cpu > 0
        check = {"block": block, "arm": arm, "cpu_seconds": cpu,
                 "wall_seconds": item["elapsed_under_backend_seconds_not_a_cost"],
                 "rss_kib": item["child_max_rss_kib_linux"]}
        if arm == "rho":
            env = item["environment"]
            assert env["KIC_RHO_DP_BITS"] == "4"
            assert env["KIC_RHO_CANON_BACKEND"] == "normal_basis"
            assert env["KIC_RHO_BATCH_CORPUS"] == spec["corpus"]
            assert Path(env["KIC_RHO_POINT_INPUT"]).name == Path(spec["points_file"]).name
            assert Path(item["command"][-6]).name == "koblitz_rho_batch_ks_v3"
            assert item["command"][-5:] == [str(spec["n"]), "0", "signed_frobenius",
                                            str(spec["L"]), str(spec["seed"])]
            data = rows(paths["stdout"])
            assert len(data) == spec["L"] + 1
            summary = data[-1]
            assert summary["kind"] == "rho_ks_batch_summary"
            assert summary["target_source"] == "public_point_jsonl"
            assert summary["quotient_mode"] == "signed_frobenius"
            assert summary["canonicalization_backend"] == "normal_basis"
            assert summary["parallel_walks"] == 32 and summary["dp_bits"] == 4
            assert summary["corpus"] == spec["corpus"] and summary["all_verified"]
            for index, (record, label) in enumerate(zip(data[:-1], labels)):
                assert record["kind"] == "rho_ks_batch_fixture"
                assert record["fixture_index"] == index
                assert record["published_fixture_scalar"] is None
                assert record["published_q"] == label["published_q"]
                scalar = record["recovered_fixture_scalar"]
                assert scalar == label["published_fixture_scalar"]
                assert curve.mul(scalar, generator) == tuple(record["published_q"])
            check["target_logs_verified"] = spec["L"]
            check["rho_walk_steps"] = summary["total_walk_steps"]
        else:
            env = item["environment"]
            assert env["KIC_S3_BATCH_WINDOW"] == "64"
            assert env["KIC_S3_PREFILTER"] == spec["prefilter"]
            assert env["KIC_QUERY_BACKEND"] == ("point_sum" if arm == "point_sum" else "s3")
            assert Path(item["command"][-5]).name == "koblitz_orbit_dlp_s3_batch"
            assert item["command"][-4] == f"construct:{spec['n']}:0:{spec['K']}"
            assert Path(item["command"][-3]).name == Path(spec["points_file"]).name
            if not relocated:
                assert item["command"][-3] == str((HERE / spec["points_file"]).resolve())
            assert item["command"][-2] == "7"
            assert Path(item["command"][-1]).name == paths["targets"].name
            assert Path(env["KIC_DUMP_BASE"]).name == paths["base"].name
            assert Path(env["KIC_DUMP_RANK"]).name == paths["rank"].name
            assert {"base", "rank", "targets"} <= paths.keys()
            rank = verify_rank(paths["rank"], paths["base"], paths["stdout"])
            assert rank["status"] == "PASS" and rank["rank"] == spec["K"]
            base, = rows(paths["base"])
            summary, = rows(paths["stdout"])
            targets = rows(paths["targets"])
            logs = rows(paths["rank"])[-1]["logs"]
            assert len(targets) == spec["L"]
            assert summary["query_backend"] == env["KIC_QUERY_BACKEND"]
            assert summary["root_prefilter_policy"] == (
                "off" if spec["prefilter"] == "off" else "blocked_bloom_512_3hash")
            assert summary["rank"] == spec["K"]
            assert summary["targets_solved"] == spec["L"] and summary["targets_failed"] == 0
            for record, label in zip(targets, labels):
                check_target(record, label, base, logs, curve, generator, checked_labels)
            check.update({"rank": rank, "target_logs_verified": spec["L"],
                          "base_hash": base["base_hash"],
                          "regular_states": summary["regular_states"],
                          "root_table_entries": summary["root_table_entries"],
                          "rank_attempts": summary["rank_attempts"],
                          "target_scalars": [row["recovered_scalar"] for row in targets],
                          "rank_logs": logs,
                          "witnesses": [row["point_indices"] for row in targets],
                          "point_index_y_bytes": summary["point_index_y_bytes"],
                          "point_index_group_additions": summary["point_index_group_additions"]})
        by_block.setdefault(block, {})[arm] = check
        checks.append(check)
    if mode == "smoke":
        return {"status": "SMOKE_PASS", "cell": cell,
                "checked": checks, "timing_eligible": False}
    assert all(set(by_block[block]) == set(ARM_ORDER) for block in range(spec["blocks"]))
    for block in range(spec["blocks"]):
        arms = by_block[block]
        a, candidate, b = (arms[name] for name in ("control_a", "point_sum", "control_b"))
        for key in ("base_hash", "regular_states", "root_table_entries", "rank_attempts",
                    "target_scalars", "rank_logs"):
            assert a[key] == b[key] == candidate[key], (block, key)
        assert a["point_index_y_bytes"] == b["point_index_y_bytes"] == 0
        assert candidate["point_index_y_bytes"] > 0
        assert candidate["point_index_group_additions"] == 2 * candidate["regular_states"]
    isolation_rows = rows(run_dir / "isolation.jsonl")
    assert len(isolation_rows) == 1
    isolation = isolation_rows[0]
    assert isolation["schema"] == "isolated-bench/1" and isolation["mode"] == "reserve"
    assert isolation["label"] == f"point-sum-cold-{cell}"
    assert report["host"]["reserved_cpu"] in isolation["reserved_cpus"]
    uncontended = (isolation["exit_status"] == 0
                   and isolation["contended_samples"] == 0
                   and not isolation["left_on_reserved"]["user_threads"])
    aa, point_over_control, point_over_rho = [], [], []
    for block in range(spec["blocks"]):
        arms = by_block[block]
        a = arms["control_a"]["cpu_seconds"]
        b = arms["control_b"]["cpu_seconds"]
        control = math.sqrt(a * b)
        aa.append(b / a)
        point_over_control.append(arms["point_sum"]["cpu_seconds"] / control)
        point_over_rho.append(arms["point_sum"]["cpu_seconds"] / arms["rho"]["cpu_seconds"])
    aa_result = paired(aa)
    aa_valid = 0.9 <= aa_result["median"] <= 1.1 and (
        aa_result["interval_95pct"][0] <= 1 <= aa_result["interval_95pct"][1])
    return {"status": "PASS", "cell": cell, "n": spec["n"], "L": spec["L"],
            "K": spec["K"], "blocks": spec["blocks"],
            "cold_run_sha256": sha(run_dir / "cold_run.json"),
            "frozen_sha256": input_receipt["frozen_sha256"],
            "checked": checks, "isolation": isolation,
            "uncontended": uncontended, "aa_valid": aa_valid,
            "timing_eligible": aa_valid and uncontended,
            "paired": {"control_b_over_a": aa_result,
                       "point_over_control": paired(point_over_control),
                       "point_over_rho": paired(point_over_rho)}}


def input_receipt_source(kind: str) -> str:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    return frozen[f"{kind}_source_sha256"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell", required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--relocated", action="store_true")
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a verification receipt"
    try:
        receipt = verify(args.cell, args.run_dir.resolve(), args.relocated)
    except BaseException as error:
        receipt = {"status": "FAIL", "error_type": type(error).__name__,
                   "error": str(error), "traceback": traceback.format_exc()}
        args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        raise
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in receipt.items() if key != "checked"},
                     sort_keys=True))


if __name__ == "__main__":
    main()
