#!/usr/bin/env python3
"""Independent Python group-law replay of a whole-process instruction cell."""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))
from run_panel import CONFIG, load_cell, sha  # noqa: E402

POINT_VERIFIER = (ROOT / "research/notes/ecc2k130/compact_orbit_point_panel_20260929"
                  / "verify_panel.py")
point_spec = importlib.util.spec_from_file_location("compact_point_verifier", POINT_VERIFIER)
assert point_spec is not None and point_spec.loader is not None
point_module = importlib.util.module_from_spec(point_spec)
point_spec.loader.exec_module(point_module)
check_fixture, check_target, rows = (point_module.check_fixture,
                                     point_module.check_target, point_module.rows)
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import verify as verify_rank  # noqa: E402

TARGET_IDENTITY_KEYS = ("target", "published_q", "recovered_scalar", "group_verified",
                        "point_indices", "x_codes", "pinned_intermediates", "probes")


def replay_ir(path: Path) -> int:
    """Parse Ir independently of the runner and reject absent totals."""
    events: tuple[str, ...] = ()
    value: int | None = None
    with path.open(errors="replace") as stream:
        for line in stream:
            if line.startswith("events:"):
                events = tuple(line.split()[1:])
                assert events.count("Ir") == 1
            elif line.startswith(("summary:", "totals:")):
                assert events
                counts = tuple(int(token) for token in line.split()[1:])
                assert len(counts) == len(events)
                value = counts[events.index("Ir")]
    assert value is not None and value > 0
    return value


def identity(record: dict) -> tuple[str, ...]:
    return tuple(json.dumps(record[key], sort_keys=True) for key in TARGET_IDENTITY_KEYS)


def verify(cell_id: str, run_dir: Path, smoke: bool = False) -> dict:
    config, cell, spec, points, fixture = load_cell(cell_id)
    report = json.loads((run_dir / "run.json").read_text())
    assert report["schema"] == "ecc2k130-compact-ir-run-v1"
    assert report["cell"] == cell and report["spec"] == spec
    assert report["config_sha256"] == sha(CONFIG)
    assert report["source_freeze_sha256"] == config["source_freeze_sha256"]
    assert report["points_sha256"] == cell["points_sha256"]
    assert report["materialization"]["schema"] == "compact-frozen-source-materialization-v1"
    source = json.loads((ROOT / config["source_freeze"]).read_text())
    assert report["materialization"]["pinned_files"] == len(source["source_sha256"])
    assert len(report["binaries"]["compact_sha256"]) == 64
    assert len(report["binaries"]["rho_sha256"]) == 64
    assert report["host"]["machine"] in ("x86_64", "aarch64", "arm64")
    assert report["backend"] == ("native" if smoke else "callgrind")
    assert report["status"] == ("SMOKE_PASS" if smoke else "PASS")
    expected = list(cell["arms"])
    expected.extend(f"{arm}_repeat" for arm in cell.get("repeat_control", []))
    if smoke:
        assert (report["sequence"] == expected or
                (len(report["sequence"]) == 1 and report["sequence"][0] in expected))
        expected = list(report["sequence"])
    assert report["sequence"] == expected
    assert [item["arm"] for item in report["runs"]] == expected
    local_spec = dict(spec, points_file=str(points), fixture_file=str(fixture))
    known, curve, generator = check_fixture(local_spec)
    assert len(known) == cell["L"]
    checked_labels: set[tuple[str, int]] = set()
    details = {}
    compact_evidence = {}
    for item in report["runs"]:
        arm, policy = item["arm"], item["policy"]
        assert policy == arm.removesuffix("_repeat")
        assert item["exit_code"] == 0 and item["stopped_for"] is None
        assert item["observed_peak_rss_bytes"] <= config["rss_limit_bytes_per_arm"]
        assert 0 < item["elapsed_under_backend_seconds_not_a_cost"] <= (
            config["timeout_seconds_per_arm"] + 2)
        assert item["environment"]["RAYON_NUM_THREADS"] == "1"
        assert not any(".fixture.jsonl" in value for value in
                       item["command"] + list(item["environment"].values()))
        stdout, stderr = run_dir / item["stdout"], run_dir / item["stderr"]
        assert sha(stdout) == item["stdout_sha256"]
        assert sha(stderr) == item["stderr_sha256"]
        if smoke:
            assert item["Ir"] is None
        else:
            profile = run_dir / item["callgrind"]
            assert sha(profile) == item["callgrind_sha256"]
            assert profile.stat().st_size == item["callgrind_bytes"]
            assert replay_ir(profile) == item["Ir"] > 0
            assert item["command"][0] == "valgrind"
            assert report["host"]["valgrind"].startswith("valgrind-")
        if policy in ("off", "blocked"):
            assert item["environment"]["KIC_S3_BATCH_WINDOW"] == "64"
            assert item["environment"]["KIC_S3_PREFILTER"] == policy
            assert item["command"][-4:-1] == [f"construct:{cell['n']}:0:{cell['k']}",
                                               str(points), "7"]
            base_path = run_dir / f"{arm}.base.jsonl"
            rank_path = run_dir / f"{arm}.rank.jsonl"
            target_path = run_dir / f"{arm}.target.jsonl"
            assert Path(item["environment"]["KIC_DUMP_BASE"]) == base_path
            assert Path(item["environment"]["KIC_DUMP_RANK"]) == rank_path
            assert Path(item["command"][-1]) == target_path
            rank = verify_rank(rank_path, base_path, stdout)
            assert rank["status"] == "PASS" and rank["rank"] == cell["k"]
            base, = rows(base_path)
            rank_rows = rows(rank_path)
            logs = rank_rows[-1]["logs"]
            summary, = rows(stdout)
            assert summary["kind"] == "compact_orbit_dlp_summary"
            assert summary["n"] == cell["n"] and summary["orbit_columns"] == cell["k"]
            assert summary["s3_batch_window"] == 64
            assert summary["root_prefilter_policy"] == (
                "blocked_bloom_512_3hash" if policy == "blocked" else "off")
            assert summary["rank"] == cell["k"]
            assert summary["targets"] == summary["targets_solved"] == cell["L"]
            assert summary["targets_failed"] == 0 and summary["base_hash"] == base["base_hash"]
            target_rows = rows(target_path)
            assert len(target_rows) == cell["L"]
            for target, label in zip(target_rows, known):
                check_target(target, label, base, logs, curve, generator, checked_labels)
            compact_evidence[arm] = (base, rank_rows,
                                     tuple(identity(target) for target in target_rows), summary)
            details[arm] = {"Ir": item["Ir"], "rank": rank,
                            "targets_verified": len(target_rows),
                            "rank_attempts": summary["rank_attempts"],
                            "attempt_floor_ratio": (summary["rank_attempts"] + cell["L"])
                                                   / (cell["k"] + cell["L"]),
                            "index_s3_counts": summary["index_s3_counts"],
                            "rank_s3_counts": summary["rank_s3_counts"],
                            "target_s3_counts": summary["target_s3_counts"]}
        else:
            assert policy == "rho"
            assert item["environment"]["KIC_RHO_CANON_BACKEND"] == "normal_basis"
            assert item["environment"]["KIC_RHO_DP_BITS"] == "4"
            assert item["environment"]["KIC_RHO_BATCH_CORPUS"] == spec["corpus"]
            assert item["environment"]["KIC_RHO_POINT_INPUT"] == str(points)
            assert item["command"][-5:] == [str(cell["n"]), "0", "signed_frobenius",
                                               str(cell["L"]), str(spec["seed"])]
            data = rows(stdout)
            assert len(data) == cell["L"] + 1
            summary = data[-1]
            assert summary["kind"] == "rho_ks_batch_summary"
            assert summary["n"] == cell["n"] and summary["fixtures"] == cell["L"]
            assert summary["all_verified"] and summary["target_source"] == "public_point_jsonl"
            assert summary["quotient_mode"] == "signed_frobenius"
            assert summary["canonicalization_backend"] == "normal_basis"
            assert summary["parallel_walks"] == 32 and summary["dp_bits"] == 4
            assert summary["corpus"] == spec["corpus"]
            assert summary["inversion_backend"] == "itoh_tsujii"
            for index, (record, label) in enumerate(zip(data[:-1], known)):
                assert record["kind"] == "rho_ks_batch_fixture"
                assert record["fixture_index"] == index
                assert record["published_fixture_scalar"] is None
                assert record["target_source"] == "public_point_jsonl"
                assert record["published_q"] == label["published_q"]
                scalar = record["recovered_fixture_scalar"]
                assert scalar == label["published_fixture_scalar"]
                assert curve.mul(scalar, generator) == tuple(record["published_q"])
            details[arm] = {"Ir": item["Ir"], "targets_verified": cell["L"],
                            "walk_steps": summary["total_walk_steps"],
                            "rho_charges": summary["charges"]}
    if "blocked" in compact_evidence:
        assert compact_evidence["blocked"][:3] == compact_evidence["off"][:3]
        before, after = compact_evidence["off"][3], compact_evidence["blocked"][3]
        for key in ("rank", "rank_attempts", "rank_relations", "rank_failures",
                    "rank_rows_without_gain", "base_hash", "regular_states",
                    "root_table_entries", "root_table_slots"):
            assert before[key] == after[key], key
    if "off_repeat" in compact_evidence:
        assert compact_evidence["off_repeat"][:3] == compact_evidence["off"][:3]
    if "rho_repeat" in details:
        assert details["rho_repeat"]["targets_verified"] == details["rho"]["targets_verified"]
    if smoke:
        return {"status": "SMOKE_PASS", "cell": cell_id, "details": details,
                "S_Ir": None, "instruction_ratio": None}
    if cell_id == "n37_L1":
        control = report["deterministic_control"]
        assert set(control) == {"off", "rho"} and all(
            0 <= value <= 0.001 for value in control.values())
    sqrt_r = math.sqrt(spec["subgroup_order"])
    s_ir = {arm: details[arm]["Ir"] / (cell["L"] * sqrt_r)
            for arm in cell["arms"]}
    ratio = {arm: details[arm]["Ir"] / details["rho"]["Ir"]
             for arm in cell["arms"] if arm != "rho"}
    return {"status": "PASS", "cell": cell_id, "n": cell["n"], "L": cell["L"],
            "k": cell["k"], "subgroup_order": spec["subgroup_order"],
            "points_sha256": sha(points), "config_sha256": sha(CONFIG),
            "S_Ir": s_ir, "instruction_ratio_IC_to_rho": ratio,
            "floor_ratio_group_equivalent": None, "details": details}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell", required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a verification receipt"
    try:
        receipt = verify(args.cell, args.run_dir.resolve(), args.smoke)
    except BaseException as error:
        receipt = {"status": "FAIL", "error_type": type(error).__name__,
                   "error": str(error), "traceback": traceback.format_exc()}
        args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        raise
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in receipt.items() if key != "details"},
                     sort_keys=True))


if __name__ == "__main__":
    main()
