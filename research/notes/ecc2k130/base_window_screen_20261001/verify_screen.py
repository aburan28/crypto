#!/usr/bin/env python3
"""Independently replay every rank, target and charged child in the window screen."""
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
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/disjoint_cold_v2_20261001"))
from check_pairing import check_target  # noqa: E402
from run_screen import (FROZEN, INPUT_RECEIPT, SOURCE_FREEZE_SHA256,
                        config, schedule)  # noqa: E402
from verify_generator import verify as verify_generator  # noqa: E402
from verify_inputs import verify as verify_inputs  # noqa: E402


def rows(path: Path) -> list:
    return [json.loads(line) for line in path.read_bytes().splitlines() if line.strip()]


def paired(values: list[float], critical: float) -> dict:
    assert len(values) == 5 and all(math.isfinite(value) and value > 0 for value in values)
    logs = [math.log(value) for value in values]
    center = statistics.mean(logs)
    half = critical * statistics.stdev(logs) / math.sqrt(5)
    return {
        "values": values,
        "median": statistics.median(values),
        "geometric_mean": math.exp(center),
        "interval_95pct": [math.exp(center - half), math.exp(center + half)],
    }


def isolation_ok(isolation: dict | None, cpu: int, cfg: dict) -> bool:
    if not isolation or isolation.get("schema") != "isolated-bench/1":
        return False
    pre = isolation.get("preflight", {})
    settle = pre.get("settle", {})
    conditions = pre.get("conditions", {})
    try:
        psi = [conditions[kind]["some"]["avg10"] for kind in ("psi_cpu", "psi_memory")]
    except (KeyError, TypeError):
        return False
    samples = isolation.get("samples", [])
    threshold = cfg["other_cpu_contention_threshold_seconds"]
    return bool(
        isolation.get("mode") == "reserve"
        and isolation.get("label") == "base-window-screen-n41"
        and cpu in isolation.get("reserved_cpus", [])
        and isolation.get("exit_status") == 0
        and "user_threads" in isolation.get("left_on_reserved", {})
        and settle.get("seconds") == 2.0
        and settle.get("other_cpu_seconds", float("inf")) <= 0.2
        and len(psi) == 2 and all(0.0 <= value <= 5.0 for value in psi)
        and samples
        and isolation.get("contended_samples") == 0
        and all(
            not sample.get("contended", True)
            and sample.get("other_cpu_seconds", float("inf")) <= threshold
            for sample in samples
        )
    )


def isolation_command_ok(isolation: dict | None, run_dir: Path, cpu: int,
                         build: dict, relocated: bool) -> bool:
    """Bind the reserve receipt to this run and these exact built binaries."""
    if not isolation:
        return False
    command = isolation.get("command")
    if not isinstance(command, list) or len(command) != 18:
        return False
    if Path(command[0]).name not in ("python", "python3"):
        return False
    script = Path(command[1])
    if script.name != "run_screen.py":
        return False
    if not relocated and script.resolve() != HERE / "run_screen.py":
        return False
    arguments = dict(zip(command[2::2], command[3::2]))
    if len(arguments) != 8 or set(arguments) != {
        "--source-root", "--generator", "--compact", "--rho",
        "--materialization", "--build-receipt", "--out", "--cpu",
    }:
        return False
    generator = Path(build["binaries"]["koblitz_base_window"]["path"])
    build_dir = generator.parents[3]
    expected = {
        "--source-root": build_dir / "source",
        "--generator": generator,
        "--compact": Path(build["binaries"]["koblitz_orbit_dlp_s3_batch"]["path"]),
        "--rho": Path(build["binaries"]["koblitz_rho_batch_ks_v3"]["path"]),
        "--materialization": build_dir / "materialization.json",
        "--build-receipt": build_dir / "BUILD_RECEIPT.json",
    }
    if any(arguments[key] != str(path) for key, path in expected.items()):
        return False
    if arguments["--cpu"] != str(cpu):
        return False
    output = Path(arguments["--out"])
    return output.name == run_dir.name if relocated else output == run_dir


def checked_files(run_dir: Path, child: dict) -> dict[str, Path]:
    assert child["command"][:2] == ["taskset", "-c"]
    assert all(".fixture.jsonl" not in word for word in child["command"])
    assert all(".fixture.jsonl" not in str(value) for value in child["environment"].values())
    result = {}
    for name, entry in child["files"].items():
        assert Path(entry["name"]).name == entry["name"]
        path = run_dir / entry["name"]
        assert path.is_file() and path.stat().st_size == entry["bytes"]
        assert sha(path) == entry["sha256"]
        result[name] = path
    assert {"stdout", "stderr"} <= result.keys()
    return result


def checked_cpu(child: dict, limit_bytes: int, limit_seconds: int) -> float:
    assert child["exit_code"] == 0 and child["stopped_for"] is None
    assert child["observed_peak_rss_bytes"] <= limit_bytes
    assert child["child_max_rss_kib_linux"] * 1024 <= limit_bytes
    assert child["elapsed_wall_seconds_not_primary_cost"] <= limit_seconds + 1
    assert 0 <= child["process_launch_wall_seconds_not_primary_cost"] <= (
        child["elapsed_wall_seconds_not_primary_cost"]
    )
    total = child["child_user_cpu_seconds"] + child["child_system_cpu_seconds"]
    assert math.isfinite(total) and total > 0
    return total


def verify(run_dir: Path, relocated: bool = False) -> dict:
    inputs = verify_inputs()
    assert inputs["status"] == "PASS"
    cfg = config()
    frozen = json.loads(FROZEN.read_text())
    report_path = run_dir / "screen_run.json"
    report = json.loads(report_path.read_text())
    assert report["schema"] == "ecc2k130-base-window-screen-run-v1"
    assert report["config"] == cfg
    assert report["frozen_sha256"] == inputs["frozen_sha256"] == sha(FROZEN)
    assert report["input_receipt_sha256"] == sha(INPUT_RECEIPT)
    assert json.loads(INPUT_RECEIPT.read_text()) == inputs
    assert report["runner_sha256"] == sha(HERE / "run_screen.py")
    assert report["source"]["source_freeze_sha256"] == SOURCE_FREEZE_SHA256
    assert report["source"]["generator_source_sha256"] == frozen["source_lock"]["sha256"][
        "examples/koblitz_base_window.rs"
    ]
    assert report["source"]["compact_source_sha256"] == cfg["compact_source_sha256"]
    assert report["source"]["rho_source_sha256"] == cfg["rho_source_sha256"]
    assert report["source"]["materialization_sha256"] == sha(run_dir / "materialization.json")
    assert report["source"]["build_receipt_sha256"] == sha(run_dir / "BUILD_RECEIPT.json")
    assert report["source"]["materialization"] == json.loads(
        (run_dir / "materialization.json").read_text()
    )
    assert report["source"]["materialization"]["pinned_files"] == 20
    build = json.loads((run_dir / "BUILD_RECEIPT.json").read_text())
    assert build["status"] == "PASS" and build["profile"] == "release"
    assert build["materialization_sha256"] == report["source"]["materialization_sha256"]
    assert build["generator_source_sha256"] == report["source"]["generator_source_sha256"]
    assert build["compact_source_sha256"] == report["source"]["compact_source_sha256"]
    assert build["rho_source_sha256"] == report["source"]["rho_source_sha256"]
    assert build["cargo_lock_sha256"] == report["source"]["cargo_lock_sha256"]
    assert build["fast_arith_sha256"] == report["source"]["frozen_fast_arith_sha256"]
    for name, key in (
        ("koblitz_base_window", "generator_binary_sha256"),
        ("koblitz_orbit_dlp_s3_batch", "compact_binary_sha256"),
        ("koblitz_rho_batch_ks_v3", "rho_binary_sha256"),
    ):
        assert build["binaries"][name]["sha256"] == report["source"][key]
    assert all(len(report["source"][key]) == 64 for key in (
        "generator_binary_sha256", "compact_binary_sha256", "rho_binary_sha256"
    ))
    assert report["host"]["machine"] == "x86_64"
    cpu = report["host"]["reserved_cpu"]
    assert isinstance(cpu, int) and cpu in report["host"]["affinity"]
    assert report["host"]["cpu_model"]
    assert report["limits"] == {
        "generator_wall_seconds": cfg["generator_wall_limit_seconds"],
        "generator_rss_bytes": cfg["generator_rss_limit_bytes"],
        "arm_wall_seconds": cfg["arm_wall_limit_seconds"],
        "arm_rss_bytes": cfg["arm_rss_limit_bytes"],
        "cell_wall_seconds": cfg["cell_wall_limit_seconds"],
    }
    expected = schedule(cfg)
    assert report["plan"] == [{"block": b, "arm": arm} for b, arm in expected]
    assert [(item["block"], item["arm"]) for item in report["runs"]] == expected[:len(report["runs"])]
    assert len(report["runs"]) <= len(expected) == 45

    curve = Curve(Field(cfg["n"], frozen["field_modulus_low_terms"]), cfg["a"])
    generator = tuple(frozen["generator"])
    assert curve.on_curve(generator) and curve.mul(cfg["subgroup_order"], generator) is None
    by_block: dict[int, dict[str, dict]] = {block: {} for block in range(cfg["blocks"])}
    window_headers: dict[int, tuple[str, list[int]]] = {}
    checked_labels: set[tuple[str, int]] = set()
    checks = []
    failed = None
    for item in report["runs"]:
        block, arm = item["block"], item["arm"]
        block_spec = frozen["blocks"][block]
        assert item["points_file"] == block_spec["points_file"]
        assert item["points_sha256"] == block_spec["points_sha256"]
        assert arm in cfg["arm_order_before_rotation"]
        children = item["children"]
        assert set(children) <= ({"rho"} if arm == "rho" else {"generator", "compact"})
        stages = ("rho",) if arm == "rho" else ("generator", "compact")
        for stage in stages:
            if stage not in children:
                continue
            child = children[stage]
            paths = checked_files(run_dir, child)
            assert child["command"][:3] == ["taskset", "-c", str(cpu)]
            executable = {
                "generator": "koblitz_base_window",
                "compact": "koblitz_orbit_dlp_s3_batch",
                "rho": "koblitz_rho_batch_ks_v3",
            }[stage]
            assert child["command"][3] == build["binaries"][executable]["path"]
            assert child["environment"]["RAYON_NUM_THREADS"] == "1"
            assert child["environment"]["LC_ALL"] == "C"
            if child["exit_code"] != 0 or child["stopped_for"] is not None:
                failed = {
                    "block": block, "arm": arm, "stage": stage,
                    "exit_code": child["exit_code"], "stopped_for": child["stopped_for"],
                }
                break
            if stage == "generator":
                assert len(child["command"]) == 12
                assert set(child["environment"]) == {"LC_ALL", "RAYON_NUM_THREADS"}
                assert {"header", "scan", "receipt"} <= paths.keys()
                window = int(arm[1])
                assert child["command"][-8:-3] == [
                    str(cfg["n"]), str(cfg["a"]), str(cfg["useful_orbit_columns"]),
                    str(window), str(cfg["raw_x_trial_cap"]),
                ]
                assert [Path(word).name for word in child["command"][-3:]] == [
                    paths[key].name for key in ("header", "scan", "receipt")
                ]
                if not relocated:
                    assert [Path(word) for word in child["command"][-3:]] == [
                        paths[key] for key in ("header", "scan", "receipt")
                    ]
                generator_check = verify_generator(
                    paths["header"], paths["scan"], paths["receipt"], protocol=True
                )
                assert generator_check["status"] == "PASS"
                gen_receipt = json.loads(paths["receipt"].read_text())
                assert len(gen_receipt["selected_orbit_keys"]) == cfg["useful_orbit_columns"]
                header = (generator_check["header_sha256"],
                          gen_receipt["selected_orbit_keys"])
                if window in window_headers:
                    assert window_headers[window] == header
                else:
                    window_headers[window] = header
                gen_cpu = checked_cpu(
                    child, cfg["generator_rss_limit_bytes"],
                    cfg["generator_wall_limit_seconds"],
                )
            elif stage == "compact":
                assert len(child["command"]) == 8
                assert set(child["environment"]) == {
                    "LC_ALL", "RAYON_NUM_THREADS", "KIC_DUMP_BASE", "KIC_DUMP_RANK",
                    "KIC_S3_BATCH_WINDOW", "KIC_S3_PREFILTER",
                }
                assert "generator" in children
                assert {"base", "rank", "targets"} <= paths.keys()
                assert Path(child["command"][3]).name == "koblitz_orbit_dlp_s3_batch"
                assert Path(child["command"][4]).name == checked_files(
                    run_dir, children["generator"]
                )["header"].name
                assert Path(child["command"][5]).name == Path(block_spec["points_file"]).name
                assert child["command"][6] == str(cfg["rank_seed"])
                assert Path(child["command"][7]).name == paths["targets"].name
                if not relocated:
                    assert Path(child["command"][4]) == checked_files(
                        run_dir, children["generator"]
                    )["header"]
                    assert Path(child["command"][5]) == (
                        HERE / block_spec["points_file"]
                    ).resolve()
                    assert Path(child["command"][7]) == paths["targets"]
                env = child["environment"]
                assert env["KIC_S3_BATCH_WINDOW"] == str(cfg["compact_batch_window"])
                assert env["KIC_S3_PREFILTER"] == cfg["compact_prefilter"]
                assert "KIC_QUERY_BACKEND" not in env
                assert Path(env["KIC_DUMP_BASE"]).name == paths["base"].name
                assert Path(env["KIC_DUMP_RANK"]).name == paths["rank"].name
                if not relocated:
                    assert Path(env["KIC_DUMP_BASE"]) == paths["base"]
                    assert Path(env["KIC_DUMP_RANK"]) == paths["rank"]
                header_path = checked_files(run_dir, children["generator"])["header"]
                assert paths["base"].read_bytes() == header_path.read_bytes()
                rank = verify_rank(paths["rank"], paths["base"], paths["stdout"])
                assert rank["status"] == "PASS" and rank["rank"] == cfg["useful_orbit_columns"]
                base, = rows(paths["base"])
                summary, = rows(paths["stdout"])
                target_rows = rows(paths["targets"])
                fixture_rows = rows(HERE / block_spec["fixture_file"])
                assert len(target_rows) == len(fixture_rows) == cfg["public_targets_per_block"]
                assert summary["rank"] == cfg["useful_orbit_columns"]
                assert summary["base_source"] == "retained_header"
                assert summary["base_hash"] == base["base_hash"]
                assert summary["targets_solved"] == len(target_rows)
                assert summary["targets_failed"] == 0
                assert summary["root_prefilter_policy"] == "blocked_bloom_512_3hash"
                logs = rows(paths["rank"])[-1]["logs"]
                recovered = []
                for index, (target, fixture) in enumerate(zip(target_rows, fixture_rows)):
                    assert fixture["fixture_index"] == index
                    label = {
                        "fixture_index": index,
                        "n": cfg["n"], "a": cfg["a"],
                        "published_q": fixture["published_q"],
                        "published_fixture_scalar": fixture["published_fixture_scalar"],
                    }
                    check_target(target, label, base, logs, curve, generator, checked_labels)
                    recovered.append(target["recovered_scalar"])
                compact_cpu = checked_cpu(
                    child, cfg["arm_rss_limit_bytes"], cfg["arm_wall_limit_seconds"]
                )
                by_block[block][arm] = {
                    "cpu_seconds": gen_cpu + compact_cpu,
                    "generator_cpu_seconds": gen_cpu,
                    "compact_cpu_seconds": compact_cpu,
                    "header_sha256": sha(header_path),
                    "base_hash": base["base_hash"],
                    "rank_logs": logs,
                    "rank_attempts": summary["rank_attempts"],
                    "target_scalars": recovered,
                    "rank_replay": rank,
                }
                checks.append({
                    "block": block, "arm": arm,
                    "cpu_seconds": gen_cpu + compact_cpu,
                    "generator_cpu_seconds": gen_cpu,
                    "compact_cpu_seconds": compact_cpu,
                    "rank": rank["rank"],
                    "target_logs_verified": len(recovered),
                    "base_hash": base["base_hash"],
                })
            else:
                assert stage == "rho" and arm == "rho"
                assert len(child["command"]) == 9
                assert set(child["environment"]) == {
                    "LC_ALL", "RAYON_NUM_THREADS", "KIC_RHO_POINT_INPUT",
                    "KIC_RHO_BATCH_CORPUS", "KIC_RHO_DP_BITS", "KIC_RHO_CANON_BACKEND",
                }
                assert Path(child["command"][3]).name == "koblitz_rho_batch_ks_v3"
                assert child["command"][-5:] == [
                    str(cfg["n"]), str(cfg["a"]), "signed_frobenius",
                    str(cfg["public_targets_per_block"]), str(block_spec["rho_seed"]),
                ]
                env = child["environment"]
                assert Path(env["KIC_RHO_POINT_INPUT"]).name == Path(
                    block_spec["points_file"]
                ).name
                if not relocated:
                    assert Path(env["KIC_RHO_POINT_INPUT"]) == (
                        HERE / block_spec["points_file"]
                    ).resolve()
                assert env["KIC_RHO_BATCH_CORPUS"] == block_spec["corpus"]
                assert env["KIC_RHO_DP_BITS"] == str(cfg["rho_dp_bits"])
                assert env["KIC_RHO_CANON_BACKEND"] == cfg["rho_canonical_backend"]
                rho_rows = rows(paths["stdout"])
                assert len(rho_rows) == cfg["public_targets_per_block"] + 1
                summary = rho_rows[-1]
                assert summary["kind"] == "rho_ks_batch_summary"
                assert summary["all_verified"] is True
                assert summary["quotient_mode"] == "signed_frobenius"
                assert summary["canonicalization_backend"] == cfg["rho_canonical_backend"]
                assert summary["parallel_walks"] == cfg["rho_walks"]
                assert summary["dp_bits"] == cfg["rho_dp_bits"]
                assert summary["corpus"] == block_spec["corpus"]
                fixture_rows = rows(HERE / block_spec["fixture_file"])
                recovered = []
                for index, (record, fixture) in enumerate(zip(rho_rows[:-1], fixture_rows)):
                    assert record["kind"] == "rho_ks_batch_fixture"
                    assert record["fixture_index"] == index
                    assert record["published_fixture_scalar"] is None
                    assert record["published_q"] == fixture["published_q"]
                    scalar = record["recovered_fixture_scalar"]
                    assert scalar == fixture["published_fixture_scalar"]
                    assert curve.mul(scalar, generator) == tuple(record["published_q"])
                    recovered.append(scalar)
                rho_cpu = checked_cpu(
                    child, cfg["arm_rss_limit_bytes"], cfg["arm_wall_limit_seconds"]
                )
                by_block[block][arm] = {
                    "cpu_seconds": rho_cpu,
                    "target_scalars": recovered,
                    "walk_steps": summary["total_walk_steps"],
                }
                checks.append({
                    "block": block, "arm": arm,
                    "cpu_seconds": rho_cpu,
                    "target_logs_verified": len(recovered),
                    "walk_steps": summary["total_walk_steps"],
                })
        if failed:
            break
    if report["status"] != "PASS":
        assert report["status"] in ("CENSORED", "FAIL")
        return {
            "schema": "ecc2k130-base-window-screen-replay-v1",
            "status": "PRODUCER_FAILURE",
            "timing_eligible": False,
            "failure": failed or report.get("failure"),
            "verified_children": len(checks),
            "checks": checks,
            "screen_run_sha256": sha(report_path),
        }
    assert failed is None and len(report["runs"]) == len(expected)
    assert set(window_headers) == set(range(4))
    for first in range(4):
        for second in range(first + 1, 4):
            assert not (
                set(window_headers[first][1]) & set(window_headers[second][1])
            )
    for block, arms in by_block.items():
        assert set(arms) == set(cfg["arm_order_before_rotation"])
        rho_scalars = arms["rho"]["target_scalars"]
        for window in range(4):
            a, b = arms[f"w{window}_a"], arms[f"w{window}_b"]
            for key in ("header_sha256", "base_hash", "rank_logs",
                        "rank_attempts", "target_scalars"):
                assert a[key] == b[key], (block, window, key)
            assert a["target_scalars"] == rho_scalars
    assert len(checks) == 45
    assert sum(check["target_logs_verified"] for check in checks if check["arm"] != "rho") == 40960
    assert sum(check["target_logs_verified"] for check in checks if check["arm"] == "rho") == 5120

    isolation_path = run_dir / "isolation.jsonl"
    isolation_rows = rows(isolation_path) if isolation_path.is_file() else []
    assert len(isolation_rows) <= 1
    isolation = isolation_rows[0] if isolation_rows else None
    uncontended = isolation_ok(isolation, cpu, cfg) and isolation_command_ok(
        isolation, run_dir, cpu, build, relocated
    )
    if isolation:
        uncontended = uncontended and isolation.get("host", {}).get("machine") == "x86_64"
        uncontended = uncontended and isolation.get("host", {}).get("cpu_model") == (
            report["host"]["cpu_model"]
        )
    aa = {}
    fixed = {}
    for window in range(4):
        aa_values, ratios = [], []
        for block in range(cfg["blocks"]):
            arms = by_block[block]
            a = arms[f"w{window}_a"]["cpu_seconds"]
            b = arms[f"w{window}_b"]["cpu_seconds"]
            rho_cpu = arms["rho"]["cpu_seconds"]
            aa_values.append(b / a)
            ratios.append(math.sqrt(a * b) / rho_cpu)
        aa[str(window)] = paired(aa_values, cfg["log_t_critical_df4"])
        fixed[str(window)] = paired(ratios, cfg["log_t_critical_df4"])
    aa_valid = {
        key: cfg["aa_ratio_min"] <= stats["median"] <= cfg["aa_ratio_max"]
        and stats["interval_95pct"][0] <= 1 <= stats["interval_95pct"][1]
        for key, stats in aa.items()
    }
    return {
        "schema": "ecc2k130-base-window-screen-replay-v1",
        "status": "PASS",
        "timing_eligible": uncontended and all(aa_valid.values()),
        "uncontended": uncontended,
        "aa_valid": aa_valid,
        "aa": aa,
        "fixed_window_over_rho": fixed,
        "window_header_sha256": {str(w): window_headers[w][0] for w in range(4)},
        "independent_rank_replays": 40,
        "verified_compact_target_logs": 40960,
        "verified_rho_target_logs": 5120,
        "screen_run_sha256": sha(report_path),
        "isolation_sha256": sha(isolation_path) if isolation else None,
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--relocated", action="store_true")
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite an independent replay"
    try:
        result = verify(args.run_dir.resolve(), args.relocated)
    except BaseException as error:
        result = {
            "schema": "ecc2k130-base-window-screen-replay-v1",
            "status": "FAIL",
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": traceback.format_exc(),
        }
        args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        raise
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        key: value for key, value in result.items() if key not in ("checks", "aa",
                                                                    "fixed_window_over_rho")
    }, sort_keys=True))


if __name__ == "__main__":
    main()
