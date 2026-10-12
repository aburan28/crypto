#!/usr/bin/env python3
"""Independently replay four cold arms on the frozen public-point corpus."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_point_panel_20260929"))
from verify_panel import check_target, rows, sha  # noqa: E402
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field, verify as verify_rank  # noqa: E402

ARMS = ("ic", "rho_v2", "rho_poly", "rho_normal")
T_CRIT_95_DF4 = 2.7764451051977987


def interval(values: list[float]) -> list[float]:
    assert len(values) == 5 and all(value > 0 for value in values)
    logs = [math.log(value) for value in values]
    mean = statistics.mean(logs)
    half = T_CRIT_95_DF4 * statistics.stdev(logs) / math.sqrt(5)
    return [math.exp(mean - half), math.exp(mean + half)]


def fixture_check(spec: dict) -> tuple[list[dict], Curve, tuple[int, int]]:
    fixture = HERE / spec["fixture_file"]
    points = HERE / spec["points_file"]
    assert sha(fixture) == spec["fixture_sha256"]
    assert sha(points) == spec["points_sha256"]
    known = rows(fixture)
    public = [json.loads(line) for line in points.read_text().splitlines()]
    assert len(known) == len(public) == spec["L"]
    assert len({tuple(q) for q in public}) == spec["L"]
    curve = Curve(Field(spec["n"], spec["field_modulus_low_terms"]), spec["a"])
    generator = tuple(spec["generator"])
    r = spec["subgroup_order"]
    assert curve.on_curve(generator) and curve.mul(r, generator) is None
    for index, (record, q) in enumerate(zip(known, public)):
        assert record["kind"] == "rho_ks_public_fixture"
        assert record["fixture_index"] == index
        assert record["n"] == spec["n"] and record["a"] == spec["a"]
        assert record["corpus"] == spec["corpus"]
        assert record["batch_seed"] == spec["seed"]
        assert record["generator"] == spec["generator"]
        assert record["subgroup_order"] == r
        assert record["field_modulus_low_terms"] == spec["field_modulus_low_terms"]
        assert q == record["published_q"]
        scalar = record["published_fixture_scalar"]
        assert 1 <= scalar < r and curve.on_curve(tuple(q))
        assert curve.mul(scalar, generator) == tuple(q)
    return known, curve, generator


def verify(n: int, length: int, run_dir: Path) -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    for filename, expected in frozen["source_sha256"].items():
        assert sha(ROOT / filename) == expected, filename
    order = tuple(tuple(block) for block in frozen["arm_order"])
    assert frozen["evaluation_repetitions"] == len(order)
    spec = frozen["specs"][f"n{n}_L{length}_eval"]
    known, curve, generator = fixture_check(spec)
    assert run_dir.name == f"n{n}_L{length}"
    host = json.loads((run_dir / "host.json").read_text())
    assert host["cpu_model"] and host["cpu_count"] > 0
    if "Linux" in host["platform"]:
        assert host["mem_total_kib"] > 0
    chosen = json.loads((run_dir / "chosen_k.json").read_text())
    k = frozen["compact_k"][str(n)][str(length)]
    assert chosen == {"n": n, "L": length, "k": k,
                      "selection": "prior_disjoint_tune_no_eval_retuning"}
    runs = json.loads((run_dir / "runs.json").read_text())
    assert len(runs) == len(order) * len(ARMS)
    measurements = {}
    checks = []
    checked_labels: set[tuple[str, int]] = set()
    for block, block_order in enumerate(order):
        group = runs[block * len(ARMS):(block + 1) * len(ARMS)]
        assert [item["arm"] for item in group] == list(block_order)
        measurements[block] = {}
        for item in group:
            arm = item["arm"]
            assert item["block"] == block and arm in ARMS
            assert item["exit_code"] == 0 and not item["timeout"]
            assert item["wall_seconds"] > 0 and item["user_seconds"] >= 0
            assert item["sys_seconds"] >= 0
            stdout = run_dir / item["stdout"]
            stderr = run_dir / item["stderr"]
            assert sha(stdout) == item["stdout_sha256"]
            assert sha(stderr) == item["stderr_sha256"]
            assert not any(".fixture.jsonl" in value
                           for value in item["command"] + list(item["environment"].values()))
            if arm == "ic":
                assert Path(item["environment"]["KIC_DUMP_BASE"]).name == f"b{block}_ic.base.jsonl"
                assert Path(item["environment"]["KIC_DUMP_RANK"]).name == f"b{block}_ic.rank.jsonl"
                base_path = run_dir / f"b{block}_ic.base.jsonl"
                rank_path = run_dir / f"b{block}_ic.rank.jsonl"
                target_path = run_dir / f"b{block}_ic.target.jsonl"
                rank = verify_rank(rank_path, base_path, stdout)
                assert rank["status"] == "PASS" and rank["rank"] == k
                base, = rows(base_path)
                assert base["subgroup_order"] == spec["subgroup_order"]
                assert base["field_modulus_low_terms"] == spec["field_modulus_low_terms"]
                logs = rows(rank_path)[-1]["logs"]
                summary, = rows(stdout)
                assert summary["rank"] == k
                assert summary["targets"] == summary["targets_solved"] == length
                assert summary["targets_failed"] == 0
                assert summary["base_hash"] == base["base_hash"]
                targets = rows(target_path)
                assert len(targets) == length
                for record, fixture in zip(targets, known):
                    check_target(record, fixture, base, logs, curve, generator,
                                 checked_labels)
                checks.append({"block": block, "arm": arm,
                               "targets_verified": length,
                               "rank_receipt": rank,
                               "rank_attempts": summary["rank_attempts"],
                               "rank_probes_mean": summary["rank_probes_mean"],
                               "counting_floor_attempt_ratio":
                               (summary["rank_attempts"] + length) / (k + length)})
            else:
                assert item["environment"]["KIC_RHO_BATCH_CORPUS"] == spec["corpus"]
                assert item["environment"]["KIC_RHO_DP_BITS"] == str(frozen["rho_dp_bits"])
                assert Path(item["environment"]["KIC_RHO_POINT_INPUT"]).name == Path(spec["points_file"]).name
                expected_backend = ("poly_xfirst" if arm == "rho_poly"
                                    else "normal_basis" if arm == "rho_normal" else None)
                assert item["environment"].get("KIC_RHO_CANON_BACKEND") == expected_backend
                data = rows(stdout)
                assert len(data) == length + 1
                summary = data[-1]
                assert summary["kind"] == "rho_ks_batch_summary"
                assert summary["target_source"] == "public_point_jsonl"
                assert summary["quotient_mode"] == "signed_frobenius"
                assert summary["fixtures"] == length and summary["all_verified"]
                assert summary["corpus"] == spec["corpus"]
                assert summary["dp_bits"] == frozen["rho_dp_bits"]
                assert summary["inversion_backend"] == "itoh_tsujii"
                if "x86_64" in host["platform"] and "pclmulqdq" in host["cpu_flags"]:
                    assert summary["field_product_backend"] == "pclmulqdq"
                if arm != "rho_v2":
                    assert summary["parallel_walks"] == 32
                    assert summary["canonicalization_backend"] == expected_backend
                    charges = summary["charges"]
                    assert charges["partition_hashes"] == summary["total_walk_steps"]
                    assert (charges["batch_inversion_inputs"]
                            + charges["batch_fallback_additions"]
                            == summary["total_walk_steps"])
                    assert charges["batch_inversion_calls"] > 0
                    if arm == "rho_normal":
                        assert charges["normal_basis_transforms"] > 0
                        assert charges["canonical_rotation_steps"] > 0
                    else:
                        assert charges["normal_basis_transforms"] == 0
                for index, (record, fixture) in enumerate(zip(data[:-1], known)):
                    assert record["kind"] == "rho_ks_batch_fixture"
                    assert record["fixture_index"] == index
                    assert record["published_fixture_scalar"] is None
                    assert record["target_source"] == "public_point_jsonl"
                    assert record["published_q"] == fixture["published_q"]
                    assert record["recovered_fixture_scalar"] == fixture["published_fixture_scalar"]
                    assert curve.mul(record["recovered_fixture_scalar"], generator) == tuple(record["published_q"])
                checks.append({"block": block, "arm": arm,
                               "targets_verified": length,
                               "walk_steps": summary["total_walk_steps"],
                               "charges": summary["charges"]})
            measurements[block][arm] = {
                "wall_seconds": item["wall_seconds"],
                "cpu_seconds": item["user_seconds"] + item["sys_seconds"],
                "rss_kib": item["max_rss_kib_linux"],
            }
    paired = {}
    for arm in ARMS:
        if arm == "rho_normal":
            continue
        wall = [measurements[b][arm]["wall_seconds"] /
                measurements[b]["rho_normal"]["wall_seconds"] for b in measurements]
        cpu = [measurements[b][arm]["cpu_seconds"] /
               measurements[b]["rho_normal"]["cpu_seconds"] for b in measurements]
        paired[arm] = {"wall_ratio_median": statistics.median(wall),
                       "wall_ratio_95pct_log_t_interval": interval(wall),
                       "cpu_ratio_median": statistics.median(cpu),
                       "cpu_ratio_95pct_log_t_interval": interval(cpu)}
    return {
        "status": "PASS", "n": n, "L": length, "k": k,
        "targets_verified_per_arm": length * len(order),
        "rank_runs_verified": len(order),
        "checks": checks, "measurements": measurements,
        "paired_to_normal_rho": paired,
        "S_ic": None, "S_rho_v2": None, "S_rho_poly": None,
        "S_rho_normal": None,
        "evidence_class": "cold_wall_and_cpu_diagnostic_no_common_operation_unit",
        "miss_nonexistence_proved": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=(37, 41, 53), required=True)
    parser.add_argument("--L", type=int, choices=(1, 1024), required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "refusing to overwrite verification receipt"
    try:
        receipt = verify(args.n, args.L, args.run_dir.resolve())
    except BaseException as error:
        receipt = {"status": "FAIL", "error_type": type(error).__name__,
                   "error": str(error), "traceback": traceback.format_exc()}
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        raise
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({name: value for name, value in receipt.items()
                      if name not in ("checks", "measurements")}, sort_keys=True))


if __name__ == "__main__":
    main()
