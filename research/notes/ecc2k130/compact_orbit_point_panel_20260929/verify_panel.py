#!/usr/bin/env python3
"""Independent full-rank, relation, and known-answer replay of a paired panel."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import statistics
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field, point, verify as verify_rank  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def check_fixture(spec: dict) -> tuple[list[dict], Curve, tuple[int, int]]:
    fixture = HERE / spec["fixture_file"]
    points = HERE / spec["points_file"]
    assert sha(fixture) == spec["fixture_sha256"]
    assert sha(points) == spec["points_sha256"]
    known = rows(fixture)
    public = [json.loads(line) for line in points.read_text().splitlines()]
    assert len(known) == len(public) == spec["L"]
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
        assert 1 <= scalar < r
        assert curve.on_curve(tuple(q))
        assert curve.mul(scalar, generator) == tuple(q)
    return known, curve, generator


def check_target(record: dict, fixture: dict, base: dict, logs: list[int],
                 curve: Curve, generator: tuple[int, int],
                 checked_labels: set[tuple[str, int]]) -> None:
    q = tuple(fixture["published_q"])
    scalar = fixture["published_fixture_scalar"]
    r = base["subgroup_order"]
    assert record["kind"] == "compact_orbit_dlp_target"
    assert record["fixture_index"] == fixture["fixture_index"]
    assert record["n"] == fixture["n"] and record["a"] == fixture["a"]
    assert record["target"] == record["published_q"] == list(q)
    assert record["generator"] == list(generator)
    assert record["published_fixture_scalar"] is None
    assert record["recovered_matches_published"] is None
    assert record["target_generation_ms_excluded"] == 0
    assert record["exit_code"] == 0
    assert record["group_verified"] is True
    assert record["recovered_scalar"] == scalar
    # The independently checked fixture has [scalar]G=Q; equality to that
    # scalar proves this recovered log for the exact same public point.
    indices = record["point_indices"]
    assert isinstance(indices, list) and len(indices) == 4
    codes = record["x_codes"]
    assert isinstance(codes, list) and len(codes) == 4
    base_points = base["factor_base_point_coordinates"]
    labels = base["factor_base_point_labels"]
    chosen = []
    recomputed_log = 0
    for index, code in zip(indices, codes):
        assert 0 <= index < len(base_points)
        selected = point(base_points[index])
        assert selected is not None and selected[0] == code
        assert curve.on_curve(selected)
        column, coefficient = labels[index]
        assert 0 <= column < len(logs)
        label_key = (base["base_hash"], index)
        if label_key not in checked_labels:
            assert curve.mul(coefficient, point(base["factor_base_representatives"][column])) == selected
            checked_labels.add(label_key)
        recomputed_log = (recomputed_log + coefficient * logs[column]) % r
        chosen.append(selected)
    left = curve.add(chosen[0], chosen[1])
    right = curve.add(chosen[2], chosen[3])
    assert left is not None and right is not None
    assert curve.add(left, right) == q
    assert set(record["pinned_intermediates"]) == {left[0], right[0]}
    assert recomputed_log == scalar


def check_run(n: int, length: int, run_dir: Path) -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    for filename, expected in frozen["source_sha256"].items():
        assert sha(ROOT / filename) == expected, filename
    spec = frozen["specs"][f"n{n}_L{length}_eval"]
    known, curve, generator = check_fixture(spec)
    assert run_dir.name == f"n{n}_L{length}"
    chosen = json.loads((run_dir / "chosen_k.json").read_text())
    k = chosen["k"]
    assert chosen["n"] == n and chosen["L"] == length
    if length == 1024:
        assert k == frozen["batch_k"][str(n)]
    else:
        assert k in frozen["single_k_grid"][str(n)]
        tune = json.loads((run_dir / "tune_scores.json").read_text())
        assert tune["chosen_k"] == k
        assert len(tune["runs"]) == len(frozen["single_k_grid"][str(n)]) * frozen["tuning_repetitions"]
        assert k == min((int(v) for v in tune["scores_cpu_seconds"]),
                        key=lambda candidate:
                        (tune["scores_cpu_seconds"][str(candidate)], candidate))
        tune_spec = frozen["specs"][f"n{n}_L1_tune"]
        check_fixture(tune_spec)
    runs = json.loads((run_dir / "runs.json").read_text())
    assert len(runs) == frozen["evaluation_repetitions"] * 2
    checks = []
    measurements = {}
    checked_labels: set[tuple[str, int]] = set()
    for block in range(frozen["evaluation_repetitions"]):
        expected_order = ("ic", "rho") if block % 2 == 0 else ("rho", "ic")
        assert [r["arm"] for r in runs[2 * block:2 * block + 2]] == list(expected_order)
        measurements[block] = {}
        for item in runs[2 * block:2 * block + 2]:
            arm = item["arm"]
            assert item["block"] == block
            assert item["exit_code"] == 0 and not item["timeout"]
            assert item["wall_seconds"] > 0 and item["user_seconds"] >= 0
            assert item["sys_seconds"] >= 0
            stdout = run_dir / item["stdout"]
            stderr = run_dir / item["stderr"]
            assert sha(stdout) == item["stdout_sha256"]
            assert sha(stderr) == item["stderr_sha256"]
            assert not any(".fixture.jsonl" in value
                           for value in item["command"] + list(item["environment"].values()))
            if arm == "rho":
                assert item["environment"]["KIC_RHO_BATCH_CORPUS"] == spec["corpus"]
                assert item["environment"]["KIC_RHO_DP_BITS"] == str(frozen["rho_dp_bits"])
                assert Path(item["environment"]["KIC_RHO_POINT_INPUT"]).name == Path(spec["points_file"]).name
                data = rows(stdout)
                assert len(data) == length + 1
                summary = data[-1]
                assert summary["kind"] == "rho_ks_batch_summary"
                assert summary["target_source"] == "public_point_jsonl"
                assert summary["quotient_mode"] == "signed_frobenius"
                assert summary["all_verified"] is True
                assert summary["fixtures"] == length
                assert summary["corpus"] == spec["corpus"]
                assert summary["dp_bits"] == frozen["rho_dp_bits"]
                for index, (record, fixture) in enumerate(zip(data[:-1], known)):
                    assert record["kind"] == "rho_ks_batch_fixture"
                    assert record["fixture_index"] == index
                    assert record["published_fixture_scalar"] is None
                    assert record["target_source"] == "public_point_jsonl"
                    assert record["published_q"] == fixture["published_q"]
                    assert record["recovered_fixture_scalar"] == fixture["published_fixture_scalar"]
                    # check_fixture independently proved the matching known-answer
                    # scalar times G equals this exact point.
                checks.append({"block": block, "arm": arm, "targets_verified": length,
                               "charges": summary["charges"]})
            else:
                assert Path(item["environment"]["KIC_DUMP_BASE"]).name == f"b{block}_ic.base.jsonl"
                assert Path(item["environment"]["KIC_DUMP_RANK"]).name == f"b{block}_ic.rank.jsonl"
                target_path = run_dir / f"b{block}_ic.target.jsonl"
                target_data = rows(target_path)
                assert len(target_data) == length
                base_path = run_dir / f"b{block}_ic.base.jsonl"
                rank_path = run_dir / f"b{block}_ic.rank.jsonl"
                rank = verify_rank(rank_path, base_path, stdout)
                assert rank["status"] == "PASS" and rank["rank"] == k
                base, = rows(base_path)
                assert base["subgroup_order"] == spec["subgroup_order"]
                assert base["field_modulus_low_terms"] == spec["field_modulus_low_terms"]
                solution = rows(rank_path)[-1]
                logs = solution["logs"]
                summary, = rows(stdout)
                assert summary["targets"] == summary["targets_solved"] == length
                assert summary["targets_failed"] == 0 and summary["rank"] == k
                assert summary["base_hash"] == base["base_hash"]
                for record, fixture in zip(target_data, known):
                    check_target(record, fixture, base, logs, curve, generator,
                                 checked_labels)
                checks.append({"block": block, "arm": arm, "targets_verified": length,
                               "rank_receipt": rank})
            measurements[block][arm] = {
                "wall_seconds": item["wall_seconds"],
                "cpu_seconds": item["user_seconds"] + item["sys_seconds"],
                "rss_kib": item["max_rss_kib_linux"],
            }
    wall = [measurements[b]["ic"]["wall_seconds"] /
            measurements[b]["rho"]["wall_seconds"] for b in measurements]
    cpu = [measurements[b]["ic"]["cpu_seconds"] /
           measurements[b]["rho"]["cpu_seconds"] for b in measurements]
    assert len(wall) == len(cpu) == frozen["evaluation_repetitions"]
    return {
        "status": "PASS", "n": n, "L": length, "k": k,
        "targets_verified_per_arm": length * frozen["evaluation_repetitions"],
        "rank_runs_verified": frozen["evaluation_repetitions"],
        "checks": checks, "measurements": measurements,
        "wall_ratio_median": statistics.median(wall),
        "wall_ratio_range": [min(wall), max(wall)],
        "cpu_ratio_median": statistics.median(cpu),
        "cpu_ratio_range": [min(cpu), max(cpu)],
        "S_ic": None, "S_rho": None, "common_operation_ratio": None,
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
        receipt = check_run(args.n, args.L, args.run_dir.resolve())
    except BaseException as error:
        receipt = {"status": "FAIL", "error_type": type(error).__name__,
                   "error": str(error), "traceback": traceback.format_exc()}
        args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        raise
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({name: value for name, value in receipt.items()
                      if name not in ("checks", "measurements")}, sort_keys=True))


if __name__ == "__main__":
    main()
