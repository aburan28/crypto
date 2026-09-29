#!/usr/bin/env python3
"""Independent replay of baseline and swap-quotient K grids against strong rho."""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
import statistics
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
POINT_VERIFIER = (ROOT / "research/notes/ecc2k130/compact_orbit_point_panel_20260929"
                  / "verify_panel.py")
point_spec = importlib.util.spec_from_file_location("compact_point_verifier", POINT_VERIFIER)
assert point_spec is not None and point_spec.loader is not None
point_module = importlib.util.module_from_spec(point_spec)
point_spec.loader.exec_module(point_module)
check_target, rows, sha = point_module.check_target, point_module.rows, point_module.sha
sys.path.insert(0, str(ROOT / "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929"))
from verify_rank import Curve, Field, verify as verify_rank  # noqa: E402

T_CRIT_95_DF8 = 2.306004135033371


def interval(values: list[float]) -> list[float]:
    assert len(values) == 9 and all(value > 0 for value in values)
    logs = [math.log(value) for value in values]
    mean = statistics.mean(logs)
    half = T_CRIT_95_DF8 * statistics.stdev(logs) / math.sqrt(9)
    return [math.exp(mean - half), math.exp(mean + half)]


def fixture_check(spec: dict) -> tuple[list[dict], Curve, tuple[int, int]]:
    fixture = HERE / spec["fixture_file"]
    points = HERE / spec["points_file"]
    assert sha(fixture) == spec["fixture_sha256"]
    assert sha(points) == spec["points_sha256"]
    known = rows(fixture)
    public = [json.loads(line) for line in points.read_text().splitlines()]
    assert len(known) == len(public) == spec["L"] == 1024
    new_points = {tuple(q) for q in public}
    assert len(new_points) == 1024
    assert spec["prior_overlap_count"] == 0
    prior_points = set()
    for filename, digest in spec["prior_point_sha256"].items():
        path = ROOT / filename
        assert sha(path) == digest
        prior_points.update(tuple(json.loads(line))
                            for line in path.read_text().splitlines())
    assert len(prior_points) == spec["prior_points_checked"]
    assert not (new_points & prior_points)
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


def verify(n: int, run_dir: Path) -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    source_lock_path = HERE / "SOURCE_FROZEN.json"
    assert sha(source_lock_path) == frozen["source_lock_sha256"]
    assert frozen["schema"] == "compact-swap-quotient-evaluation-freeze-v1"
    for key, value in json.loads(source_lock_path.read_text()).items():
        if key != "schema":
            assert frozen[key] == value, key
    for filename, expected in frozen["source_sha256"].items():
        assert sha(ROOT / filename) == expected, filename
    spec = frozen["specs"][f"n{n}_L1024_eval"]
    known, curve, generator = fixture_check(spec)
    assert run_dir.name == f"n{n}_L1024"
    host = json.loads((run_dir / "host.json").read_text())
    assert host["cpu_model"] and host["cpu_count"] > 0
    if "Linux" in host["platform"]:
        assert host["mem_total_kib"] > 0
    grid = frozen["k_grid"][str(n)]
    assert json.loads((run_dir / "grid.json").read_text()) == {
        "n": n, "L": 1024, "k_grid": grid,
        "selection": "pre_registered_swap_quotient_on_prior_disjoint_corpus",
        "policies": ["base", "swap", "rho_normal"]}
    names = tuple(f"{policy}_k{k}" for k in grid
                  for policy in ("base", "swap")) + ("rho_normal",)
    assert len(names) == 9
    runs = json.loads((run_dir / "runs.json").read_text())
    assert len(runs) == 9 * len(names)
    measurements = {}
    checks = []
    bases = {}
    summaries = {}
    checked_labels: set[tuple[str, int]] = set()
    for block, offset in enumerate(frozen["arm_order_offsets"]):
        order = names[offset:] + names[:offset]
        group = runs[block * len(names):(block + 1) * len(names)]
        assert [item["arm"] for item in group] == list(order)
        measurements[block] = {}
        for item in group:
            arm = item["arm"]
            assert item["block"] == block and arm in names
            assert item["exit_code"] == 0 and not item["timeout"]
            assert item["wall_seconds"] > 0 and item["user_seconds"] >= 0
            assert item["sys_seconds"] >= 0
            if "Linux" in host["platform"]:
                assert item["max_rss_kib_linux"] > 0
            stdout = run_dir / item["stdout"]
            stderr = run_dir / item["stderr"]
            assert sha(stdout) == item["stdout_sha256"]
            assert sha(stderr) == item["stderr_sha256"]
            assert not any(".fixture.jsonl" in value for value in
                           item["command"] + list(item["environment"].values()))
            if arm != "rho_normal":
                policy, k_text = arm.split("_k")
                assert policy in ("base", "swap")
                k = int(k_text)
                assert k in grid
                binary_digest = host["baseline_binary_sha256"] if policy == "base" else host["swap_binary_sha256"]
                assert len(binary_digest) == 64
                assert Path(item["command"][0]).name == (
                    "koblitz_orbit_dlp_fast" if policy == "base"
                    else "koblitz_orbit_dlp_swap")
                assert item["command"][1] == f"construct:{n}:0:{k}"
                base_path = run_dir / f"b{block}_{arm}.base.jsonl"
                rank_path = run_dir / f"b{block}_{arm}.rank.jsonl"
                target_path = run_dir / f"b{block}_{arm}.target.jsonl"
                assert Path(item["environment"]["KIC_DUMP_BASE"]).name == base_path.name
                assert Path(item["environment"]["KIC_DUMP_RANK"]).name == rank_path.name
                rank = verify_rank(rank_path, base_path, stdout)
                assert rank["status"] == "PASS" and rank["rank"] == k
                base, = rows(base_path)
                bases[(block, k, policy)] = base
                logs = rows(rank_path)[-1]["logs"]
                summary, = rows(stdout)
                summaries[(block, k, policy)] = summary
                assert summary["rank"] == k
                if policy == "swap":
                    assert summary["index_policy"] == "swap_frobenius_quotient"
                    assert summary["ordered_state_candidates"] == k * k * n
                    assert summary["representative_state_candidates"] == (k * k * n + k) // 2
                    assert summary["root_table_slots"] > summary["root_table_entries"]
                else:
                    assert "index_policy" not in summary
                assert summary["targets"] == summary["targets_solved"] == 1024
                assert summary["targets_failed"] == 0
                assert summary["base_hash"] == base["base_hash"]
                assert summary["factor_base_points"] == k * 2 * n
                targets = rows(target_path)
                assert len(targets) == 1024
                assert all(isinstance(record["probes"], int) and record["probes"] >= 0
                           for record in targets)
                for record, fixture in zip(targets, known):
                    check_target(record, fixture, base, logs, curve, generator,
                                 checked_labels)
                checks.append({
                    "block": block, "arm": arm, "policy": policy, "rank": k,
                    "targets_verified": 1024, "rank_receipt": rank,
                    "rank_attempts": summary["rank_attempts"],
                    "rank_probes_mean": summary["rank_probes_mean"],
                    "target_probes_mean": sum(record["probes"] for record in targets) / 1024,
                    "counting_floor_attempt_ratio":
                        (summary["rank_attempts"] + 1024) / (k + 1024),
                    "timing_ms": summary["timing_ms"],
                    "regular_states": summary["regular_states"],
                    "root_table_entries": summary["root_table_entries"],
                    "root_table_slots": summary.get("root_table_slots"),
                    "factor_base_points": summary["factor_base_points"],
                })
            else:
                assert arm == "rho_normal"
                assert item["environment"]["KIC_RHO_CANON_BACKEND"] == "normal_basis"
                assert item["environment"]["KIC_RHO_BATCH_CORPUS"] == spec["corpus"]
                assert item["environment"]["KIC_RHO_DP_BITS"] == str(frozen["rho_dp_bits"])
                assert Path(item["environment"]["KIC_RHO_POINT_INPUT"]).name == Path(spec["points_file"]).name
                data = rows(stdout)
                assert len(data) == 1025
                summary = data[-1]
                assert summary["kind"] == "rho_ks_batch_summary"
                assert summary["target_source"] == "public_point_jsonl"
                assert summary["quotient_mode"] == "signed_frobenius"
                assert summary["canonicalization_backend"] == "normal_basis"
                assert summary["parallel_walks"] == 32
                assert summary["fixtures"] == 1024 and summary["all_verified"]
                assert summary["corpus"] == spec["corpus"]
                assert summary["dp_bits"] == frozen["rho_dp_bits"]
                assert summary["inversion_backend"] == "itoh_tsujii"
                if "x86_64" in host["platform"] and "pclmulqdq" in host["cpu_flags"]:
                    assert summary["field_product_backend"] == "pclmulqdq"
                charges = summary["charges"]
                assert charges["partition_hashes"] == summary["total_walk_steps"]
                assert (charges["batch_inversion_inputs"]
                        + charges["batch_fallback_additions"]
                        == summary["total_walk_steps"])
                assert charges["batch_inversion_calls"] > 0
                assert charges["normal_basis_transforms"] > 0
                for index, (record, fixture) in enumerate(zip(data[:-1], known)):
                    assert record["kind"] == "rho_ks_batch_fixture"
                    assert record["fixture_index"] == index
                    assert record["published_fixture_scalar"] is None
                    assert record["target_source"] == "public_point_jsonl"
                    assert record["published_q"] == fixture["published_q"]
                    assert record["recovered_fixture_scalar"] == fixture["published_fixture_scalar"]
                    assert curve.mul(record["recovered_fixture_scalar"], generator) == tuple(record["published_q"])
                checks.append({"block": block, "arm": arm,
                               "targets_verified": 1024,
                               "walk_steps": summary["total_walk_steps"],
                               "charges": charges})
            measurements[block][arm] = {
                "wall_seconds": item["wall_seconds"],
                "cpu_seconds": item["user_seconds"] + item["sys_seconds"],
                "rss_kib": item["max_rss_kib_linux"],
            }
    # Equal-useful-size base construction and the root-key cardinality must
    # agree for each same-K pair. Only the enumerated representation changes.
    for block in range(9):
        for k in grid:
            assert bases[(block, k, "base")] == bases[(block, k, "swap")]
            baseline = summaries[(block, k, "base")]
            quotient = summaries[(block, k, "swap")]
            assert baseline["root_table_entries"] == quotient["root_table_entries"]
            fixed_regular = 2 * quotient["regular_states"] - baseline["regular_states"]
            assert 0 <= fixed_regular <= k
            assert baseline["rank"] == quotient["rank"] == k

    paired = {}
    for k in grid:
        base_arm, swap_arm = f"base_k{k}", f"swap_k{k}"
        for label, numerator, denominator in (
            (f"base_k{k}_to_rho", base_arm, "rho_normal"),
            (f"swap_k{k}_to_rho", swap_arm, "rho_normal"),
            (f"swap_k{k}_to_base", swap_arm, base_arm),
        ):
            cpu = [measurements[b][numerator]["cpu_seconds"] /
                   measurements[b][denominator]["cpu_seconds"] for b in range(9)]
            wall = [measurements[b][numerator]["wall_seconds"] /
                    measurements[b][denominator]["wall_seconds"] for b in range(9)]
            paired[label] = {
                "cpu_ratio_median": statistics.median(cpu),
                "cpu_ratio_95pct_log_t_interval": interval(cpu),
                "wall_ratio_median": statistics.median(wall),
                "wall_ratio_95pct_log_t_interval": interval(wall),
            }
    return {
        "status": "PASS", "n": n, "L": 1024, "k_grid": grid,
        "targets_verified_per_arm": 9216,
        "rank_runs_verified": 72,
        "checks": checks, "measurements": measurements,
        "paired": paired,
        "S_base": None, "S_swap": None, "S_rho_normal": None,
        "evidence_class": "cold_full_process_diagnostic_no_common_operation_unit",
        "miss_nonexistence_proved": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, choices=(41, 53), required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a verification receipt"
    try:
        receipt = verify(args.n, args.run_dir.resolve())
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
