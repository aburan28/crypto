#!/usr/bin/env python3
"""Independently replay every rank and public-point log in the filter panel."""
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

ARMS = ("off_a", "filter", "off_b", "rho_normal")
COUNT_KEYS = (
    "calls", "regular_calls", "exceptional_calls", "returned_pairs", "no_root_returns",
    "scalar_inversions", "batch_inversions", "batch_multiplications",
    "consumed_candidates", "discarded_candidates", "table_lookups", "table_lookup_probes",
    "table_insert_probes", "lift_attempts", "root_keys_considered", "prefilter_checks",
    "prefilter_definite_misses", "prefilter_false_positives", "prefilter_true_positives",
)
TARGET_IDENTITY_KEYS = (
    "target", "published_q", "recovered_scalar", "group_verified", "point_indices",
    "x_codes", "pinned_intermediates", "probes",
)
T_CRIT_95_DF4 = 2.7764451051977987


def interval(ratios: list[float]) -> list[float]:
    assert len(ratios) == 5 and all(value > 0 for value in ratios)
    logs = [math.log(value) for value in ratios]
    half = T_CRIT_95_DF4 * statistics.stdev(logs) / math.sqrt(5)
    mean = statistics.mean(logs)
    return [math.exp(mean - half), math.exp(mean + half)]


def check_counts(counts: dict, filtered: bool) -> None:
    assert set(counts) == set(COUNT_KEYS)
    assert all(isinstance(counts[key], int) and counts[key] >= 0 for key in COUNT_KEYS)
    assert counts["calls"] == counts["regular_calls"] + counts["exceptional_calls"]
    assert counts["calls"] == counts["returned_pairs"] + counts["no_root_returns"]
    assert counts["calls"] == counts["consumed_candidates"] + counts["discarded_candidates"]
    assert counts["table_lookup_probes"] >= counts["table_lookups"]
    assert counts["root_keys_considered"] == (counts["table_lookups"]
            + counts["prefilter_definite_misses"])
    assert counts["prefilter_checks"] == (counts["prefilter_definite_misses"]
            + counts["prefilter_false_positives"] + counts["prefilter_true_positives"])
    assert counts["scalar_inversions"] == 0
    assert counts["batch_multiplications"] == 3 * counts["regular_calls"]
    if counts["regular_calls"]:
        assert 0 < counts["batch_inversions"] <= counts["regular_calls"]
    else:
        assert counts["batch_inversions"] == 0
    if filtered:
        assert counts["prefilter_checks"] == counts["root_keys_considered"]
        assert (counts["prefilter_false_positives"] + counts["prefilter_true_positives"]
                == counts["table_lookups"])
    else:
        assert counts["prefilter_checks"] == counts["prefilter_definite_misses"] == 0
        assert counts["prefilter_false_positives"] == counts["prefilter_true_positives"] == 0


def check_fixture(spec: dict) -> tuple[list[dict], Curve, tuple[int, int]]:
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
        prior_points.update(tuple(json.loads(line)) for line in path.read_text().splitlines())
    assert len(prior_points) == spec["prior_points_checked"]
    assert not new_points & prior_points
    curve = Curve(Field(spec["n"], spec["field_modulus_low_terms"]), spec["a"])
    generator = tuple(spec["generator"])
    r = spec["subgroup_order"]
    assert curve.on_curve(generator) and curve.mul(r, generator) is None
    for index, (record, q) in enumerate(zip(known, public)):
        assert record["kind"] == "rho_ks_public_fixture"
        assert record["fixture_index"] == index
        assert record["n"] == spec["n"] and record["a"] == spec["a"]
        assert record["corpus"] == spec["corpus"] and record["batch_seed"] == spec["seed"]
        assert record["generator"] == spec["generator"]
        assert record["subgroup_order"] == r
        assert record["field_modulus_low_terms"] == spec["field_modulus_low_terms"]
        assert record["published_q"] == q
        scalar = record["published_fixture_scalar"]
        assert 1 <= scalar < r and curve.on_curve(tuple(q))
        assert curve.mul(scalar, generator) == tuple(q)
    return known, curve, generator


def freeze() -> dict:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    lock = HERE / "SOURCE_FROZEN.json"
    assert sha(lock) == frozen["source_lock_sha256"]
    assert frozen["schema"] == "compact-s3-prefilter-evaluation-freeze-v1"
    for key, value in json.loads(lock.read_text()).items():
        if key != "schema":
            assert frozen[key] == value, key
    for filename, digest in frozen["source_sha256"].items():
        assert sha(ROOT / filename) == digest, filename
    return frozen


def target_identity(record: dict) -> tuple[str, ...]:
    return tuple(json.dumps(record[key], sort_keys=True) for key in TARGET_IDENTITY_KEYS)


def paired(measurements: dict, numerator: str, denominator: str) -> dict:
    cpu = [measurements[b][numerator]["cpu_seconds"] / measurements[b][denominator]["cpu_seconds"]
           for b in range(5)]
    wall = [measurements[b][numerator]["wall_seconds"] / measurements[b][denominator]["wall_seconds"]
            for b in range(5)]
    return {"cpu_ratio_median": statistics.median(cpu),
            "cpu_ratio_95pct_log_t_interval": interval(cpu),
            "wall_ratio_median": statistics.median(wall),
            "wall_ratio_95pct_log_t_interval": interval(wall)}


def verify(n: int, run_dir: Path) -> dict:
    frozen = freeze()
    spec = frozen["specs"][f"n{n}_L1024_eval"]
    known, curve, generator = check_fixture(spec)
    assert run_dir.name == f"n{n}_L1024"
    host = json.loads((run_dir / "host.json").read_text())
    assert host["cpu_model"] and host["cpu_count"] > 0
    isolation, = rows(run_dir / "isolation.jsonl")
    assert isolation["schema"] == "isolated-bench/1" and isolation["mode"] == "reserve"
    assert isolation["exit_status"] == 0
    assert isolation["reserved_cpus"] == [host["reserved_cpu"]]
    assert isolation["command"][-2:] == ["--cpu", str(host["reserved_cpu"])]
    assert json.loads((run_dir / "grid.json").read_text()) == {
        "n": n, "L": 1024, "k": frozen["k_by_n"][str(n)],
        "selection": "pre_registered_W64_incumbent_on_new_disjoint_corpus",
        "policies": list(ARMS),
    }
    k = frozen["k_by_n"][str(n)]
    runs = json.loads((run_dir / "runs.json").read_text())
    assert len(runs) == 5 * len(ARMS)
    measurements: dict[int, dict] = {}
    checks = []
    evidence: dict[tuple[int, str], tuple] = {}
    checked_labels: set[tuple[str, int]] = set()
    for block, offset in enumerate(frozen["arm_order_offsets"]):
        order = ARMS[offset:] + ARMS[:offset]
        group = runs[block * len(ARMS):(block + 1) * len(ARMS)]
        assert [item["arm"] for item in group] == list(order)
        measurements[block] = {}
        for item in group:
            arm = item["arm"]
            assert item["block"] == block and arm in ARMS
            assert item["exit_code"] == 0 and not item["timeout"], (block, arm, item["exit_code"])
            assert 0 < item["wall_seconds"] <= frozen["per_arm_timeout_seconds"] + 2
            assert item["user_seconds"] >= 0 and item["sys_seconds"] >= 0
            assert item["max_rss_kib_linux"] > 0
            assert item["max_rss_kib_linux"] * 1024 < frozen["per_arm_address_space_limit_bytes"]
            command = item["command"]
            assert command[:3] == ["taskset", "-c", str(host["reserved_cpu"])]
            assert not any(".fixture.jsonl" in value for value in
                           command + list(item["environment"].values()))
            stdout, stderr = run_dir / item["stdout"], run_dir / item["stderr"]
            assert sha(stdout) == item["stdout_sha256"]
            assert sha(stderr) == item["stderr_sha256"]
            assert not stderr.read_text().strip(), (block, arm, stderr.read_text()[:1000])
            if arm != "rho_normal":
                filtered = arm == "filter"
                assert Path(command[3]).name == "koblitz_orbit_dlp_s3_batch"
                assert command[4:7] == [f"construct:{n}:0:{k}",
                        str((HERE / spec["points_file"]).resolve()), "7"]
                assert item["environment"]["KIC_S3_BATCH_WINDOW"] == "64"
                assert item["environment"]["KIC_S3_PREFILTER"] == ("blocked" if filtered else "off")
                base_path = run_dir / f"b{block}_{arm}.base.jsonl"
                rank_path = run_dir / f"b{block}_{arm}.rank.jsonl"
                target_path = run_dir / f"b{block}_{arm}.target.jsonl"
                assert Path(item["environment"]["KIC_DUMP_BASE"]).name == base_path.name
                assert Path(item["environment"]["KIC_DUMP_RANK"]).name == rank_path.name
                assert Path(command[7]).name == target_path.name
                rank = verify_rank(rank_path, base_path, stdout)
                assert rank["status"] == "PASS" and rank["rank"] == k
                base, = rows(base_path)
                rank_rows = rows(rank_path)
                logs = rank_rows[-1]["logs"]
                summary, = rows(stdout)
                assert summary["kind"] == "compact_orbit_dlp_summary"
                assert summary["rank"] == k and summary["s3_batch_window"] == 64
                assert summary["index_policy"] == "swap_frobenius_quotient"
                assert summary["root_prefilter_policy"] == (
                    "blocked_bloom_512_3hash" if filtered else "off")
                assert (summary["root_prefilter_bytes"] > 0) == filtered
                assert summary["root_prefilter_bytes"] == 64 * summary["root_prefilter_blocks"]
                assert summary["ordered_state_candidates"] == k * k * n
                assert summary["representative_state_candidates"] == (k * k * n + k) // 2
                assert summary["root_table_slots"] > summary["root_table_entries"]
                assert summary["targets"] == summary["targets_solved"] == 1024
                assert summary["targets_failed"] == 0 and summary["base_hash"] == base["base_hash"]
                assert summary["factor_base_points"] == k * 2 * n
                for phase in ("index", "rank", "target"):
                    check_counts(summary[f"{phase}_s3_counts"], filtered)
                assert summary["index_s3_counts"]["calls"] == summary["representative_state_candidates"]
                assert summary["index_s3_counts"]["returned_pairs"] == summary["regular_states"]
                targets = rows(target_path)
                assert len(targets) == 1024
                aggregate = {key: 0 for key in COUNT_KEYS}
                for record, fixture in zip(targets, known):
                    check_target(record, fixture, base, logs, curve, generator, checked_labels)
                    counts = record["s3_counts"]
                    check_counts(counts, filtered)
                    assert record["probes"] == counts["root_keys_considered"]
                    for key in COUNT_KEYS:
                        aggregate[key] += counts[key]
                assert aggregate == summary["target_s3_counts"]
                evidence[(block, arm)] = (base, rank_rows,
                                          tuple(target_identity(record) for record in targets),
                                          summary)
                checks.append({"block": block, "arm": arm, "rank_receipt": rank,
                               "targets_verified": 1024, "base_sha256": sha(base_path),
                               "rank_sha256": sha(rank_path), "target_sha256": sha(target_path),
                               "timing_ms": summary["timing_ms"],
                               "root_prefilter_bytes": summary["root_prefilter_bytes"],
                               "rank_s3_counts": summary["rank_s3_counts"],
                               "target_s3_counts": summary["target_s3_counts"]})
            else:
                assert Path(command[3]).name == "koblitz_rho_batch_ks_v3"
                assert command[4:] == [str(n), "0", "signed_frobenius", "1024", str(spec["seed"])]
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
                assert summary["parallel_walks"] == 32 and summary["fixtures"] == 1024
                assert summary["all_verified"] and summary["corpus"] == spec["corpus"]
                assert summary["dp_bits"] == frozen["rho_dp_bits"]
                assert summary["inversion_backend"] == "itoh_tsujii"
                if "x86_64" in host["platform"] and "pclmulqdq" in host["cpu_flags"]:
                    assert summary["field_product_backend"] == "pclmulqdq"
                charges = summary["charges"]
                assert charges["partition_hashes"] == summary["total_walk_steps"]
                assert (charges["batch_inversion_inputs"] + charges["batch_fallback_additions"]
                        == summary["total_walk_steps"])
                for index, (record, fixture) in enumerate(zip(data[:-1], known)):
                    assert record["kind"] == "rho_ks_batch_fixture"
                    assert record["fixture_index"] == index
                    assert record["published_fixture_scalar"] is None
                    assert record["target_source"] == "public_point_jsonl"
                    assert record["published_q"] == fixture["published_q"]
                    assert record["recovered_fixture_scalar"] == fixture["published_fixture_scalar"]
                    assert curve.mul(record["recovered_fixture_scalar"], generator) == tuple(record["published_q"])
                checks.append({"block": block, "arm": arm, "targets_verified": 1024,
                               "walk_steps": summary["total_walk_steps"], "charges": charges})
            measurements[block][arm] = {
                "wall_seconds": item["wall_seconds"],
                "cpu_seconds": item["user_seconds"] + item["sys_seconds"],
                "rss_kib": item["max_rss_kib_linux"],
            }
    for block in range(5):
        baseline = evidence[(block, "off_a")]
        for arm in ("filter", "off_b"):
            candidate = evidence[(block, arm)]
            assert candidate[:3] == baseline[:3], (block, arm)
            for key in ("regular_states", "root_table_entries", "root_table_slots",
                        "rank", "rank_attempts", "rank_relations", "rank_failures",
                        "rank_rows_without_gain", "base_hash", "rank_probes_mean"):
                assert candidate[3][key] == baseline[3][key], (block, arm, key)
            for phase in ("index", "rank", "target"):
                for key in ("calls", "regular_calls", "exceptional_calls", "returned_pairs",
                            "no_root_returns", "consumed_candidates", "discarded_candidates",
                            "root_keys_considered", "lift_attempts"):
                    assert (candidate[3][f"{phase}_s3_counts"][key]
                            == baseline[3][f"{phase}_s3_counts"][key]), (block, arm, phase, key)
    for block in range(5):
        off_a = measurements[block]["off_a"]
        off_b = measurements[block]["off_b"]
        measurements[block]["off_geo"] = {
            key: math.sqrt(off_a[key] * off_b[key])
            for key in ("cpu_seconds", "wall_seconds")}
    paired_ratios = {
        "off_b_to_off_a": paired(measurements, "off_b", "off_a"),
        "filter_to_off_geo": paired(measurements, "filter", "off_geo"),
        "filter_to_rho": paired(measurements, "filter", "rho_normal"),
        "off_geo_to_rho": paired(measurements, "off_geo", "rho_normal"),
    }
    aa = paired_ratios["off_b_to_off_a"]
    aa_valid = (0.90 <= aa["cpu_ratio_median"] <= 1.10
                and aa["cpu_ratio_95pct_log_t_interval"][0] <= 1
                <= aa["cpu_ratio_95pct_log_t_interval"][1])
    uncontended = isolation["contended_samples"] == 0
    timing_eligible = aa_valid and uncontended
    engineering_win = (timing_eligible
                       and paired_ratios["filter_to_off_geo"]["cpu_ratio_95pct_log_t_interval"][1] < 1)
    timing_crossover_cell = (timing_eligible
                            and paired_ratios["filter_to_rho"]["cpu_ratio_95pct_log_t_interval"][1] < 1)
    return {
        "status": "PASS", "n": n, "L": 1024, "k": k,
        "targets_verified_per_arm": 5 * 1024, "rank_runs_verified": 5 * 3,
        "checks": checks, "measurements": measurements,
        "paired": paired_ratios, "aa_timing_valid": aa_valid,
        "uncontended": uncontended, "timing_eligible": timing_eligible,
        "engineering_win": engineering_win, "timing_crossover_cell": timing_crossover_cell,
        "S_off": None, "S_filter": None, "S_rho_normal": None,
        "evidence_class": "cold_full_process_diagnostic_no_common_operation_unit",
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
