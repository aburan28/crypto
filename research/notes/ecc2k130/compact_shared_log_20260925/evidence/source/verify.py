#!/usr/bin/env python3
"""Independent group-law/rank replay of the frozen shared-log control."""
from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import types

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SPEC = HERE / "input_spec.json"
PR747 = (HERE / "evidence/source/pr747_verify.py" if
         (HERE / "evidence/source/pr747_verify.py").exists() else
         HERE.parent / "compact_base_sweep_20260925/verify.py")
PR747_SHA = "dd334ff96c88c5489884c395ca81b415541547ec7673285e03c817334d46de30"
PR737 = (HERE / "evidence/source/pr737_verify.py" if
         (HERE / "evidence/source/pr737_verify.py").exists() else
         HERE.parent / "paired_fullrank_20260925/verify.py")
PR737_SHA = "ed804a6bdace125cc41e514d93b1102763581f4feeaad60ba277cace592230bf"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_module(path: Path, expected_sha: str, name: str):
    data = path.read_bytes()
    assert sha(data) == expected_sha, path
    module = types.ModuleType(name)
    module.__file__ = str(path)
    exec(compile(data, str(path), "exec"), module.__dict__)
    return module


def reference():
    return (load_module(PR747, PR747_SHA, "pr747_compact_verify"),
            load_module(PR737, PR737_SHA, "pr737_group_law"))


def read_stdout(run: Path) -> bytes:
    raw = run / "producer.stdout.jsonl"
    if raw.exists():
        return raw.read_bytes()
    with gzip.open(run / "producer.stdout.jsonl.gz", "rb") as stream:
        return stream.read()


def read_run(run: Path) -> tuple[dict, dict, list[dict]]:
    manifest_data = (run / "manifest.json").read_bytes()
    manifest = json.loads(manifest_data)
    receipt = json.loads((run / "receipt.json").read_bytes())
    assert receipt["returncode"] == 0 and receipt["termination"] is None
    assert sha(manifest_data) == receipt["manifest_sha256"]
    stdout = read_stdout(run)
    assert sha(stdout) == receipt["stdout_sha256"]
    assert sha((run / "producer.stderr.txt").read_bytes()) == receipt["stderr_sha256"]
    spec = json.loads(SPEC.read_bytes())
    assert manifest["input_spec_sha256"] == sha(SPEC.read_bytes())
    source_path = manifest["source_path"]
    assert manifest["source_sha256"] == spec["source_sha256"][source_path]
    assert sha((run / manifest["input_file"]).read_bytes()) == manifest["input_sha256"]
    snapshot = HERE / "evidence/source" / (Path(source_path).name + ".gz")
    if snapshot.exists():
        with gzip.open(snapshot, "rb") as stream:
            source_data = stream.read()
        assert sha(source_data) == manifest["source_sha256"]
    else:
        assert sha((REPO / source_path).read_bytes()) == manifest["source_sha256"]
    rows = [json.loads(line) for line in stdout.splitlines() if line.strip()]
    return manifest, receipt, rows


def stats(batch: dict, n: int, r: int) -> dict:
    regular = int(batch["regular_states"])
    assert 0 < regular <= n*r*r
    entries = int(batch["index_entries"])
    assert 0 < entries <= 2*regular
    queries = batch["query_observations"]
    failed_s3 = 2*n*regular
    for query in queries:
        assert 0 <= query["s3_calls"] <= failed_s3
        assert query["indexed_partner_hits"] <= query["partner_roots"] <= 2*query["s3_calls"]
        assert query["indexed_partner_hits"] == query["group_lift_attempts"]
        if not query["hit"]:
            assert query["s3_calls"] == failed_s3
    return {"regular_states": regular, "exceptional_states": n*r*r-regular,
            "root_candidates": 2*regular, "unique_indexed_roots": entries,
            "root_occupancy": entries/(2*regular), "failed_query_s3_calls": failed_s3,
            "total_s3_calls": sum(row["s3_calls"] for row in queries),
            "total_partner_roots": sum(row["partner_roots"] for row in queries),
            "total_indexed_partner_hits": sum(row["indexed_partner_hits"] for row in queries),
            "total_group_lift_attempts": sum(row["group_lift_attempts"] for row in queries)}


def training(run: Path) -> dict:
    prior, math = reference()
    spec = json.loads(SPEC.read_bytes())
    manifest, receipt, rows = read_run(run)
    assert manifest["mode"] == "train" and manifest["block"] is None
    assert len(rows) == 1
    observed = rows[0]
    n, r = manifest["n"], manifest["R"]
    arm = spec["arms"][str(n)]
    assert r == arm["R"] and observed["orbit_columns"] == r
    expected = arm["training_scalars_validator_only"]
    assert len(expected) == 128 and len(set(expected)) == 128
    assert [int(x) for x in (run / "target_scalars.txt").read_text().splitlines()] == expected
    header = observed["compact_orbit_base_header"]
    curve, generator, reps, by_point, by_x, _, point_set_sha = prior.verify_base(header, math)
    q = int(header["subgroup_order"])
    assert q == arm["q"] and header["orbit_columns"] == r
    assert [int(x) for x in header["generator"]] == arm["generator"]
    assert header["field_modulus_low_terms"] == arm["field_modulus_low_terms"]
    batch = observed["compact_orbit_batch"]
    assert observed["models_examined"] == observed["sat_solver_invocations"] == 0
    assert batch["targets_requested"] == 128 and not batch["sat_verification_included"]
    observations = batch["query_observations"]
    assert [row["scalar"] for row in observations] == expected
    relations = batch["relations"]
    assert [row["scalar"] for row in relations] == [row["scalar"] for row in observations if row["hit"]]
    assert batch["failed_target_scalars"] == [row["scalar"] for row in observations if not row["hit"]]
    assert batch["targets_extracted"] == len(relations)
    counters = stats(batch, n, r)
    incremental = math.IncrementalRank(r, q)
    matrix, rhs = [], []
    full_rank_at = None
    frozen_solution = None
    post_rank_predictions = 0
    for index, relation in enumerate(relations, start=1):
        scalar = int(relation["scalar"])
        target = curve.scalar(generator, scalar)
        lift = prior.check_witness(curve, by_x, target,
                                   relation["x_codes"], relation["pinned_intermediates"])
        row = [0]*r
        for point in lift:
            column, coefficient = by_point[point]
            row[column] = (row[column]+coefficient) % q
        if frozen_solution is not None:
            assert sum(x*y for x, y in zip(row, frozen_solution)) % q == scalar
            post_rank_predictions += 1
        matrix.append(row)
        rhs.append(scalar)
        rank = incremental.add(row)
        if rank == r and full_rank_at is None:
            full_rank_at = index
            solved_rank, frozen_solution = math.rank_and_solve(matrix, rhs, q, r)
            assert solved_rank == r and frozen_solution is not None
    rank, solution = math.rank_and_solve(matrix, rhs, q, r)
    assert rank == len(incremental.pivots) and solution == frozen_solution
    if solution is not None:
        for rep, log in zip(reps, solution):
            assert curve.scalar(generator, log) == rep
    return {"schema_version": "1.0", "classification": "FULL_RANK" if rank == r else "RANK_CENSORED",
            "n": n, "R": r, "F": arm["F"], "q": q,
            "run_manifest_sha256": sha((run / "manifest.json").read_bytes()),
            "producer_stdout_sha256": receipt["stdout_sha256"],
            "point_set_sha256": point_set_sha, "training_queries": 128,
            "relations": len(relations), "rank": rank, "full_rank_at_relation": full_rank_at,
            "post_rank_scalar_predictions": post_rank_predictions,
            "factor_base_log_solution": solution,
            "factor_base_log_solution_verified": solution is not None,
            "producer_base_ms": observed["base_ms"],
            "producer_scan_ms": batch["regular_state_scan_ms"],
            "producer_index_ms": batch["root_index_build_ms"],
            "producer_query_ms": batch["batch_loop_wall_ms"],
            "process_wall_ms": receipt["wall_ms"],
            "process_cpu_s": receipt["user_cpu_s"]+receipt["system_cpu_s"],
            "peak_rss_bytes": receipt["peak_rss_bytes"], **counters}


def compact(train: Path, run: Path, train_report: dict) -> dict:
    prior, math = reference()
    spec = json.loads(SPEC.read_bytes())
    train_manifest, _, train_rows = read_run(train)
    assert train_report["run_manifest_sha256"] == sha((train / "manifest.json").read_bytes())
    assert train_report["factor_base_log_solution_verified"]
    manifest, receipt, rows = read_run(run)
    assert manifest["mode"] == "compact" and len(rows) == 1
    n, r, block, length = manifest["n"], manifest["R"], manifest["block"], manifest["length"]
    arm = spec["arms"][str(n)]
    assert n == train_report["n"] == train_manifest["n"]
    assert r == train_report["R"] == train_manifest["R"]
    assert manifest["input_sha256"] == arm["blocks"][block]["files"][str(length)]["sha256"]
    expected = arm["blocks"][block]
    targets = [tuple(json.loads(line)) for line in (run / "target_points.jsonl").read_bytes().splitlines()]
    assert [list(point) for point in targets] == expected["points"][:length]
    assert len(set(targets)) == length
    observed = rows[0]
    assert observed["n"] == n and observed["a"] == 0 and observed["orbit_columns"] == r
    header = observed["compact_orbit_base_header"]
    assert header == train_rows[0]["compact_orbit_base_header"]
    curve, generator, _, by_point, by_x, _, point_set_sha = prior.verify_base(header, math)
    assert point_set_sha == train_report["point_set_sha256"]
    batch = observed["compact_orbit_point_batch"]
    assert batch["targets_requested"] == length and not batch["sat_verification_included"]
    observations = batch["query_observations"]
    assert [tuple(row["target_point"]) for row in observations] == targets
    relations = batch["relations"]
    assert [tuple(row["target_point"]) for row in relations] == [targets[i] for i,row in enumerate(observations) if row["hit"]]
    assert [tuple(point) for point in batch["failed_target_points"]] == [targets[i] for i,row in enumerate(observations) if not row["hit"]]
    assert batch["targets_extracted"] == len(relations)
    counters = stats(batch, n, r)
    logs = train_report["factor_base_log_solution"]
    recovered = []
    for relation in relations:
        point = tuple(relation["target_point"])
        index = targets.index(point)
        assert curve.is_on_curve(point) and curve.scalar(point, arm["q"]) is None
        lift = prior.check_witness(curve, by_x, point,
                                   relation["x_codes"], relation["pinned_intermediates"])
        scalar = sum(by_point[p][1]*logs[by_point[p][0]] for p in lift) % arm["q"]
        assert curve.scalar(generator, scalar) == point
        assert scalar == expected["validator_scalars_only"][index]
        recovered.append({"index": index, "target_point": list(point), "scalar": scalar})
    return {"schema_version": "1.0", "classification": "COMPLETE" if len(recovered) == length else "QUERY_CENSORED",
            "n": n, "R": r, "block": block, "length": length,
            "run_manifest_sha256": sha((run / "manifest.json").read_bytes()),
            "producer_stdout_sha256": receipt["stdout_sha256"],
            "point_set_sha256": point_set_sha, "relations": len(relations),
            "recovered": recovered,
            "producer_base_ms": observed["base_ms"],
            "producer_scan_ms": batch["regular_state_scan_ms"],
            "producer_index_ms": batch["root_index_build_ms"],
            "producer_query_ms": batch["batch_loop_wall_ms"],
            "process_wall_ms": receipt["wall_ms"],
            "process_cpu_s": receipt["user_cpu_s"]+receipt["system_cpu_s"],
            "peak_rss_bytes": receipt["peak_rss_bytes"], **counters}


def rho(run: Path) -> dict:
    _, math = reference()
    spec = json.loads(SPEC.read_bytes())
    manifest, receipt, rows = read_run(run)
    assert manifest["mode"] == "rho"
    n, block, length = manifest["n"], manifest["block"], manifest["length"]
    arm = spec["arms"][str(n)]
    expected = arm["blocks"][block]
    assert manifest["input_sha256"] == expected["files"][str(length)]["sha256"]
    targets = [tuple(json.loads(line)) for line in (run / "target_points.jsonl").read_bytes().splitlines()]
    assert [list(point) for point in targets] == expected["points"][:length]
    assert len(rows) == length+1
    summary = rows[-1]
    assert summary["kind"] == "rho_ks_batch_summary"
    assert summary["n"] == n and summary["a"] == 0
    assert summary["fixtures"] == length and summary["automorphism_size"] == 2*n
    assert summary["quotient_mode"] == "signed_frobenius" and summary["dp_bits"] == 4
    assert summary["precompute_walks"] == 0 and summary["precompute_steps"] == 0
    assert summary["target_source"] == "explicit_public_points" and summary["all_verified"]
    assert summary["batch_seed"] == 202609250000 + 100*n + block
    assert summary["corpus"] == f"sparse-shared-n{n}-b{block}-L{length}"
    curve = math.Curve(n, 0, arm["field_modulus_low_terms"])
    generator = math.pt(arm["generator"])
    assert curve.scalar(generator, arm["q"]) is None
    cross, steps, last_table = 0, 0, 0
    recovered = []
    for index, row in enumerate(rows[:-1]):
        assert row["kind"] == "rho_ks_batch_fixture"
        assert row["n"] == n and row["a"] == 0 and row["fixture_index"] == index
        assert row["batch_seed"] == summary["batch_seed"]
        assert row["automorphism_size"] == 2*n and row["quotient_mode"] == "signed_frobenius"
        assert row["target_source"] == "explicit_public_points"
        assert row["published_fixture_scalar"] is None and row["verified"]
        point = targets[index]
        scalar = int(row["recovered_fixture_scalar"])
        assert tuple(row["published_q"]) == point
        assert curve.is_on_curve(point) and curve.scalar(point, arm["q"]) is None
        assert curve.scalar(generator, scalar) == point
        assert scalar == expected["validator_scalars_only"][index]
        assert row["table_entries_before"] >= last_table
        assert row["table_entries_after"] >= row["table_entries_before"]
        last_table = row["table_entries_after"]
        cross += int(row["cross_target_solve"])
        steps += row["walk_steps"]
        recovered.append(scalar)
    assert summary["cross_target_solves"] == cross
    assert summary["total_walk_steps"] == steps
    assert summary["table_entries"] == last_table
    assert summary["charges"]["group_additions"] >= steps
    return {"schema_version": "1.0", "classification": "COMPLETE",
            "n": n, "block": block, "length": length,
            "run_manifest_sha256": sha((run / "manifest.json").read_bytes()),
            "producer_stdout_sha256": receipt["stdout_sha256"],
            "recovered_scalars": recovered,
            "total_walk_steps": steps, "cross_target_solves": cross,
            "table_entries": last_table, "charges": summary["charges"],
            "producer_setup_ms": summary["setup_ms"],
            "producer_in_process_ms": summary["in_process_ms"],
            "process_wall_ms": receipt["wall_ms"],
            "process_cpu_s": receipt["user_cpu_s"]+receipt["system_cpu_s"],
            "peak_rss_bytes": receipt["peak_rss_bytes"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("train", "compact", "rho"), required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--training", type=Path)
    parser.add_argument("--training-report", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.mode == "train":
        assert args.training is None and args.training_report is None
        report = training(args.run)
    elif args.mode == "compact":
        assert args.training and args.training_report
        report = compact(args.training, args.run, json.loads(args.training_report.read_bytes()))
    else:
        assert args.training is None and args.training_report is None
        report = rho(args.run)
    data = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(data)
    print(json.dumps({key: report[key] for key in ("classification", "n", "process_wall_ms")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
