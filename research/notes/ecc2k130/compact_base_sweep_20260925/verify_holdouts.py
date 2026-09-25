#!/usr/bin/env python3
"""Replay frozen point-only compact queries and identical-Q rho receipts."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
VERIFY = HERE / "verify.py"
SPEC = HERE / "input_spec.json"


def load_verify():
    spec = importlib.util.spec_from_file_location("sparse_base_verify", VERIFY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_holdout(training: Path, holdout: Path) -> dict:
    v = load_verify()
    spec = json.loads(SPEC.read_bytes())
    train = v.verify_training(training)
    train_manifest, _, train_observation = v.read_run(training)
    manifest, receipt, observed = v.read_run(holdout)
    n, r = train["n"], train["R"]
    assert manifest["mode"] == "holdout" and manifest["n"] == n and manifest["R"] == r
    assert train_manifest["n"] == n and train_manifest["R"] == r
    assert observed["n"] == n and observed["a"] == 0
    assert observed["orbit_columns"] == r
    header = observed["compact_orbit_base_header"]
    assert header == train_observation["compact_orbit_base_header"]
    math = v.load_reference_math()
    curve, generator, reps, by_point, by_x, _, point_set_sha = v.verify_base(header, math)
    assert point_set_sha == train["point_set_sha256"]
    points_bytes = (holdout / "target_points.jsonl").read_bytes()
    assert v.digest(points_bytes) == manifest["input_sha256"]
    targets = [tuple(json.loads(line)) for line in points_bytes.splitlines()]
    expected = spec["holdouts"][str(n)]
    assert [list(point) for point in targets] == [x["q"] for x in expected]
    assert len(set(targets)) == len(targets) == 3
    batch = observed["compact_orbit_point_batch"]
    assert batch["targets_requested"] == 3
    assert not batch["sat_verification_included"]
    queries = batch["query_observations"]
    assert [tuple(x["target_point"]) for x in queries] == targets
    assert [tuple(x) for x in batch["failed_target_points"]] == [targets[i] for i,x in enumerate(queries) if not x["hit"]]
    relations = batch["relations"]
    assert [tuple(x["target_point"]) for x in relations] == [targets[i] for i,x in enumerate(queries) if x["hit"]]
    failed_s3 = 2*n*batch["regular_states"]
    for query in queries:
        assert query["partner_roots"] >= query["indexed_partner_hits"]
        assert query["indexed_partner_hits"] == query["group_lift_attempts"]
        if query["hit"]:
            assert query["s3_calls"] <= failed_s3
        else:
            assert query["s3_calls"] == failed_s3
    recovered = []
    logs = train["factor_base_log_solution"]
    for relation in relations:
        point = tuple(relation["target_point"])
        index = targets.index(point)
        assert curve.is_on_curve(point) and curve.scalar(point, train["q"]) is None
        lift = v.check_witness(curve, by_x, point, relation["x_codes"], relation["pinned_intermediates"])
        record = {"seed": expected[index]["seed"], "target_point": list(point),
                  "group_lift_verified": True, "scalar_recovered": None}
        if logs is not None:
            scalar = sum(by_point[p][1] * logs[by_point[p][0]] for p in lift) % train["q"]
            assert curve.scalar(generator, scalar) == point
            assert scalar == expected[index]["scalar_validator_only"]
            record["scalar_recovered"] = scalar
        recovered.append(record)
    assert len(reps) == r
    return {"schema_version": "1.0", "classification": "FROZEN_POINT_ONLY_HOLDOUT",
            "n": n, "R": r, "training_count": train["count"],
            "training_rank": train["rank"], "point_set_sha256": point_set_sha,
            "targets": 3, "relations": len(relations),
            "scalar_logs_recovered": sum(x["scalar_recovered"] is not None for x in recovered),
            "failed_query_s3_calls_each": failed_s3,
            "total_s3_calls": sum(x["s3_calls"] for x in queries),
            "total_partner_roots": sum(x["partner_roots"] for x in queries),
            "total_indexed_partner_hits": sum(x["indexed_partner_hits"] for x in queries),
            "process_wall_ms": receipt["wall_ms"],
            "process_cpu_s": receipt["user_cpu_s"]+receipt["system_cpu_s"],
            "peak_rss_bytes": receipt["peak_rss_bytes"],
            "producer_stdout_sha256": receipt["stdout_sha256"], "rows": recovered}


def check_rho(rho_run: Path) -> dict:
    v = load_verify()
    manifest, receipt, observed = v.read_run(rho_run)
    spec = json.loads(SPEC.read_bytes())
    assert manifest["mode"] == "rho"
    n = manifest["n"]
    expected = spec["holdouts"][str(n)][manifest["seed_index"]]
    assert observed["public_hash_seed"] == expected["seed"]
    assert observed["published_q"] == expected["q"]
    assert observed["recovered_fixture_scalar"] == expected["scalar_validator_only"]
    assert observed["automorphism_size"] == 2*n
    assert observed["verified"] and observed["reference_group_validation"]
    math = v.load_reference_math()
    curve = math.Curve(n, 0, observed["field_modulus_low_terms"])
    generator = math.pt(observed["generator"])
    target = math.pt(observed["published_q"])
    order = int(observed["subgroup_order"])
    assert curve.is_on_curve(generator) and curve.is_on_curve(target)
    assert curve.scalar(generator, order) is None and curve.scalar(target, order) is None
    assert curve.scalar(generator, observed["recovered_fixture_scalar"]) == target
    return {"schema_version": "1.0", "classification": "SAME_Q_L1_RHO_CONTROL",
            "n": n, "seed": expected["seed"], "q": list(target),
            "recovered_scalar": observed["recovered_fixture_scalar"],
            "walk_steps": observed["walk_steps"],
            "ideal_steps": observed["ideal_steps"],
            "process_wall_ms": receipt["wall_ms"],
            "process_cpu_s": receipt["user_cpu_s"]+receipt["system_cpu_s"],
            "peak_rss_bytes": receipt["peak_rss_bytes"],
            "producer_stdout_sha256": receipt["stdout_sha256"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--holdout", type=Path)
    group.add_argument("--rho", type=Path)
    parser.add_argument("--training", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.holdout:
        assert args.training is not None
        report = check_holdout(args.training, args.holdout)
    else:
        assert args.training is None
        report = check_rho(args.rho)
    data = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(data)
    print(json.dumps({k: report[k] for k in ("classification", "n", "process_wall_ms")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
