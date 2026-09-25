#!/usr/bin/env python3
"""Independent affine/group/rank replay of a compact-base sweep run."""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
PR737_VERIFY = (HERE / "evidence/source/pr737_independent_math.py"
                if (HERE / "evidence/source/pr737_independent_math.py").exists()
                else HERE.parent / "paired_fullrank_20260925/verify.py")
SPEC = HERE / "input_spec.json"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_reference_math():
    spec = importlib.util.spec_from_file_location("paired_fullrank_verify", PR737_VERIFY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_stdout(run: Path) -> bytes:
    raw = run / "producer.stdout.jsonl"
    if raw.exists():
        return raw.read_bytes()
    with gzip.open(run / "producer.stdout.jsonl.gz", "rb") as stream:
        return stream.read()


def read_run(run: Path) -> tuple[dict, dict, dict]:
    manifest = json.loads((run / "manifest.json").read_bytes())
    receipt = json.loads((run / "receipt.json").read_bytes())
    stdout = read_stdout(run)
    assert receipt["returncode"] == 0 and not receipt["timed_out"]
    assert digest(stdout) == receipt["stdout_sha256"]
    assert digest((run / "producer.stderr.txt").read_bytes()) == receipt["stderr_sha256"]
    assert digest((run / "manifest.json").read_bytes()) == receipt["manifest_sha256"]
    assert manifest["input_spec_sha256"] == digest(SPEC.read_bytes())
    source_snapshot = HERE / "evidence/source" / (Path(manifest["source_path"]).name + ".gz")
    if source_snapshot.exists():
        with gzip.open(source_snapshot, "rb") as stream:
            source_bytes = stream.read()
    else:
        source_bytes = (REPO / manifest["source_path"]).read_bytes()
    assert digest(source_bytes) == manifest["source_sha256"]
    assert len(stdout.splitlines()) == 1
    return manifest, receipt, json.loads(stdout)


def check_witness(curve, points_by_x, target, codes, intermediates):
    assert len(codes) == 4 and len(intermediates) == 2
    assert all(code in points_by_x for code in codes)
    u, v = intermediates

    def s3(x, y, z):
        cross = curve.mul(x, y) ^ curve.mul(x, z) ^ curve.mul(y, z)
        return curve.square(cross) ^ curve.mul(curve.mul(x, y), z) ^ 1

    assert s3(codes[0], codes[1], u) == 0
    assert s3(codes[2], codes[3], v) == 0
    assert s3(u, v, target[0]) == 0
    for choice in itertools.product(*(points_by_x[x] for x in codes)):
        left, right = curve.add(choice[0], choice[1]), curve.add(choice[2], choice[3])
        if left is not None and right is not None and left[0] == u and right[0] == v:
            if curve.add(left, right) == target:
                return choice
    raise AssertionError("no independently valid point lift for S3 witness")


def verify_base(header: dict, math):
    n, a, order = int(header["n"]), int(header["a"]), int(header["subgroup_order"])
    curve = math.Curve(n, a, header["field_modulus_low_terms"])
    assert a == 0 and header["kind"] == "point_defined_factor_base"
    assert header["signed_automorphism_size"] == 2 * n
    assert header["factor_base_points"] == 2 * n * header["orbit_columns"]
    generator = math.pt(header["generator"])
    assert curve.is_on_curve(generator) and curve.scalar(generator, order) is None
    points = [math.pt(raw) for raw in header["factor_base_point_coordinates"]]
    labels = [tuple(raw) for raw in header["factor_base_point_labels"]]
    reps = [math.pt(raw) for raw in header["factor_base_representatives"]]
    assert len(points) == len(labels) == header["factor_base_points"]
    assert len(reps) == header["orbit_columns"]
    assert len(set(points)) == len(points)
    assert all(curve.is_on_curve(p) for p in points)
    assert all(curve.is_on_curve(p) and curve.scalar(p, order) is None for p in reps)
    by_point = dict(zip(points, labels))
    first_frob = (curve.square(reps[0][0]), curve.square(reps[0][1]))
    assert first_frob in by_point and by_point[first_frob][0] == 0
    lam = by_point[first_frob][1]
    assert curve.scalar(generator, lam) == (curve.square(generator[0]), curve.square(generator[1]))
    seen = set()
    for column, rep in enumerate(reps):
        current, coefficient = rep, 1
        for _ in range(n):
            negative = (current[0], current[1] ^ current[0])
            assert by_point[current] == (column, coefficient)
            assert by_point[negative] == (column, (-coefficient) % order)
            seen.add(current)
            seen.add(negative)
            current = (curve.square(current[0]), curve.square(current[1]))
            coefficient = coefficient * lam % order
        assert current == rep and coefficient == 1
    assert seen == set(points)
    points_by_x = defaultdict(list)
    for point in points:
        points_by_x[point[0]].append(point)
    assert all(len(group) == 2 and group[0][1] ^ group[1][1] == x
               for x, group in points_by_x.items())
    point_set_sha = digest(json.dumps(header["factor_base_point_coordinates"],
                                      separators=(",", ":")).encode())
    return curve, generator, reps, by_point, points_by_x, lam, point_set_sha


def verify_training(run: Path) -> dict:
    math = load_reference_math()
    spec = json.loads(SPEC.read_bytes())
    manifest, receipt, observed = read_run(run)
    assert manifest["mode"] == "train"
    n, r, count = manifest["n"], manifest["R"], manifest["count"]
    assert manifest["input_sha256"] == digest((run / "target_scalars.txt").read_bytes())
    expected_scalars = [int(x) for x in (HERE / f"training_scalars_n{n}.txt").read_text().splitlines()[:count]]
    scalars = [int(x) for x in (run / "target_scalars.txt").read_text().splitlines()]
    assert scalars == expected_scalars and len(set(scalars)) == count
    header = observed["compact_orbit_base_header"]
    curve, generator, reps, by_point, by_x, lam, point_set_sha = verify_base(header, math)
    order = int(header["subgroup_order"])
    assert order == spec["arms"][str(n)][0]["subgroup_order"]
    assert len(reps) == r == observed["orbit_columns"]
    assert header["factor_base_points"] == observed["factor_base_points"]
    batch = observed["compact_orbit_batch"]
    assert observed["models_examined"] == observed["sat_solver_invocations"] == 0
    assert batch["sat_verification_included"] is False
    assert batch["targets_requested"] == count
    observations = batch["query_observations"]
    assert len(observations) == count
    assert [x["scalar"] for x in observations] == scalars
    misses = [x["scalar"] for x in observations if not x["hit"]]
    assert misses == batch["failed_target_scalars"]
    relations = batch["relations"]
    assert [x["scalar"] for x in relations] == [x["scalar"] for x in observations if x["hit"]]
    assert len(relations) + len(misses) == count
    regular = int(batch["regular_states"])
    candidate_roots = 2 * regular
    indexed_roots = int(batch["index_entries"])
    assert 0 < indexed_roots <= candidate_roots
    assert regular <= n * r * r
    failed_s3 = 2 * n * regular
    for x in observations:
        assert x["partner_roots"] >= x["indexed_partner_hits"]
        assert x["indexed_partner_hits"] == x["group_lift_attempts"]
        assert 0 <= x["s3_calls"] <= failed_s3
        if x["hit"]:
            assert x["group_lift_attempts"] > 0
        else:
            assert x["s3_calls"] == failed_s3
    incremental = math.IncrementalRank(r, order)
    rows, rhs = [], []
    first_full_rank = None
    heldout_predictions = 0
    frozen_solution = None
    for index, rel in enumerate(relations, start=1):
        scalar = rel["scalar"]
        target = curve.scalar(generator, scalar)
        assert target is not None and curve.is_on_curve(target)
        lift = check_witness(curve, by_x, target, rel["x_codes"], rel["pinned_intermediates"])
        row = [0] * r
        for point in lift:
            column, coefficient = by_point[point]
            row[column] = (row[column] + coefficient) % order
        if frozen_solution is not None:
            assert sum(x*y for x, y in zip(row, frozen_solution)) % order == scalar
            heldout_predictions += 1
        rows.append(row)
        rhs.append(scalar)
        rank = incremental.add(row)
        if rank == r and first_full_rank is None:
            first_full_rank = index
            solved_rank, frozen_solution = math.rank_and_solve(rows, rhs, order, r)
            assert solved_rank == r and frozen_solution is not None
    rank, solution = math.rank_and_solve(rows, rhs, order, r)
    assert rank == len(incremental.pivots)
    assert solution == frozen_solution
    if solution is not None:
        for rep, log in zip(reps, solution):
            assert curve.scalar(generator, log) == rep
    report = {"schema_version": "1.0", "classification": "SPARSE_COMPACT_BASE_DIAGNOSTIC",
              "run_manifest_sha256": digest((run / "manifest.json").read_bytes()),
              "producer_stdout_sha256": receipt["stdout_sha256"],
              "n": n, "R": r, "F": 2*n*r, "q": order, "count": count,
              "hits": len(relations), "misses": len(misses),
              "rank": rank, "full_rank_at_relation": first_full_rank,
              "post_rank_scalar_predictions": heldout_predictions,
              "factor_base_log_solution": solution,
              "factor_base_log_solution_verified": solution is not None,
              "point_set_sha256": point_set_sha, "scanned_x": header["scanned_x"],
              "regular_states": regular, "exceptional_states": n*r*r-regular,
              "root_candidates": candidate_roots, "unique_indexed_roots": indexed_roots,
              "root_occupancy": indexed_roots/candidate_roots,
              "failed_query_s3_calls_each": failed_s3,
              "total_s3_calls": sum(x["s3_calls"] for x in observations),
              "total_partner_roots": sum(x["partner_roots"] for x in observations),
              "total_indexed_partner_hits": sum(x["indexed_partner_hits"] for x in observations),
              "total_group_lift_attempts": sum(x["group_lift_attempts"] for x in observations),
              "process_wall_ms": receipt["wall_ms"],
              "process_cpu_s": receipt["user_cpu_s"]+receipt["system_cpu_s"],
              "peak_rss_bytes": receipt["peak_rss_bytes"],
              "producer_base_ms": observed["base_ms"],
              "producer_scan_ms": batch["regular_state_scan_ms"],
              "producer_index_ms": batch["root_index_build_ms"],
              "producer_target_loop_ms": batch["batch_loop_wall_ms"],
              "frobenius_eigenvalue": lam}
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = verify_training(args.training)
    data = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(data)
    print(json.dumps({key: report[key] for key in ("n", "R", "count", "hits", "rank",
                                                     "full_rank_at_relation", "root_occupancy",
                                                     "failed_query_s3_calls_each")}, sort_keys=True))


if __name__ == "__main__":
    main()
