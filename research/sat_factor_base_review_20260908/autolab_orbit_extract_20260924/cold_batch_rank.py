#!/usr/bin/env python3
"""Cold reused-index n=53 batch, independent group replay, and orbit rank.

This is a stage/rank probe. Targets are published synthetic multiples of G;
it does not constitute an unknown-target discrete-log recovery or a rho race.
"""

import argparse
from collections import defaultdict
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import subprocess
import time


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BASE_GZ = HERE / "independent_replay_20260924_codex/base_header.jsonl.gz"
REPLAY = HERE / "independent_replay_20260924_codex/replay.py"
EXE = REPO / "target/release/examples/koblitz_s5_sat_instance"
SOURCE = REPO / "examples/koblitz_s5_sat_instance.rs"
DOMAIN = b"ECC2K53-COMPACT-COLD-BATCH-20260924-v1/"
BASE_HASH = "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def load_verifier():
    spec = importlib.util.spec_from_file_location("independent_replay", REPLAY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def target_schedule(count, order):
    values = []
    used = set()
    counter = 0
    while len(values) < count:
        candidate = 1 + int.from_bytes(
            hashlib.sha256(DOMAIN + str(counter).encode()).digest(), "big"
        ) % (order - 1)
        counter += 1
        if candidate not in used:
            used.add(candidate)
            values.append(candidate)
    return values


def verify_orbit_labels(curve, header):
    points = [tuple(point) for point in header["factor_base_point_coordinates"]]
    labels = [tuple(label) for label in header["factor_base_point_labels"]]
    reps = [tuple(point) for point in header["factor_base_representatives"]]
    assert len(points) == len(labels) == len(reps) * 106
    by_point = dict(zip(points, labels))
    assert len(by_point) == len(points)
    first = reps[0]
    squared = (curve.multiply(first[0], first[0]), curve.multiply(first[1], first[1]))
    column, lam = by_point[squared]
    assert column == 0
    assert curve.scalar(curve.generator, lam) == (
        curve.multiply(curve.generator[0], curve.generator[0]),
        curve.multiply(curve.generator[1], curve.generator[1]),
    )
    assert pow(lam, curve.n, curve.order) == 1
    seen = set()
    for column, rep in enumerate(reps):
        assert curve.on_curve(rep)
        assert curve.scalar(rep, curve.order) is None
        point = rep
        coefficient = 1
        for _ in range(curve.n):
            negative = (point[0], point[0] ^ point[1])
            assert by_point[point] == (column, coefficient)
            assert by_point[negative] == (column, (-coefficient) % curve.order)
            seen.add(point)
            seen.add(negative)
            point = (curve.multiply(point[0], point[0]), curve.multiply(point[1], point[1]))
            coefficient = coefficient * lam % curve.order
        assert point == rep and coefficient == 1
    assert seen == set(points)
    return by_point, reps, lam


class Echelon:
    def __init__(self, columns, modulus):
        self.columns = columns
        self.modulus = modulus
        self.pivots = {}

    def insert(self, row, rhs):
        modulus = self.modulus
        row = {col: value % modulus for col, value in row.items() if value % modulus}
        rhs %= modulus
        while row:
            pivot = min(row)
            if pivot not in self.pivots:
                inverse = pow(row[pivot], -1, modulus)
                self.pivots[pivot] = (
                    {col: value * inverse % modulus for col, value in row.items()},
                    rhs * inverse % modulus,
                )
                return True
            old, old_rhs = self.pivots[pivot]
            factor = row[pivot]
            for col, value in old.items():
                updated = (row.get(col, 0) - factor * value) % modulus
                if updated:
                    row[col] = updated
                else:
                    row.pop(col, None)
            rhs = (rhs - factor * old_rhs) % modulus
        assert rhs == 0, "inconsistent relation rows"
        return False

    def solution(self):
        if len(self.pivots) != self.columns:
            return None
        result = [0] * self.columns
        for pivot in sorted(self.pivots, reverse=True):
            row, rhs = self.pivots[pivot]
            result[pivot] = (rhs - sum(
                value * result[col] for col, value in row.items() if col > pivot
            )) % self.modulus
        return result


def run(args):
    out = args.out
    if args.replay_only:
        assert out.is_dir()
    else:
        out.mkdir(parents=True, exist_ok=False)
    with gzip.open(BASE_GZ, "rb") as stream:
        base = stream.read()
    assert base.count(b"\n") == 1
    header = json.loads(base)
    assert header["base_hash"] == BASE_HASH
    assert header["n"] == 53 and header["a"] == 0
    assert header["orbit_columns"] == 220
    order = header["subgroup_order"]
    scalars = target_schedule(args.targets, order)
    target_bytes = "".join(f"{scalar}\n" for scalar in scalars).encode()
    if args.replay_only:
        if (out / "base_header.jsonl").exists():
            assert (out / "base_header.jsonl").read_bytes() == base
        assert (out / "target_scalars.txt").read_bytes() == target_bytes
    else:
        (out / "base_header.jsonl").write_bytes(base)
        (out / "target_scalars.txt").write_bytes(target_bytes)
    recorded_manifest = (json.loads((out / "manifest.json").read_text())
                         if args.replay_only else None)
    manifest = {
        "schema_version": "1.0", "experiment": "cold_compact_orbit_rank_probe",
        "scope": "Public synthetic n=53 known-scalar targets; no unknown-target DLP or rho comparison",
        "domain": DOMAIN.decode(), "targets": len(scalars),
        "subgroup_order": order, "factor_base_hash": BASE_HASH,
        "base_gzip_sha256": digest(BASE_GZ.read_bytes()),
        "base_header_sha256": digest(base), "target_schedule_sha256": digest(target_bytes),
        "producer_source_sha256": (recorded_manifest["producer_source_sha256"]
                                   if args.replay_only else digest(SOURCE.read_bytes())),
        "producer_executable_sha256": (recorded_manifest["producer_executable_sha256"]
                                       if args.replay_only else digest(EXE.read_bytes())),
        "verifier_sha256": digest(REPLAY.read_bytes()),
        "command": [str(EXE), "53", "0", "1", "10", "natural", "1", "2000", "1", "internal"],
        "environment": {
            "KIC_ALGEBRA_ENCODING": "orbit_factorized",
            "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
            "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
            "KIC_ORBIT_REP_ENCODING": "one_hot",
            "KIC_ORBIT_BATCH_ONLY": "1",
            "KIC_FACTOR_BASE_JSONL": "base_header.jsonl",
            "KIC_ORBIT_TARGET_SCALARS": "target_scalars.txt",
        },
        "stopping_rule": "Run exactly the first N schedule targets; count every failure and relation rank",
    }
    if args.replay_only:
        recorded = json.loads((out / "manifest.json").read_text())
        for key in ("schema_version", "experiment", "domain", "targets", "subgroup_order",
                    "factor_base_hash", "base_gzip_sha256", "base_header_sha256",
                    "target_schedule_sha256", "verifier_sha256", "environment"):
            assert recorded[key] == manifest[key], key
        manifest = recorded
    else:
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    env = {key: val for key, val in os.environ.items() if not key.startswith("KIC_")}
    env.update(manifest["environment"])
    env["KIC_FACTOR_BASE_JSONL"] = str(out / "base_header.jsonl")
    env["KIC_ORBIT_TARGET_SCALARS"] = str(out / "target_scalars.txt")
    env["KIC_TASK_ID"] = "TASK-IC-COMPACT-COLD-RANK-20260924"
    if args.replay_only:
        stdout = (out / "producer.stdout.jsonl").read_text()
        prior = json.loads((out / "validation.json").read_text())
        wall_ms = prior["process_wall_ms"]
        child_peak_rss_raw = prior["child_peak_rss_raw"]
    else:
        started = time.perf_counter()
        process = subprocess.run(manifest["command"], cwd=REPO, env=env,
                                 capture_output=True, text=True, timeout=args.timeout)
        wall_ms = (time.perf_counter() - started) * 1000
        stdout = process.stdout
        (out / "producer.stdout.jsonl").write_text(stdout)
        (out / "producer.stderr.txt").write_text(process.stderr)
        assert process.returncode == 0, (process.returncode, process.stderr[-1000:])
        child_peak_rss_raw = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    assert len(stdout.splitlines()) == 1
    observation = json.loads(stdout)
    batch = observation["compact_orbit_batch"]
    assert observation["factor_base_input_hash"] == BASE_HASH
    assert batch["targets_requested"] == len(scalars)
    assert batch["targets_extracted"] + len(batch["failed_target_scalars"]) == len(scalars)
    assert observation["models_examined"] == 0
    assert observation["decomposition_verdict"] == "UNKNOWN"
    # A batch-only run never invokes SAT. Older producer receipts falsely
    # claimed it did; retain those raw receipts, but interpret them here.
    assert batch["sat_verification_included"] in (True, False)
    verifier = load_verifier()
    curve = verifier.Curve(header)
    by_point, reps, lam = verify_orbit_labels(curve, header)
    by_x = defaultdict(list)
    for point in by_point:
        by_x[point[0]].append(point)
    matrix = Echelon(len(reps), order)
    rank_gains = []
    first_full_rank_at_extracted = None
    full_rank_solution = None
    heldout_predictions_verified = 0
    for relation_index, relation in enumerate(batch["relations"], start=1):
        scalar = relation["scalar"]
        target = curve.scalar(curve.generator, scalar)
        assert target is not None
        lift = verifier.check_witness(curve, by_x, target,
                                      relation["x_codes"], relation["pinned_intermediates"])
        row = {}
        for point in map(tuple, lift):
            column, coefficient = by_point[point]
            row[column] = (row.get(column, 0) + coefficient) % order
        if full_rank_solution is not None:
            prediction = sum(coefficient * full_rank_solution[column]
                             for column, coefficient in row.items()) % order
            assert prediction == scalar, "post-rank held-out scalar prediction failed"
            heldout_predictions_verified += 1
        if matrix.insert(row, scalar):
            rank_gains.append(scalar)
            if len(matrix.pivots) == len(reps) and first_full_rank_at_extracted is None:
                first_full_rank_at_extracted = relation_index
                full_rank_solution = matrix.solution()
    solution = matrix.solution()
    assert solution == full_rank_solution
    if solution is not None:
        for rep, log in zip(reps, solution):
            assert curve.scalar(curve.generator, log) == rep
    report = {
        "schema_version": "1.0", "classification": "STAGE_AND_RANK_ONLY",
        "targets_requested": len(scalars),
        "targets_extracted": len(batch["relations"]),
        "failed_target_scalars": batch["failed_target_scalars"],
        "rank": len(matrix.pivots), "columns": len(reps),
        "rank_gain_scalars": rank_gains,
        "full_rank_at_extracted": first_full_rank_at_extracted,
        "independent_holdout_relations_after_full_rank": heldout_predictions_verified,
        "all_orbit_labels_independently_verified": True,
        "all_relations_independently_group_verified": True,
        "factor_base_log_solution_verified": solution is not None,
        "factor_base_log_solution": solution,
        "frobenius_eigenvalue": lam,
        "process_wall_ms": wall_ms,
        "child_peak_rss_raw": child_peak_rss_raw,
        "producer_charged_stage_ms": batch["charged_total_ms"],
        "producer_scan_ms": batch["regular_state_scan_ms"],
        "producer_root_index_ms": batch["root_index_build_ms"],
        "producer_query_ms_sum": batch["query_ms_sum"],
        "sat_verification_actually_included": False,
        "producer_stdout_sha256": digest(stdout.encode()),
        "manifest_sha256": digest((out / "manifest.json").read_bytes()),
    }
    (out / "validation.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: report[key] for key in (
        "targets_requested", "targets_extracted", "rank", "columns",
        "factor_base_log_solution_verified", "process_wall_ms", "producer_charged_stage_ms",
        "child_peak_rss_raw")}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--targets", type=int, default=512)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--replay-only", action="store_true")
    arguments = parser.parse_args()
    assert 1 <= arguments.targets <= 1024
    run(arguments)
