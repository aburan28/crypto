#!/usr/bin/env sage -python
"""Independently replay the frozen N53 S3 panel through checked Sage."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

from sage.all import EllipticCurve, GF, PolynomialRing, matrix, vector
from sage.env import SAGE_VERSION

from run_panel import canonical, semantic_digest


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CRYPTO = ROOT.parent / "crypto"
BASE_FILE = HERE / "fb244_preflight.json"
ROWS_FILE = HERE / "measurement_rows_raw.jsonl"
RECEIPT_FILE = HERE / "independent_sage_replay.json"
VERIFIED_FILE = HERE / "measurement_rows_verified.jsonl"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def coefficient_row(indices: list[int], labels: list[tuple[int, int]],
                    columns: int, r: int) -> list[int]:
    row = [0] * columns
    for index in indices:
        column, coefficient = labels[index]
        row[column] = (row[column] + coefficient) % r
    return row


def streaming_span(pivots: dict[int, list[int]], target: list[int],
                   columns: int, r: int) -> int | None:
    remainder = target[:]
    value = 0
    for column in range(columns):
        factor = remainder[column]
        if factor == 0:
            continue
        pivot = pivots.get(column)
        if pivot is None:
            return None
        for j in range(column, columns):
            remainder[j] = (remainder[j] - factor * pivot[j]) % r
        value = (value + factor * pivot[columns]) % r
    return value if not any(remainder) else None


def main() -> None:
    if RECEIPT_FILE.exists() or VERIFIED_FILE.exists():
        raise FileExistsError("independent replay artifacts already exist")
    started = time.perf_counter_ns()
    rows = [json.loads(line) for line in ROWS_FILE.read_text().splitlines() if line]
    assert len(rows) == 10 and len({row["run_id"] for row in rows}) == 10
    workload = load(HERE / "workload.json")
    summary = load(HERE / "paired_summary_exploratory.json")
    freeze = load(HERE / "freeze_receipt.json")
    assert summary["measurement_rows_sha256"] == sha(ROWS_FILE)
    assert freeze["workload_sha256"] == sha(HERE / "workload.json")
    assert workload["identity_sha256"] == hashlib.sha256(
        canonical(workload["record"])
    ).hexdigest()
    assert sha(HERE / "sage_runtime_info.json")
    manifests = {name: load(HERE / f"{name}_manifest.json")
                 for name in ("baseline", "candidate")}
    for name, manifest in manifests.items():
        assert sha(HERE / f"{name}_manifest.json") == (
            freeze["variants"][name]["manifest_sha256"]
        )
        assert hashlib.sha256(canonical(manifest["identity_record"])).hexdigest() == (
            manifest["identity_sha256"]
        )
        candidate = manifest["candidate_freeze"]
        assert sha(Path(candidate["source_path"])) == candidate["source_sha256"]
        assert sha(Path(candidate["binary_path"])) == candidate["binary_sha256"]
        assert sha(Path(candidate["target_path"])) == candidate["target_sha256"]

    raw_runs = []
    for row in rows:
        manifest = manifests[row["variant"]]
        assert row["candidate_id"] == manifest["candidate_id"]
        assert row["workload_id"] == workload["workload_id"]
        assert row["run_id"] == (
            f"{row['candidate_id']}W{row['workload_id']}R{row['repetition']}"
        )
        path = Path(row["artifacts"]["raw"])
        assert sha(path) == row["artifacts"]["raw_sha256"]
        assert sha(Path(row["artifacts"]["stdout"])) == (
            row["artifacts"]["stdout_sha256"]
        )
        assert sha(Path(row["artifacts"]["stderr"])) == (
            row["artifacts"]["stderr_sha256"]
        )
        run = load(path)
        assert row["status"] == "reported_success_pending_sage"
        assert semantic_digest(run) == (
            row["semantic_digest_excluding_clocks_memory_and_s3_counters"]
        )
        assert run["group_verified"] is True
        assert run["target"] == workload["record"]["targets"][0]
        timings = run["timing_ms"]
        assert abs(row["online_ms"] - timings["target_online_phase_sum"]) < 0.02
        assert abs(sum(row["online_exclusive_phases_ms"].values())
                   - row["online_ms"]) < 0.02
        raw_runs.append(run)
    semantic_hashes = {semantic_digest(run) for run in raw_runs}
    assert len(semantic_hashes) == 1
    run = raw_runs[0]

    record = manifests["candidate"]["identity_record"]
    field = record["field"]
    curve = record["curve"]
    n = field["degree"]
    r = curve["subgroup_order"]
    assert n == 53 and field["characteristic"] == 2
    R = PolynomialRing(GF(2), "t")
    t = R.gen()
    F = GF(2**n, name="z", modulus=t**53 + t**6 + t**2 + t + 1)
    z = F.gen()
    powers = [z**i for i in range(n)]

    def from_word(value: int):
        # Build directly from the declared little-endian polynomial bits.
        value = int(value)
        assert 0 <= value < 1 << n
        result = F(0)
        while value:
            bit = (value & -value).bit_length() - 1
            result += powers[bit]
            value &= value - 1
        return result

    assert all(from_word(1 << i) == z**i for i in range(n))
    E = EllipticCurve(F, [F(1), F(curve["coefficients"]["a"]),
                          F(0), F(0), F(curve["coefficients"]["b"])])
    K = GF(r)

    def point(words: list[int]):
        return E(from_word(words[0]), from_word(words[1]))

    G = point(curve["generator"])
    Q = point(workload["record"]["targets"][0])
    assert G != E(0) and r * G == E(0)
    assert Q != E(0) and r * Q == E(0)
    lam = record["endomorphism"]["frobenius_eigenvalue_mod_r"]
    assert pow(lam, n, r) == 1
    assert all(pow(lam, k, r) != 1 for k in range(1, n))
    assert E(from_word(curve["generator"][0])**2,
             from_word(curve["generator"][1])**2) == lam * G

    base = load(BASE_FILE)
    points = [point(coords) for coords in base["factor_base_point_coordinates"]]
    representatives = [point(coords) for coords in base["representative_points"]]
    labels = [tuple(map(int, label)) for label in base["factor_base_point_labels"]]
    columns = base["orbit_columns"]
    assert len(points) == len(labels) == 25864
    assert columns == len(representatives) == 244
    assert base["base_hash"] == record["factor_base"]["enumerated_set_digest"]["value"]
    assert run["factor_base_digest"] == base["base_hash"]
    orbit_coefficients = {pow(lam, j, r) for j in range(n)} | {
        (-pow(lam, j, r)) % r for j in range(n)
    }
    assert len(orbit_coefficients) == 106
    for rep in representatives:
        assert rep != E(0) and r * rep == E(0)
    for value, (column, coefficient) in zip(points, labels):
        assert 0 <= column < columns and coefficient in orbit_coefficients
        assert value == coefficient * representatives[column]
        assert value != E(0) and r * value == E(0)

    target_indices = list(map(int, run["target_relation_indices"]))
    assert len(target_indices) == 4
    assert sum((points[i] for i in target_indices), E(0)) == Q
    target_coefficients = coefficient_row(target_indices, labels, columns, r)
    witnesses = run["rank_relation_witnesses"]
    assert len(witnesses) == run["rank_attempts_completed"] == 238
    assert run["rank_attempts"] == 238 and run["rank_failures"] == 0
    rows_mod_r = []
    rhs_values = []
    pivots: dict[int, list[int]] = {}
    gains = 0
    first_span = None
    first_scalar = None
    for prefix, witness in enumerate(witnesses, start=1):
        indices = list(map(int, witness["point_indices"]))
        scalar = int(witness["relation_scalar"])
        assert len(indices) == 4
        assert sum((points[i] for i in indices), E(0)) == scalar * G
        coefficients = coefficient_row(indices, labels, columns, r)
        rows_mod_r.append(coefficients)
        rhs_values.append(scalar)
        work = coefficients + [scalar]
        gained = False
        for column in range(columns):
            factor = work[column]
            if factor == 0:
                continue
            pivot = pivots.get(column)
            if pivot is None:
                inverse = pow(factor, -1, r)
                pivots[column] = [(x * inverse) % r for x in work]
                gained = True
                break
            for j in range(column, columns + 1):
                work[j] = (work[j] - factor * pivot[j]) % r
        assert gained == bool(witness["rank_gain"])
        gains += gained
        if first_span is None and gained:
            value = streaming_span(pivots, target_coefficients, columns, r)
            if value is not None:
                first_span = prefix
                first_scalar = value
    assert gains == run["rank"] == run["rank_new_rows"] == 237
    assert first_span == run["target_span_stop_relation_prefix"] == 237
    assert first_scalar == run["recovered_scalar"] == (
        run["target_span_recovered_scalar"]
    )
    M = matrix(K, rows_mod_r)
    b = vector(K, rhs_values)
    sage_rank = M.rank()
    assert sage_rank == gains
    combination = M.transpose().solve_right(vector(K, target_coefficients))
    recovered = int(sum(alpha * rhs for alpha, rhs in zip(combination, b))) % r
    assert recovered == first_scalar and recovered * G == Q
    assert all(raw["recovered_scalar"] == recovered for raw in raw_runs)

    receipt = {
        "kind": "independent_sage_s3_pair_query_replay",
        "candidate_ids": {name: manifest["candidate_id"]
                          for name, manifest in manifests.items()},
        "workload_id": workload["workload_id"],
        "curve_id": curve["curve_id"],
        "sage_version": SAGE_VERSION,
        "sage_runtime_info_sha256": sha(HERE / "sage_runtime_info.json"),
        "replay_script_sha256": sha(Path(__file__)),
        "measurement_rows_raw_sha256": sha(ROWS_FILE),
        "base_fixture_sha256": sha(BASE_FILE),
        "run_count": len(raw_runs),
        "unique_semantic_trace_count": len(semantic_hashes),
        "semantic_trace_sha256": next(iter(semantic_hashes)),
        "all_run_artifact_hashes_checked": True,
        "all_runs_semantically_identical_excluding_clocks_memory_and_s3_counters": True,
        "factor_base_points_checked": len(points),
        "factor_base_columns": columns,
        "relation_witnesses_group_checked": len(witnesses),
        "matrix_rows": M.nrows(),
        "matrix_columns": M.ncols(),
        "matrix_rank": int(sage_rank),
        "first_target_span_relation_prefix": first_span,
        "target_relation_point_sum_checked": True,
        "recovered_scalar_from_independent_matrix_solve": recovered,
        "recovered_scalar_group_replay_checked": True,
        "replay_wall_ns": time.perf_counter_ns() - started,
        "run_ids": [row["run_id"] for row in rows],
    }
    RECEIPT_FILE.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    verified = []
    for row in rows:
        value = dict(row)
        value["status"] = "success"
        value["sage_verified"] = True
        value["sage_replay_receipt"] = str(RECEIPT_FILE)
        value["sage_replay_receipt_sha256"] = sha(RECEIPT_FILE)
        verified.append(value)
    VERIFIED_FILE.write_text("".join(json.dumps(row, sort_keys=True) + "\n"
                                     for row in verified))
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
