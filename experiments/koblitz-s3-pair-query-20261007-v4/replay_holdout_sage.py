#!/usr/bin/env sage -python
"""Independently verify the twelve fresh one-target S3 workloads."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import statistics
import time

from sage.all import EllipticCurve, GF, PolynomialRing, matrix, vector
from sage.env import SAGE_VERSION

from replay_sage import coefficient_row, streaming_span
from run_panel import canonical, semantic_digest


HERE = Path(__file__).resolve().parent
HOLDOUT = HERE / "holdout"
RAW_ROWS = HOLDOUT / "measurement_rows_raw.jsonl"
RECEIPT = HOLDOUT / "independent_sage_replay.json"
VERIFIED_ROWS = HOLDOUT / "measurement_rows_verified.jsonl"
FINAL_SUMMARY = HOLDOUT / "final_summary_exploratory.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def main() -> None:
    for path in (RECEIPT, VERIFIED_ROWS, FINAL_SUMMARY):
        if path.exists():
            raise FileExistsError(f"refusing to replace replay artifact: {path}")
    started = time.perf_counter_ns()
    fixtures_path = HOLDOUT / "fixtures.json"
    fixtures = load(fixtures_path)
    summary = load(HOLDOUT / "paired_summary_exploratory.json")
    raw_rows = [json.loads(line) for line in RAW_ROWS.read_text().splitlines() if line]
    assert len(fixtures["fixtures"]) == 12 and len(raw_rows) == 24
    assert len({row["run_id"] for row in raw_rows}) == 24
    assert summary["fixtures_sha256"] == sha(fixtures_path)
    assert summary["measurement_rows_raw_sha256"] == sha(RAW_ROWS)
    assert summary["all_pairs_semantically_equal"] is True
    primary_replay = load(HERE / "independent_sage_replay.json")
    manifests = {v: load(HERE / f"{v}_manifest.json")
                 for v in ("baseline", "candidate")}
    frozen = load(HERE / "freeze_receipt.json")
    for variant, manifest in manifests.items():
        candidate = manifest["candidate_freeze"]
        assert sha(HERE / f"{variant}_manifest.json") == (
            frozen["variants"][variant]["manifest_sha256"]
        )
        assert sha(Path(candidate["source_path"])) == candidate["source_sha256"]
        assert sha(Path(candidate["binary_path"])) == candidate["binary_sha256"]
        assert manifest["candidate_id"] == (
            fixtures["source_candidate_ids"][variant]
        )
    for row in raw_rows:
        manifest = manifests[row["variant"]]
        assert row["candidate_id"] == manifest["candidate_id"]
        workload = load(HOLDOUT / row["target_label"] / "workload.json")
        assert row["workload_id"] == workload["workload_id"]
        assert row["run_id"] == f"{row['candidate_id']}W{row['workload_id']}R1"
        assert workload["identity_sha256"] == hashlib.sha256(
            canonical(workload["record"])
        ).hexdigest()
        assert workload["record"]["targets"] == [row["target"]]
        raw_path = Path(row["artifacts"]["raw"])
        assert sha(raw_path) == row["artifacts"]["raw_sha256"]
        for name in ("stdout", "stderr"):
            assert sha(Path(row["artifacts"][name])) == (
                row["artifacts"][f"{name}_sha256"]
            )
        run = load(raw_path)
        assert semantic_digest(run) == (
            row["semantic_digest_excluding_clocks_memory_and_s3_counters"]
        )
        assert row["status"] == "reported_success_pending_sage"
        assert run["group_verified"] is True and run["target"] == row["target"]
        assert abs(sum(row["online_exclusive_phases_ms"].values())
                   - row["online_ms"]) < 0.02

    record = manifests["candidate"]["identity_record"]
    curve = record["curve"]
    n = record["field"]["degree"]
    r = curve["subgroup_order"]
    assert n == 53
    R = PolynomialRing(GF(2), "t")
    t = R.gen()
    F = GF(2**n, name="z", modulus=t**53 + t**6 + t**2 + t + 1)
    powers = [F.gen() ** i for i in range(n)]

    def from_word(word: int):
        word = int(word)
        assert 0 <= word < 1 << n
        value = F(0)
        while word:
            bit = (word & -word).bit_length() - 1
            value += powers[bit]
            word &= word - 1
        return value

    E = EllipticCurve(F, [F(1), F(0), F(0), F(0), F(1)])
    K = GF(r)

    def point(words: list[int]):
        return E(from_word(words[0]), from_word(words[1]))

    G = point(curve["generator"])
    assert G != E(0) and r * G == E(0)
    lam = record["endomorphism"]["frobenius_eigenvalue_mod_r"]
    assert pow(lam, n, r) == 1
    assert E(from_word(curve["generator"][0])**2,
             from_word(curve["generator"][1])**2) == lam * G
    base_path = HERE / "fb244_preflight.json"
    base = load(base_path)
    assert sha(base_path) == primary_replay["base_fixture_sha256"]
    points = [point(coords) for coords in base["factor_base_point_coordinates"]]
    representatives = [point(coords) for coords in base["representative_points"]]
    labels = [tuple(map(int, label)) for label in base["factor_base_point_labels"]]
    columns = int(base["orbit_columns"])
    assert len(points) == len(labels) == 25864
    assert len(representatives) == columns == 244
    assert base["base_hash"] == record["factor_base"]["enumerated_set_digest"]["value"]
    orbit_coefficients = {pow(lam, j, r) for j in range(n)} | {
        (-pow(lam, j, r)) % r for j in range(n)
    }
    for rep in representatives:
        assert rep != E(0) and r * rep == E(0)
    for value, (column, coefficient) in zip(points, labels):
        assert 0 <= column < columns and coefficient in orbit_coefficients
        assert value == coefficient * representatives[column]
        assert value != E(0) and r * value == E(0)

    run_by_key = {(row["target_label"], row["variant"]): row for row in raw_rows}
    results = []
    for fixture in fixtures["fixtures"]:
        label = fixture["label"]
        baseline_row = run_by_key[(label, "baseline")]
        candidate_row = run_by_key[(label, "candidate")]
        baseline_run = load(Path(baseline_row["artifacts"]["raw"]))
        candidate_run = load(Path(candidate_row["artifacts"]["raw"]))
        assert semantic_digest(baseline_run) == semantic_digest(candidate_run)
        run = candidate_run
        Q = point(fixture["public_point"])
        assert Q != E(0) and r * Q == E(0)
        assert fixture["fixture_scalar"] * G == Q
        target_indices = list(map(int, run["target_relation_indices"]))
        assert len(target_indices) == 4
        assert sum((points[i] for i in target_indices), E(0)) == Q
        target_row = coefficient_row(target_indices, labels, columns, r)
        witnesses = run["rank_relation_witnesses"]
        assert len(witnesses) == run["rank_attempts_completed"] == run["rank_attempts"]
        assert run["rank_failures"] == 0
        matrix_rows = []
        rhs_values = []
        pivots: dict[int, list[int]] = {}
        gained_rows = 0
        first_span = None
        first_scalar = None
        for prefix, witness in enumerate(witnesses, 1):
            indices = list(map(int, witness["point_indices"]))
            scalar = int(witness["relation_scalar"])
            assert len(indices) == 4
            assert sum((points[i] for i in indices), E(0)) == scalar * G
            coefficients = coefficient_row(indices, labels, columns, r)
            matrix_rows.append(coefficients)
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
            gained_rows += gained
            if first_span is None and gained:
                value = streaming_span(pivots, target_row, columns, r)
                if value is not None:
                    first_span = prefix
                    first_scalar = value
        assert gained_rows == run["rank"] == run["rank_new_rows"]
        assert first_span == run["target_span_stop_relation_prefix"]
        assert first_scalar == run["recovered_scalar"] == fixture["fixture_scalar"]
        M = matrix(K, matrix_rows)
        assert M.rank() == gained_rows
        b = vector(K, rhs_values)
        combination = M.transpose().solve_right(vector(K, target_row))
        recovered = int(sum(a * v for a, v in zip(combination, b))) % r
        assert recovered == first_scalar and recovered * G == Q
        results.append({
            "target_label": label,
            "workload_id": fixture["workload_id"],
            "baseline_run_id": baseline_row["run_id"],
            "candidate_run_id": candidate_row["run_id"],
            "baseline_raw_sha256": baseline_row["artifacts"]["raw_sha256"],
            "candidate_raw_sha256": candidate_row["artifacts"]["raw_sha256"],
            "relation_witnesses_group_checked": len(witnesses),
            "matrix_rank": int(M.rank()),
            "first_target_span_relation_prefix": first_span,
            "fixture_scalar": fixture["fixture_scalar"],
            "recovered_scalar_from_independent_matrix_solve": recovered,
            "target_point_replay_checked": True,
            "semantic_pair_equal": True,
        })
        print(json.dumps({"target_label": label, "rank": int(M.rank()),
                          "witnesses": len(witnesses), "verified": True}), flush=True)

    receipt = {
        "kind": "independent_sage_fresh_target_s3_replay",
        "sage_version": SAGE_VERSION,
        "sage_runtime_info_sha256": sha(HERE / "sage_runtime_info.json"),
        "replay_script_sha256": sha(Path(__file__)),
        "fixtures_sha256": sha(fixtures_path),
        "measurement_rows_raw_sha256": sha(RAW_ROWS),
        "base_fixture_sha256": sha(base_path),
        "factor_base_points_checked": len(points),
        "target_count": len(results),
        "raw_run_count": len(raw_rows),
        "relation_witnesses_group_checked": sum(
            item["relation_witnesses_group_checked"] for item in results
        ),
        "all_target_scalars_independently_recovered": True,
        "all_target_point_replays_checked": True,
        "all_semantic_pairs_equal": True,
        "replay_wall_ns": time.perf_counter_ns() - started,
        "targets": results,
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    verified = []
    for row in raw_rows:
        value = dict(row)
        value["status"] = "success"
        value["sage_verified"] = True
        value["sage_replay_receipt"] = str(RECEIPT)
        value["sage_replay_receipt_sha256"] = sha(RECEIPT)
        verified.append(value)
    VERIFIED_ROWS.write_text("".join(json.dumps(row, sort_keys=True) + "\n"
                                     for row in verified))
    ratios = [pair["exploratory_online_ratio"] for pair in summary["pairs"]]
    final = {
        "kind": "s3_pair_fresh_one_target_panel_verified_exploratory",
        "question": fixtures["question"],
        "target_count": len(results),
        "successful_verified_pairs": len(results),
        "verified_raw_runs": len(verified),
        "baseline_faster_pair_count": sum(ratio < 1 for ratio in ratios),
        "candidate_faster_pair_count": sum(ratio > 1 for ratio in ratios),
        "exploratory_paired_ratio_median": statistics.median(ratios),
        "exploratory_paired_ratio_range": [min(ratios), max(ratios)],
        "controlled_speedup": None,
        "claim_limit": (
            "No auditable host-wide CPU isolation receipt; timing ratios are "
            "exploratory. Twelve deterministic targets are a small correctness "
            "sample, not a demonstrated general success rate."
        ),
        "independent_sage_replay_sha256": sha(RECEIPT),
        "measurement_rows_verified_sha256": sha(VERIFIED_ROWS),
    }
    FINAL_SUMMARY.write_text(json.dumps(final, indent=2, sort_keys=True) + "\n")
    print(json.dumps(final, sort_keys=True))


if __name__ == "__main__":
    main()
