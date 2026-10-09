"""Independently replay the paired prepared-S3 panel with checked Sage.

Run only through /Volumes/SSD990/cryptanalysis/sage -python after saving
the launcher's --runtime-info receipt. The prior checked Sage receipt covers
every factor-base point and label with these exact enumerated-set digests.
This replay checks all new target witnesses and distinct relation traces.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

from sage.all import EllipticCurve, GF, PolynomialRing, matrix, vector
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONFIG = {
    41: {"r": 549756390943, "generator": [2056947637384, 1635505394702],
         "modulus_powers": [41, 3, 0], "base": HERE / "n41_base_dump.json",
         "expected_points": 20008},
    53: {"r": 21044858204113, "generator": [198217578752339, 7929897206038174],
         "modulus_powers": [53, 6, 2, 1, 0],
         "base": ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/fb244_preflight.json",
         "expected_points": 25864},
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def row_coefficients(indices: list[int], labels: list[tuple[int, int]], r: int) -> list[int]:
    row = [0] * 244
    for index in indices:
        column, coefficient = labels[index]
        row[column] = (row[column] + coefficient) % r
    return row


def setup_curve(n: int, prior: dict):
    cfg = CONFIG[n]
    base_path = cfg["base"]
    base = json.loads(base_path.read_text())
    old = next(item for item in prior["curves"] if item["n"] == n)
    assert prior["all_factor_base_points_and_labels_checked"] is True
    assert old["factor_base_points_checked"] == cfg["expected_points"]
    assert old["factor_base_digest"] == base["base_hash"]
    assert base["n"] == n and base["a"] == 0 and base["subgroup_order"] == cfg["r"]
    assert base["orbit_columns"] == 244
    assert len(base["factor_base_point_coordinates"]) == cfg["expected_points"]
    assert len(base["factor_base_point_labels"]) == cfg["expected_points"]

    polynomial = PolynomialRing(GF(2), "t")
    t = polynomial.gen()
    field = GF(2**n, name="z", modulus=sum(t**p for p in cfg["modulus_powers"]))
    powers = [field.gen() ** bit for bit in range(n)]

    def from_word(word: int):
        assert isinstance(word, int) and 0 <= word < 1 << n
        result = field(0)
        while word:
            bit = (word & -word).bit_length() - 1
            result += powers[bit]
            word &= word - 1
        return result

    curve = EllipticCurve(field, [field(1), field(0), field(0), field(0), field(1)])

    def point(words: list[int]):
        assert len(words) == 2
        return curve(from_word(words[0]), from_word(words[1]))

    generator = point(cfg["generator"])
    assert generator != curve(0) and cfg["r"] * generator == curve(0)
    representatives = [point(pair) for pair in base["representative_points"]]
    points = [point(pair) for pair in base["factor_base_point_coordinates"]]
    labels = [tuple(map(int, label)) for label in base["factor_base_point_labels"]]
    assert len(representatives) == 244 and len(points) == len(labels)
    return cfg, base, curve, generator, representatives, points, labels, point, sha(base_path.read_bytes())


def replay_curve(n: int, rows: list[dict], prior: dict) -> dict:
    (cfg, base, curve, generator, representatives, points, labels, point,
     base_file_sha) = setup_curve(n, prior)
    r = cfg["r"]
    K = GF(r)
    fixtures_file = ROOT / f"experiments/koblitz-s3-fidelity-20261007/pilot/n{n}/fixtures.json"
    fixtures = {item["label"]: item for item in json.loads(fixtures_file.read_text())["fixtures"]}
    assert len(fixtures) == 12 and len(rows) == 24
    traces: dict[str, tuple] = {}
    prepared_logs: dict[str, object] = {}
    checked_point_labels = set()
    checked_targets = []
    checked_relations = 0
    for entry in rows:
        fixture = fixtures[entry["label"]]
        report = entry["report"]
        assert entry["verified_native_and_fixture"] is True and entry["phase_ok"] is True
        assert report["group_verified"] is True and report["target_count"] == 1
        assert report["factor_base_digest"] == base["base_hash"]
        assert report["factor_base_points"] == len(points) and report["orbit_columns"] == 244
        assert report["target"] == fixture["public_point"]
        assert report["recovered_scalar"] == int(fixture["fixture_scalar"])
        target = point(fixture["public_point"])
        assert int(fixture["fixture_scalar"]) * generator == target
        target_indices = list(map(int, report["target_relation_indices"]))
        assert len(target_indices) == 4
        assert sum((points[index] for index in target_indices), curve(0)) == target
        checked_point_labels.update(target_indices)
        target_row = row_coefficients(target_indices, labels, r)
        witnesses = report["rank_relation_witnesses"]
        trace_hash = sha(json.dumps(witnesses, sort_keys=True, separators=(",", ":")).encode())
        if trace_hash not in traces:
            assert len(witnesses) == report["rank_verified_relations"]
            matrix_rows = []
            right_hand_side = []
            for witness in witnesses:
                indices = list(map(int, witness["point_indices"]))
                scalar = int(witness["relation_scalar"])
                assert len(indices) == 4
                assert sum((points[index] for index in indices), curve(0)) == scalar * generator
                checked_point_labels.update(indices)
                matrix_rows.append(row_coefficients(indices, labels, r))
                right_hand_side.append(scalar)
                checked_relations += 1
            M = matrix(K, matrix_rows)
            rank = int(M.rank())
            assert rank == report["rank"] == report["rank_new_rows"]
            traces[trace_hash] = (M, vector(K, right_hand_side), rank)
        M, rhs, rank = traces[trace_hash]
        if entry["arm"] == "prepared":
            assert report["rank_complete"] is True and rank == 244
            if trace_hash not in prepared_logs:
                prepared_logs[trace_hash] = M.solve_right(rhs)
            independent_logs = prepared_logs[trace_hash]
            assert [int(value) for value in independent_logs] == report["solved_base_logs"]
            recovered = int(vector(K, target_row).dot_product(independent_logs)) % r
        else:
            weights = M.transpose().solve_right(vector(K, target_row))
            recovered = int(weights.dot_product(rhs)) % r
        assert recovered == int(fixture["fixture_scalar"])
        assert recovered * generator == target
        checked_targets.append({"label": entry["label"], "arm": entry["arm"],
                                "raw_sha256": entry["raw_sha256"], "trace_sha256": trace_hash,
                                "rank": rank, "recovered_scalar": recovered})
        print(json.dumps({"n": n, "target": entry["label"], "arm": entry["arm"],
                          "rank": rank, "independent_recovery": True}), flush=True)
    for index in checked_point_labels:
        column, coefficient = labels[index]
        assert points[index] == coefficient * representatives[column]
    assert len({(item["label"], item["arm"]) for item in checked_targets}) == 24
    return {"n": n, "field_modulus_powers": cfg["modulus_powers"],
            "subgroup_order": r, "generator": cfg["generator"],
            "base_fixture_sha256": base_file_sha, "base_digest": base["base_hash"],
            "factor_base_points_prior_checked": len(points),
            "point_labels_rechecked": len(checked_point_labels),
            "distinct_relation_traces": len(traces),
            "relation_witnesses_group_checked": checked_relations,
            "targets": checked_targets}


def main(ledger_path: Path, runtime_path: Path, output_path: Path) -> None:
    if output_path.exists():
        raise FileExistsError(output_path)
    runtime = json.loads(runtime_path.read_text())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    prior_path = HERE / "prior_base_replay.json"
    prior = json.loads(prior_path.read_text())
    assert prior["target_count"] == 24 and prior["all_target_scalars_independently_recovered"] is True
    started = time.perf_counter_ns()
    ledger_bytes = ledger_path.read_bytes()
    rows = [json.loads(line) for line in ledger_bytes.splitlines()]
    assert len(rows) == 48
    for entry in rows:
        report = entry["report"]
        assert report is not None and entry["error"] is None and entry["exit_code"] == 0
        raw_path = ledger_path.parent / f"n{entry['n']}_{entry['label']}_{entry['arm']}.jsonl"
        assert sha(raw_path.read_bytes()) == entry["raw_sha256"]
        assert json.loads(raw_path.read_text()) == report
    curves = [replay_curve(n, [row for row in rows if row["n"] == n], prior) for n in (41, 53)]
    result = {"kind": "s3_prepared_local_independent_sage_replay_v1",
              "sage_version": SAGE_VERSION,
              "sage_runtime_info_sha256": sha(runtime_path.read_bytes()),
              "replay_script_sha256": sha(Path(__file__).read_bytes()),
              "prior_base_replay_sha256": sha(prior_path.read_bytes()),
              "ledger_sha256": sha(ledger_bytes), "run_count": len(rows),
              "target_count": 24, "all_scalars_independently_recovered": True,
              "relation_witnesses_group_checked": sum(c["relation_witnesses_group_checked"] for c in curves),
              "curves": curves, "replay_wall_ns": time.perf_counter_ns() - started}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(output_path), "runs": len(rows),
                      "witnesses": result["relation_witnesses_group_checked"]}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--runtime-info", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    main(args.ledger, args.runtime_info, args.receipt)
