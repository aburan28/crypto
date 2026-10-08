"""Independently replay all distinct N41/N53 Modal pilot targets in Sage.

Save /Volumes/SSD990/cryptanalysis/sage --runtime-info to
modal/sage_runtime_info.json before launching this file through the checked
repository Sage launcher. This validates arithmetic and relations, not host
timing or physical CPU isolation.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import time

from sage.all import EllipticCurve, GF, PolynomialRing, matrix, vector
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
MODAL = HERE / "modal"
PANEL = MODAL / "pilot_raw_v2.json.gz"
RUNTIME = MODAL / "sage_runtime_info.json"
RECEIPT = MODAL / "independent_sage_replay.json"
OLD_N53_BASE = HERE.parent / "koblitz-s3-pair-query-20261007-v4/fb244_preflight.json"
CONFIG = {
    41: {"r": 549756390943, "generator": [2056947637384, 1635505394702],
         "modulus_powers": [41, 3, 0], "base": MODAL / "n41_base_dump.json",
         "expected_points": 20008},
    53: {"r": 21044858204113, "generator": [198217578752339, 7929897206038174],
         "modulus_powers": [53, 6, 2, 1, 0], "base": OLD_N53_BASE,
         "expected_points": 25864},
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def coefficient_row(indices: list[int], labels: list[tuple[int, int]],
                    columns: int, modulus: int) -> list[int]:
    row = [0] * columns
    for index in indices:
        column, coefficient = labels[index]
        row[column] = (row[column] + coefficient) % modulus
    return row


def streaming_span(pivots: dict[int, list[int]], target: list[int],
                   columns: int, modulus: int) -> int | None:
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
            remainder[j] = (remainder[j] - factor * pivot[j]) % modulus
        value = (value + factor * pivot[columns]) % modulus
    return value if not any(remainder) else None


def semantic_witness(run: dict) -> tuple:
    return (
        run["factor_base_digest"], run["factor_base_points"], run["orbit_columns"],
        run["rank"], run["rank_attempts"], run["rank_new_rows"],
        run["target_relation_indices"], run["target_span_stop_relation_prefix"],
        run["recovered_scalar"],
        [(row["point_indices"], row["rank_gain"], row["relation_scalar"], row["x_codes"])
         for row in run["rank_relation_witnesses"]],
    )


def replay_curve(n: int, panel: dict) -> dict:
    cfg = CONFIG[n]
    modulus = cfg["r"]
    P = PolynomialRing(GF(2), "t")
    t = P.gen()
    field = GF(2**n, name="z", modulus=sum(t**power for power in cfg["modulus_powers"]))
    powers = [field.gen() ** bit for bit in range(n)]

    def from_word(word: int):
        assert isinstance(word, int) and 0 <= word < 1 << n
        value = field(0)
        while word:
            bit = (word & -word).bit_length() - 1
            value += powers[bit]
            word &= word - 1
        return value

    curve = EllipticCurve(field, [field(1), field(0), field(0), field(0), field(1)])

    def point(words: list[int]):
        assert len(words) == 2
        return curve(from_word(words[0]), from_word(words[1]))

    generator = point(cfg["generator"])
    assert generator != curve(0) and modulus * generator == curve(0)
    base_file = cfg["base"]
    base = load(base_file)
    assert base["n"] == n and base["a"] == 0 and base["subgroup_order"] == modulus
    assert base["orbit_columns"] == 244
    coordinates = base["factor_base_point_coordinates"]
    labels = [tuple(map(int, row)) for row in base["factor_base_point_labels"]]
    assert len(coordinates) == len(labels) == cfg["expected_points"]
    representatives = [point(row) for row in base["representative_points"]]
    assert len(representatives) == 244
    points = [point(row) for row in coordinates]
    assert len(set(map(tuple, coordinates))) == len(coordinates)
    for rep in representatives:
        assert rep != curve(0) and modulus * rep == curve(0)
    for value, (column, coefficient) in zip(points, labels):
        assert 0 <= column < 244 and value != curve(0)
        assert value == coefficient * representatives[column]
        assert modulus * value == curve(0)

    fixtures_file = HERE / f"pilot/n{n}/fixtures.json"
    fixtures = load(fixtures_file)["fixtures"]
    blocks = {block["label"]: block for block in panel["blocks"] if block["n"] == n}
    runs = {run["tag"]: run for run in panel["runs"]}
    assert len(fixtures) == len(blocks) == 12
    results = []
    for fixture in fixtures:
        block = blocks[fixture["label"]]
        target_file = HERE / f"pilot/n{n}/{fixture['label']}/public_target.json"
        target_file_sha256 = sha(target_file.read_bytes())
        assert load(target_file) == fixture["public_point"]
        target = point(fixture["public_point"])
        assert target != curve(0) and modulus * target == curve(0)
        assert int(fixture["fixture_scalar"]) * generator == target
        block_runs = [runs[tag] for tag in block["run_tags"]]
        assert len(block_runs) == 6
        assert all(run["verified_native_and_fixture"] and run["exit_code"] == 0
                   for run in block_runs)
        assert all(run["target_file_sha256"] == target_file_sha256 for run in block_runs)
        assert all(run["report"]["target"] == fixture["public_point"] and
                   run["report"]["recovered_scalar"] == int(fixture["fixture_scalar"]) and
                   run["report"]["factor_base_digest"] == base["base_hash"] and
                   run["report"]["factor_base_points"] == len(points)
                   for run in block_runs)
        assert block["semantic_equal"]
        assert all(semantic_witness(run["report"]) == semantic_witness(block_runs[0]["report"])
                   for run in block_runs)
        run = block_runs[0]["report"]
        target_indices = list(map(int, run["target_relation_indices"]))
        assert len(target_indices) == 4
        assert sum((points[i] for i in target_indices), curve(0)) == target
        target_row = coefficient_row(target_indices, labels, 244, modulus)
        witnesses = run["rank_relation_witnesses"]
        assert len(witnesses) == run["rank_attempts_completed"] == run["rank_attempts"]
        assert run["rank_failures"] == 0
        rows = []
        rhs = []
        pivots: dict[int, list[int]] = {}
        gained_rows = 0
        first_span = None
        first_scalar = None
        for prefix, witness in enumerate(witnesses, 1):
            indices = list(map(int, witness["point_indices"]))
            scalar = int(witness["relation_scalar"])
            assert len(indices) == 4
            assert sum((points[i] for i in indices), curve(0)) == scalar * generator
            coeffs = coefficient_row(indices, labels, 244, modulus)
            rows.append(coeffs)
            rhs.append(scalar)
            work = coeffs + [scalar]
            gained = False
            for column in range(244):
                factor = work[column]
                if factor == 0:
                    continue
                pivot = pivots.get(column)
                if pivot is None:
                    inverse = pow(factor, -1, modulus)
                    pivots[column] = [(value * inverse) % modulus for value in work]
                    gained = True
                    break
                for j in range(column, 245):
                    work[j] = (work[j] - factor * pivot[j]) % modulus
            assert gained == bool(witness["rank_gain"])
            gained_rows += gained
            if first_span is None and gained:
                maybe_scalar = streaming_span(pivots, target_row, 244, modulus)
                if maybe_scalar is not None:
                    first_span = prefix
                    first_scalar = maybe_scalar
        assert gained_rows == run["rank"] == run["rank_new_rows"]
        assert first_span == run["target_span_stop_relation_prefix"]
        assert first_scalar == run["recovered_scalar"] == int(fixture["fixture_scalar"])
        K = GF(modulus)
        mat = matrix(K, rows)
        assert mat.rank() == gained_rows
        combination = mat.transpose().solve_right(vector(K, target_row))
        independently_recovered = int(sum(c * v for c, v in zip(combination, vector(K, rhs)))) % modulus
        assert independently_recovered == first_scalar
        assert independently_recovered * generator == target
        results.append({"label": fixture["label"],
                        "public_target_sha256": target_file_sha256,
                        "representative_raw_sha256": block_runs[0]["raw_sha256"],
                        "relation_witnesses_group_checked": len(witnesses),
                        "matrix_rank": int(mat.rank()),
                        "first_target_span_relation_prefix": first_span,
                        "recovered_scalar": independently_recovered,
                        "all_six_runs_semantically_equal": True})
        print(json.dumps({"n": n, "target": fixture["label"],
                          "relations": len(witnesses), "recovered": True}), flush=True)
    return {"n": n, "field_modulus_powers": cfg["modulus_powers"],
            "subgroup_order": modulus, "generator": cfg["generator"],
            "factor_base_fixture": str(base_file),
            "factor_base_fixture_sha256": sha(base_file.read_bytes()),
            "factor_base_digest": base["base_hash"],
            "factor_base_points_checked": len(points), "folded_columns_checked": 244,
            "fixtures_sha256": sha(fixtures_file.read_bytes()),
            "relation_witnesses_group_checked": sum(row["relation_witnesses_group_checked"]
                                                    for row in results),
            "targets": results}


def main() -> None:
    if RECEIPT.exists():
        raise FileExistsError(RECEIPT)
    runtime = load(RUNTIME)
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    started = time.perf_counter_ns()
    panel_bytes = PANEL.read_bytes()
    panel = json.loads(gzip.decompress(panel_bytes))
    assert len(panel["blocks"]) == 24 and len(panel["runs"]) == 146
    assert all(run["raw_text"] is not None and sha(run["raw_text"].encode()) == run["raw_sha256"]
               and json.loads(run["raw_text"]) == run["report"]
               for run in panel["runs"])
    curves = [replay_curve(n, panel) for n in (41, 53)]
    receipt = {"kind": "modal_two_curve_independent_sage_replay_v1",
               "scope": "all distinct public points, full factor bases, one semantic relation trace per target, rank and target scalar; no host timing certification",
               "sage_version": SAGE_VERSION,
               "sage_runtime_info_sha256": sha(RUNTIME.read_bytes()),
               "replay_script_sha256": sha(Path(__file__).read_bytes()),
               "panel_sha256": sha(panel_bytes),
               "target_count": sum(len(c["targets"]) for c in curves),
               "run_count_with_exact_raw_hash_checked": len(panel["runs"]),
               "relation_witnesses_group_checked": sum(c["relation_witnesses_group_checked"]
                                                       for c in curves),
               "all_target_scalars_independently_recovered": True,
               "all_factor_base_points_and_labels_checked": True,
               "curves": curves, "replay_wall_ns": time.perf_counter_ns() - started}
    RECEIPT.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(RECEIPT), "targets": receipt["target_count"],
                      "relations": receipt["relation_witnesses_group_checked"]}), flush=True)


if __name__ == "__main__":
    main()
