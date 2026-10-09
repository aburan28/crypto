"""Independently replay 396/552-column S3 relations and logs in Sage.

Launch through /Volumes/SSD990/cryptanalysis/sage -python after saving its
--runtime-info receipt. The code rebuilds each factor base from encoded
points, checks its Frobenius labels, then solves each target relation system.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, matrix, vector
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
PRIOR_BASE = ROOT / "experiments/koblitz-s3-base-size-20261009/local/n53_base_c396.json.gz"
N = 53
R = 21044858204113
G_WORDS = (198217578752339, 7929897206038174)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def base_record(columns: int) -> tuple[dict, str]:
    if columns == 396:
        raw = gzip.decompress(PRIOR_BASE.read_bytes())
    else:
        assert columns == 552
        raw = gzip.decompress((HERE / "local" / f"n53_base_c{columns}.json.gz").read_bytes())
    return json.loads(raw), sha(raw)


def load_curve():
    polynomial = PolynomialRing(GF(2), "t")
    t = polynomial.gen()
    field = GF(2**N, name="z", modulus=t**53 + t**6 + t**2 + t + 1)
    powers = [field.gen() ** bit for bit in range(N)]

    def from_word(word: int):
        assert isinstance(word, int) and 0 <= word < 1 << N
        value = field(0)
        while word:
            bit = (word & -word).bit_length() - 1
            value += powers[bit]
            word &= word - 1
        return value

    curve = EllipticCurve(field, [field(1), field(0), field(0), field(0), field(1)])

    def point(words):
        assert len(words) == 2
        return curve(from_word(int(words[0])), from_word(int(words[1])))

    generator = point(G_WORDS)
    assert generator != curve(0) and R * generator == curve(0)
    return curve, generator, point


def check_base(columns: int, curve, generator, point, expected: dict) -> dict:
    base, base_sha = base_record(columns)
    assert base["n"] == N and base["a"] == 0 and base["subgroup_order"] == R
    assert base["orbit_columns"] == columns
    assert base["base_hash"] == expected["factor_base"]["enumerated_set_digest"]
    assert base["factor_base_scanned_x"] == expected["factor_base"]["last_scanned_x"]
    encoded = base["factor_base_point_coordinates"]
    labels = [tuple(map(int, item)) for item in base["factor_base_point_labels"]]
    reps = base["representative_points"]
    assert len(reps) == columns and len(encoded) == len(labels) == 2 * N * columns
    lam = labels[2][1]
    assert 0 < lam < R and pow(lam, N, R) == 1
    frob = lambda P: curve(P[0] ** 2, P[1] ** 2)
    assert lam * generator == frob(generator)
    points = []
    for column, pair in enumerate(reps):
        representative = point(pair)
        assert R * representative == curve(0)
        current = representative
        coefficient = 1
        for shift in range(N):
            index = 2 * (column * N + shift)
            direct = point(encoded[index])
            negative = point(encoded[index + 1])
            assert direct == current and negative == -current
            assert labels[index] == (column, coefficient)
            assert labels[index + 1] == (column, (-coefficient) % R)
            points.extend((direct, negative))
            current = frob(current)
            coefficient = coefficient * lam % R
        assert current == representative and coefficient == 1
    return {"base_sha256": base_sha, "base_digest": base["base_hash"],
            "points": points, "labels": labels,
            "base_point_count_checked": len(points),
            "representatives_subgroup_checked": columns,
            "frobenius_eigenvalue": lam}


def row_coefficients(indices, labels, columns: int):
    row = [0] * columns
    for index in indices:
        column, coefficient = labels[index]
        row[column] = (row[column] + coefficient) % R
    return row


def replay(ledger: Path, runtime: Path, receipt: Path) -> dict:
    if receipt.exists():
        raise FileExistsError(receipt)
    checked_runtime = json.loads(runtime.read_text())
    assert checked_runtime["status"] == "verified"
    assert checked_runtime["sage_version"] == SAGE_VERSION
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    assert plan["parent_plan_sha256"] == sha(parent_bytes)
    parent = json.loads(parent_bytes)
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    raw = ledger.read_bytes()
    ledger_bytes = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
    rows = [json.loads(line) for line in ledger_bytes.splitlines()]
    assert len(rows) == 32
    curve, generator, point = load_curve()
    bases = {}
    for arm in ("baseline", "expanded"):
        frozen = plan["arms"][arm]
        columns = frozen["columns"]
        candidate = frozen["candidate_id"]
        path = HERE / "identities/candidates" / f"{candidate}.json"
        expected = json.loads(path.read_text())
        assert expected.pop("candidate_id") == candidate
        bases[columns] = check_base(columns, curve, generator, point, expected)
        assert bases[columns]["base_digest"] == frozen["factor_base_digest"]
        assert len(bases[columns]["points"]) == frozen["factor_base_points"]
    assert bases[552]["base_sha256"] == plan["base552_raw_sha256"]

    field = GF(R)
    traces = {}
    seen_witnesses = set()
    checked = []
    for entry in rows:
        frozen = plan["arms"][entry["arm"]]
        columns = frozen["columns"]
        assert entry["candidate_id"] == frozen["candidate_id"]
        assert entry["source_sha256"] == frozen["source_sha256"]
        assert entry["binary_sha256"] == frozen["binary_sha256"]
        base = bases[columns]
        target = targets[(N, entry["label"])]
        assert entry["run_id"] == f"{frozen['candidate_id']}W{target['workload_id']}R{entry['repetition']}"
        report = entry["report"]
        assert entry["verified_native_and_fixture"] and entry["phase_ok"]
        assert report["group_verified"] and report["target_count"] == 1
        assert report["factor_base_digest"] == base["base_digest"]
        assert report["factor_base_points"] == len(base["points"])
        assert report["target"] == target["public_point"]
        assert report["recovered_scalar"] == target["fixture_scalar"]
        public_point = point(target["public_point"])
        assert int(target["fixture_scalar"]) * generator == public_point
        indices = list(map(int, report["target_relation_indices"]))
        assert len(indices) == 4
        assert sum((base["points"][index] for index in indices), curve(0)) == public_point
        target_row = vector(field, row_coefficients(indices, base["labels"], columns))
        witnesses = report["rank_relation_witnesses"] + report["fallback_relation_witnesses"]
        trace_hash = sha(json.dumps(witnesses, sort_keys=True, separators=(",", ":")).encode())
        key = (columns, trace_hash)
        if key not in traces:
            coefficients = []
            rhs = []
            for witness in witnesses:
                witness_indices = tuple(map(int, witness["point_indices"]))
                scalar = int(witness["relation_scalar"])
                assert len(witness_indices) == 4
                witness_key = (columns, scalar, witness_indices)
                if witness_key not in seen_witnesses:
                    assert sum((base["points"][index] for index in witness_indices),
                               curve(0)) == scalar * generator
                    seen_witnesses.add(witness_key)
                coefficients.append(row_coefficients(witness_indices, base["labels"], columns))
                rhs.append(scalar)
            relation_matrix = matrix(field, coefficients)
            relation_rhs = vector(field, rhs)
            rank = relation_matrix.rank()
            assert rank == report["final_rank"]
            traces[key] = (relation_matrix, relation_rhs, rank)
        relation_matrix, relation_rhs, rank = traces[key]
        if rank == columns:
            logs = relation_matrix.solve_right(relation_rhs)
            recovered = int(target_row.dot_product(logs)) % R
        else:
            weights = relation_matrix.transpose().solve_right(target_row)
            recovered = int(weights.dot_product(relation_rhs)) % R
        assert recovered == target["fixture_scalar"]
        assert recovered * generator == public_point
        checked.append({"columns": columns, "label": entry["label"],
                        "run_id": entry["run_id"], "rank": int(rank),
                        "independent_recovery": True})
    assert len(checked) == len(rows)
    result = {"kind": "s3_base552_independent_sage_replay_v1",
              "sage_version": SAGE_VERSION,
              "runtime_info_sha256": sha(runtime.read_bytes()),
              "freeze_sha256": sha(plan_bytes),
              "parent_plan_sha256": sha(parent_bytes),
              "ledger_sha256": sha(ledger_bytes),
              "run_count": len(rows),
              "target_count": len({(row["n"], row["label"]) for row in rows}),
              "factor_bases": {str(columns): {key: value for key, value in base.items()
                                             if key not in ("points", "labels")}
                               for columns, base in bases.items()},
              "distinct_relation_witnesses_group_checked": len(seen_witnesses),
              "distinct_relation_traces_rank_checked": len(traces),
              "all_scalars_independently_recovered": True,
              "runs": checked}
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(receipt), "runs": len(rows),
                      "witnesses": len(seen_witnesses), "traces": len(traces)}), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=HERE / "local/raw_paired32.jsonl.gz")
    parser.add_argument("--runtime-info", type=Path,
                        default=HERE / "local/sage_runtime_info.json")
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    replay(args.ledger, args.runtime_info, args.receipt)
