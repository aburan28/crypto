"""Independently replay both-curve dense-root relations and target scalars."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, matrix, vector
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
BASES = {
    41: ROOT / "experiments/koblitz-s3-n41-c800-fullrank-20261009/local/n41_base_c800.json.gz",
    53: ROOT / "experiments/koblitz-s3-target-conjugate-20261009/local/n53_base_c1200.json.gz",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def load_curve(n: int, parent: dict):
    curve_info = parent["curves"][str(n)]
    polynomial = PolynomialRing(GF(2), "t")
    t = polynomial.gen()
    modulus = sum((t**power for power in curve_info["modulus_powers"]), polynomial(0))
    field = GF(2**n, name="z", modulus=modulus)
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

    def point(words):
        assert len(words) == 2
        return curve(from_word(int(words[0])), from_word(int(words[1])))

    generator = point(curve_info["generator"])
    r = curve_info["r"]
    assert generator != curve(0) and r * generator == curve(0)
    return curve, generator, point, r


def load_base(n: int, columns: int, curve, generator, point, r: int, expected_digest: str):
    raw = gzip.decompress(BASES[n].read_bytes())
    base = json.loads(raw)
    assert base["n"] == n and base["a"] == 0 and base["subgroup_order"] == r
    assert base["orbit_columns"] == columns and base["base_hash"] == expected_digest
    encoded = base["factor_base_point_coordinates"]
    labels = [tuple(map(int, item)) for item in base["factor_base_point_labels"]]
    reps = base["representative_points"]
    assert len(reps) == columns and len(encoded) == len(labels) == 2 * n * columns
    lam = labels[2][1]
    frob = lambda P: curve(P[0] ** 2, P[1] ** 2)
    assert 0 < lam < r and pow(lam, n, r) == 1
    assert lam * generator == frob(generator)
    points = []
    for column, words in enumerate(reps):
        representative = point(words)
        assert r * representative == curve(0)
        current = representative
        coefficient = 1
        for shift in range(n):
            index = 2 * (column * n + shift)
            direct, negative = point(encoded[index]), point(encoded[index + 1])
            assert direct == current and negative == -current
            assert labels[index] == (column, coefficient)
            assert labels[index + 1] == (column, (-coefficient) % r)
            points.extend((direct, negative))
            current = frob(current)
            coefficient = coefficient * lam % r
        assert current == representative and coefficient == 1
    return points, labels, sha(raw)


def coefficients(indices, labels, columns: int, r: int):
    row = [0] * columns
    for index in indices:
        column, coefficient = labels[index]
        row[column] = (row[column] + coefficient) % r
    return row


def replay(ledger: Path, runtime: Path, receipt: Path) -> dict:
    if receipt.exists():
        raise FileExistsError(receipt)
    checked_runtime = json.loads(runtime.read_bytes())
    assert checked_runtime["status"] == "verified"
    assert checked_runtime["sage_version"] == SAGE_VERSION
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    assert plan["parent_freeze_sha256"] == sha(parent_bytes)
    parent = json.loads(parent_bytes)
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    payload = ledger.read_bytes()
    raw = gzip.decompress(payload) if payload[:2] == b"\x1f\x8b" else payload
    rows = [json.loads(line) for line in raw.splitlines()]
    assert len(rows) == 32
    checked = []
    base_receipts = {}
    witness_count = 0
    trace_count = 0
    for n in (41, 53):
        curve, generator, point, r = load_curve(n, parent)
        frozen = plan["curves"][str(n)]["dense"]
        columns = frozen["columns"]
        points, labels, base_sha = load_base(
            n, columns, curve, generator, point, r, frozen["factor_base_digest"])
        assert len(points) == frozen["factor_base_points"]
        field = GF(r)
        traces = {}
        seen_witnesses = set()
        for entry in (row for row in rows if row["n"] == n):
            arm = entry["arm"]
            assert entry["candidate_id"] == plan["curves"][str(n)][arm]["candidate_id"]
            assert entry["verified_native_and_fixture"] and entry["phase_ok"]
            target = targets[(n, entry["label"])]
            assert entry["workload_id"] == target["workload_id"]
            report = entry["report"]
            assert report["target"] == target["public_point"]
            assert report["recovered_scalar"] == target["fixture_scalar"]
            public = point(target["public_point"])
            assert target["fixture_scalar"] * generator == public
            target_indices = tuple(map(int, report["target_relation_indices"]))
            assert len(target_indices) == 4
            assert sum((points[index] for index in target_indices), curve(0)) == public
            target_row = vector(field, coefficients(target_indices, labels, columns, r))
            witnesses = report["rank_relation_witnesses"] + report["fallback_relation_witnesses"]
            trace_hash = sha(canonical(witnesses))
            if trace_hash not in traces:
                relation_rows, rhs = [], []
                for witness in witnesses:
                    indices = tuple(map(int, witness["point_indices"]))
                    scalar = int(witness["relation_scalar"])
                    key = (scalar, indices)
                    assert len(indices) == 4
                    if key not in seen_witnesses:
                        assert sum((points[index] for index in indices), curve(0)) == scalar * generator
                        seen_witnesses.add(key)
                    relation_rows.append(coefficients(indices, labels, columns, r))
                    rhs.append(scalar)
                relation_matrix = matrix(field, relation_rows)
                assert relation_matrix.rank() == report["final_rank"] == columns
                logs = relation_matrix.solve_right(vector(field, rhs))
                traces[trace_hash] = logs
            recovered = int(target_row.dot_product(traces[trace_hash])) % r
            assert recovered == target["fixture_scalar"]
            assert recovered * generator == public
            checked.append({"n": n, "label": entry["label"],
                            "run_id": entry["run_id"], "independent_recovery": True})
        witness_count += len(seen_witnesses)
        trace_count += len(traces)
        base_receipts[str(n)] = {"encoded_base_sha256": base_sha,
                                 "factor_base_points_checked": len(points),
                                 "prepared_relation_traces_checked": len(traces),
                                 "distinct_relation_witnesses_group_checked": len(seen_witnesses)}
    assert len(checked) == len(rows)
    result = {"kind": "s3_dense_root_paired_screen_independent_sage_replay_v1",
              "sage_version": SAGE_VERSION,
              "runtime_info_sha256": sha(runtime.read_bytes()),
              "freeze_sha256": sha(plan_bytes),
              "ledger_sha256": sha(raw),
              "run_count": len(rows),
              "target_count": len({(row["n"], row["label"]) for row in rows}),
              "factor_bases": base_receipts,
              "distinct_relation_witnesses_group_checked": witness_count,
              "distinct_prepared_relation_traces_rank_checked": trace_count,
              "all_scalars_independently_recovered": True,
              "runs": checked}
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"runs": len(rows), "witnesses": witness_count,
                      "traces": trace_count}, sort_keys=True), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--runtime-info", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    replay(args.ledger, args.runtime_info, args.receipt)
