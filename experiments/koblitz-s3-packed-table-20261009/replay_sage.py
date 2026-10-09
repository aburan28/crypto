"""Independently replay packed-table S3 target relations in checked Sage."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

from sage.all import GF, matrix, vector
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
BASE_REPLAY = ROOT / "experiments/koblitz-s3-base-size-20261009/replay_sage.py"
spec = importlib.util.spec_from_file_location("s3_base_size_replay", BASE_REPLAY)
base_replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base_replay)


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def replay(ledger_path: Path, runtime_path: Path, receipt_path: Path) -> dict:
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    runtime = json.loads(runtime_path.read_text())
    assert runtime["status"] == "verified" and runtime["sage_version"] == SAGE_VERSION
    plan_bytes = (HERE / "FREEZE.json").read_bytes()
    plan = json.loads(plan_bytes)
    parent_bytes = PARENT.read_bytes()
    assert plan["parent_plan_sha256"] == sha(parent_bytes)
    targets = {(item["n"], item["label"]): item
               for item in json.loads(parent_bytes)["targets"]}
    payload = ledger_path.read_bytes()
    raw = gzip.decompress(payload) if payload[:2] == b"\x1f\x8b" else payload
    rows = [json.loads(line) for line in raw.splitlines()]
    assert len(rows) == 32 and plan["factor_base_columns"] == 396
    curve, generator, point = base_replay.load_curve()
    candidate = plan["arms"]["packed"]["candidate_id"]
    record = json.loads((HERE / "identities/candidates" / f"{candidate}.json").read_text())
    assert record.pop("candidate_id") == candidate
    base = base_replay.check_base(396, curve, generator, point, record)
    assert base["base_digest"] == plan["factor_base_digest"]
    field = GF(base_replay.R)
    traces = {}
    seen_witnesses = set()
    checked = []
    for entry in rows:
        target = targets[(53, entry["label"])]
        arm = entry["arm"]
        frozen = plan["arms"][arm]
        assert entry["candidate_id"] == frozen["candidate_id"]
        assert entry["run_id"] == (
            f"{frozen['candidate_id']}W{target['workload_id']}R{entry['repetition']}")
        assert entry["source_sha256"] == frozen["source_sha256"]
        assert entry["binary_sha256"] == frozen["binary_sha256"]
        report = entry["report"]
        assert entry["verified_native_and_fixture"] and entry["phase_ok"]
        assert report["group_verified"] and report["target_count"] == 1
        assert report["factor_base_digest"] == base["base_digest"]
        assert report["factor_base_points"] == len(base["points"])
        assert report["target"] == target["public_point"]
        assert report["recovered_scalar"] == target["fixture_scalar"]
        public_point = point(target["public_point"])
        assert int(target["fixture_scalar"]) * generator == public_point
        indices = tuple(map(int, report["target_relation_indices"]))
        assert len(indices) == 4
        assert sum((base["points"][index] for index in indices), curve(0)) == public_point
        target_row = vector(field, base_replay.row_coefficients(
            indices, base["labels"], 396))
        witnesses = report["rank_relation_witnesses"] + report["fallback_relation_witnesses"]
        trace_hash = sha(json.dumps(witnesses, sort_keys=True, separators=(",", ":")).encode())
        if trace_hash not in traces:
            coefficients = []
            rhs = []
            for witness in witnesses:
                witness_indices = tuple(map(int, witness["point_indices"]))
                scalar = int(witness["relation_scalar"])
                assert len(witness_indices) == 4
                key = (scalar, witness_indices)
                if key not in seen_witnesses:
                    assert sum((base["points"][index] for index in witness_indices),
                               curve(0)) == scalar * generator
                    seen_witnesses.add(key)
                coefficients.append(base_replay.row_coefficients(
                    witness_indices, base["labels"], 396))
                rhs.append(scalar)
            relation_matrix = matrix(field, coefficients)
            relation_rhs = vector(field, rhs)
            assert relation_matrix.rank() == report["final_rank"] == 396
            logs = relation_matrix.solve_right(relation_rhs)
            traces[trace_hash] = logs
        recovered = int(target_row.dot_product(traces[trace_hash])) % base_replay.R
        assert recovered == target["fixture_scalar"]
        assert recovered * generator == public_point
        checked.append({"label": entry["label"], "arm": arm,
                        "run_id": entry["run_id"], "rank": 396,
                        "independent_recovery": True})
    assert len(checked) == 32 and len(traces) == 1
    result = {"kind": "s3_packed_root_table_independent_sage_replay_v1",
              "sage_version": SAGE_VERSION,
              "runtime_info_sha256": sha(runtime_path.read_bytes()),
              "plan_sha256": sha(plan_bytes),
              "parent_plan_sha256": sha(parent_bytes),
              "ledger_sha256": sha(raw),
              "base_sha256": base["base_sha256"],
              "factor_base_points_and_labels_checked": base["base_point_count_checked"],
              "representatives_subgroup_checked": base["representatives_subgroup_checked"],
              "distinct_relation_witnesses_group_checked": len(seen_witnesses),
              "distinct_relation_traces_rank_checked": len(traces),
              "run_count": len(checked), "target_count": 8,
              "all_scalars_independently_recovered": True,
              "runs": checked}
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(receipt_path), "runs": len(checked),
                      "witnesses": len(seen_witnesses)}), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=HERE / "local/raw_paired32.jsonl.gz")
    parser.add_argument("--runtime-info", type=Path,
                        default=HERE / "local/sage_runtime_info.json")
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    replay(args.ledger, args.runtime_info, args.receipt)
