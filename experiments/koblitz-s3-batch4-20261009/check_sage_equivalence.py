"""Match batch4 equations to the independently Sage-replayed prepared pilot."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tarfile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / "experiments/koblitz-s3-crosscurve-20261009/local"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def equations(report: dict) -> list[tuple[int, tuple[int, ...]]]:
    witnesses = report["rank_relation_witnesses"] + report["fallback_relation_witnesses"]
    return [(int(item["relation_scalar"]), tuple(item["point_indices"]))
            for item in witnesses]


def main(ledger: Path, receipt_path: Path) -> None:
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    prior_receipt_bytes = (PRIOR / "sage_replay_pilot24.json").read_bytes()
    prior_receipt = json.loads(prior_receipt_bytes)
    assert prior_receipt["all_scalars_independently_recovered"]
    assert prior_receipt["run_count"] == 192
    with tarfile.open(PRIOR / "raw_pilot24.tar.gz", "r:gz") as archive:
        old_bytes = archive.extractfile("raw.jsonl").read()
    assert sha(old_bytes) == prior_receipt["ledger_sha256"]
    old_rows = [json.loads(line) for line in old_bytes.splitlines()]
    old = {(row["n"], row["label"]): row["report"] for row in old_rows
           if row["arm"] == "prepared" and row["run_id"].endswith("R1")}
    new_bytes = ledger.read_bytes()
    new_rows = [json.loads(line) for line in new_bytes.splitlines()]
    selected = [row for row in new_rows if row["arm"] == "batch4"]
    assert selected and len(selected) * 2 == len(new_rows)
    checked = []
    for row in selected:
        assert row["verified_native_and_fixture"] and row["phase_ok"]
        key = row["n"], row["label"]
        before, after = old[key], row["report"]
        assert before["factor_base_digest"] == after["factor_base_digest"]
        assert before["target"] == after["target"]
        assert before["target_relation_indices"] == after["target_relation_indices"]
        assert before["recovered_scalar"] == after["recovered_scalar"]
        assert equations(before) == equations(after)
        checked.append({"n": row["n"], "label": row["label"],
                        "run_id": row["run_id"],
                        "rank_relation_count": len(after["rank_relation_witnesses"]),
                        "fallback_relation_count": len(after["fallback_relation_witnesses"])})
    assert len({(row["n"], row["label"], row["run_id"]) for row in checked}) == len(checked)
    receipt = {"kind": "s3_batch4_equation_equivalence_to_sage_replayed_pilot_v1",
               "prior_sage_receipt_sha256": sha(prior_receipt_bytes),
               "prior_ledger_sha256": sha(old_bytes),
               "batch4_ledger_sha256": sha(new_bytes),
               "batch4_runs_checked": len(checked),
               "distinct_public_points": len({(row["n"], row["label"]) for row in checked}),
               "same_target_relations_and_rank_equations": True,
               "runs": checked}
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"receipt": str(receipt_path), "runs": len(checked),
                      "points": receipt["distinct_public_points"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger, args.receipt)
