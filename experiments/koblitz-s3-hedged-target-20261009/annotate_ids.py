"""Attach frozen IC1 candidate and run IDs to completed tuning probe rows."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main(input_path: Path, output_path: Path) -> None:
    if output_path.exists():
        raise FileExistsError(output_path)
    plan = json.loads((HERE / "FREEZE.json").read_bytes())
    parent = json.loads((HERE.parents[1] /
                         "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json").read_bytes())
    assert plan["parent_plan_sha256"] == sha((HERE.parents[1] /
        "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json").read_bytes())
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    rows = [json.loads(line) for line in input_path.read_bytes().splitlines()]
    assert len(rows) == 8 * plan["target_count"]
    seen = set()
    for row in rows:
        assert row["source_sha256"] == plan["source_sha256"]
        assert row["hedged_binary_sha256"] == plan["local_binary_sha256"]
        target = targets[row["n"], row["label"]]
        assert row["target_sha256"] == target["point_sha256"]
        assert row["workload_id"] == target["workload_id"]
        candidate = (plan["reference_candidate_id"] if row["delay_us"] is None
                     else plan["hedged_candidate_ids_by_delay_us"][str(row["delay_us"])])
        run_id = f"{candidate}W{row['workload_id']}R{row['repetition']}"
        assert run_id not in seen
        seen.add(run_id)
        row["candidate_id"] = candidate
        row["run_id"] = run_id
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(json.dumps(row, separators=(",", ":")) + "\n"
                                   for row in rows))
    print(json.dumps({"output": str(output_path), "rows": len(rows),
                      "source_ledger_sha256": sha(input_path.read_bytes()),
                      "canonical_ledger_sha256": sha(output_path.read_bytes())}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.input, args.output)
