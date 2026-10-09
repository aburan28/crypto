"""Attach arm labels from immutable candidate IDs to completed paired rows."""

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
    by_id = {candidate: arm for arm, candidate in
             plan["candidate_ids_by_arm"].items()}
    assert len(by_id) == 3
    source = input_path.read_bytes()
    rows = [json.loads(line) for line in source.splitlines()]
    assert len(rows) == 6 * len(plan["targets"])
    seen = set()
    for row in rows:
        arm = by_id[row["candidate_id"]]
        assert row.get("arm", arm) == arm
        row["arm"] = arm
        assert row["run_id"] == (
            f"{row['candidate_id']}W{row['workload_id']}R{row['repetition']}")
        assert row["run_id"] not in seen
        seen.add(row["run_id"])
    output_path.write_text("".join(json.dumps(row, separators=(",", ":")) + "\n"
                                for row in rows))
    print(json.dumps({"rows": len(rows), "output": str(output_path),
                      "source_sha256": sha(source),
                      "canonical_sha256": sha(output_path.read_bytes())}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.input, args.output)
