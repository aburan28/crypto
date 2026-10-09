"""Check that tiled queries choose the same earliest relation as serial S3."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
GROUPLIFT = ROOT / "experiments/koblitz-s3-grouplift-20261009/FREEZE.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(ledger: Path, binary: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    parent = json.loads(PARENT.read_bytes())
    old = json.loads(GROUPLIFT.read_bytes())
    plan = json.loads((HERE / "FREEZE.json").read_bytes())
    assert old["parent_plan_sha256"] == sha(PARENT)
    assert plan["parent_plan_sha256"] == sha(PARENT)
    assert sha(binary) == old["local_binary_sha256"]
    targets = {(item["n"], item["label"]): item for item in parent["targets"]}
    rows = [json.loads(line) for line in ledger.read_bytes().splitlines()]
    by_label = defaultdict(list)
    for row in rows:
        by_label[row["label"]].append(row)
    controls = []
    scratch = output.parent / "serial_control_raw"
    scratch.mkdir(parents=True, exist_ok=True)
    for label, selected in sorted(by_label.items()):
        target = targets[53, label]
        point_file = PARENT.parent / target["point_file_relative"]
        raw_file = scratch / f"n53_{label}.json"
        argv = [str(binary), "53", "0", "244", str(parent["relation_seed"]),
                str(point_file), str(raw_file), "14", "243"]
        run = subprocess.run(argv, capture_output=True, timeout=180)
        assert run.returncode == 0
        raw = raw_file.read_bytes()
        assert run.stdout == raw
        report = json.loads(raw)
        assert report["group_verified"]
        assert report["recovered_scalar"] == target["fixture_scalar"]
        assert report["target"] == target["public_point"]
        expected = report["target_relation_indices"]
        tile_rows = [row for row in selected if row["tile_states"] is not None]
        assert len(tile_rows) == 6 and all(row["verified_native_and_fixture"] for row in tile_rows)
        assert all(row["report"]["target_relation_indices"] == expected
                   for row in tile_rows)
        controls.append({"label": label, "workload_id": target["workload_id"],
                         "point_sha256": target["point_sha256"],
                         "serial_candidate_id": old["curves"]["53"]["grouplift_candidate_id"],
                         "serial_raw_sha256": hashlib.sha256(raw).hexdigest(),
                         "serial_first_relation_indices": expected,
                         "tiled_runs_matching_serial_first_relation": len(tile_rows)})
        print(json.dumps({"label": label, "matched_tiled_runs": len(tile_rows)}), flush=True)
    receipt = {"kind": "s3_tiled_serial_first_relation_control_v1",
               "plan_sha256": sha(HERE / "FREEZE.json"),
               "ledger_sha256": sha(ledger),
               "serial_binary_sha256": sha(binary),
               "target_count": len(controls), "all_tiled_runs_match": True,
               "controls": controls}
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger, args.binary, args.output)
