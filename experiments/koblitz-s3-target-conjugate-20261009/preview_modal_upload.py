"""Freeze the exact local files a C1200 Modal VM panel would receive.

This command only reads local files and writes a local manifest. It does not
connect to Modal or transfer any repository content.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
HOLDOUT = HERE / "MODAL_HOLDOUT_FREEZE.json"
PILOT = HERE / "C1200_FULLRANK_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
EXAMPLES = {
    "baseline": ROOT / "examples/s3_partial_prepared_parallel_targetspan.rs",
    "expanded": ROOT / "examples/s3_partial_prepared_fullrank_targetspan.rs",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def upload_files() -> list[Path]:
    holdout = json.loads(HOLDOUT.read_text())
    directories = [ROOT / "src"]
    files = [ROOT / "Cargo.toml", ROOT / "Cargo.lock", *EXAMPLES.values(),
             HOLDOUT, PILOT, PARENT, HERE / "modal_c1200_fullrank_panel.py",
             ROOT / "experiments/koblitz-s3-prepared-20261009/modal_prepared_panel.py",
             *(ROOT / name for name in (
                 "docs/ic/boundary_targets.json",
                 "docs/ic/calibration.json",
                 "docs/ecbench/schema.sql",
                 "docs/curves/registry.json"))]
    files.extend(HERE / "identities/candidates" / f"{arm['candidate_id']}.json"
                 for arm in holdout["arms"].values())
    files.extend(ROOT / "experiments/koblitz-s3-crosscurve-20261009" /
                 target["point_file_relative"] for target in holdout["targets"])
    warmup = next(target for target in json.loads(PARENT.read_text())["targets"]
                  if target["n"] == 53 and target["label"] == "T025")
    files.append(ROOT / "experiments/koblitz-s3-crosscurve-20261009" /
                 warmup["point_file_relative"])
    for directory in directories:
        files.extend(path for path in directory.rglob("*") if path.is_file())
    assert len(files) == len(set(files))
    assert all(path.is_file() for path in files)
    return sorted(files)


def preview() -> dict:
    holdout = json.loads(HOLDOUT.read_text())
    pilot = json.loads(PILOT.read_text())
    parent = json.loads(PARENT.read_text())
    assert holdout["kind"] == "s3_fullrank_modal_holdout_v1"
    assert holdout["pilot_plan_sha256"] == sha(PILOT.read_bytes())
    assert holdout["parent_plan_sha256"] == sha(PARENT.read_bytes())
    assert holdout["arms"] == pilot["arms"]
    targets = {(row["n"], row["label"]): row for row in parent["targets"]}
    assert len(holdout["targets"]) == 24
    assert [row["label"] for row in holdout["targets"]] == [
        f"T{i:03}" for i in range(33, 57)]
    for arm, source in EXAMPLES.items():
        assert sha(source.read_bytes()) == holdout["arms"][arm]["source_sha256"]
    for row in holdout["targets"]:
        parent_row = targets[(53, row["label"])]
        assert all(row[key] == parent_row[key] for key in
                   ("workload_id", "point_sha256", "public_point", "fixture_scalar"))
        point = (ROOT / "experiments/koblitz-s3-crosscurve-20261009"
                 / row["point_file_relative"])
        assert sha(point.read_bytes()) == row["point_sha256"]
    warmup = targets[(53, "T025")]
    warmup_point = (ROOT / "experiments/koblitz-s3-crosscurve-20261009"
                    / warmup["point_file_relative"])
    assert sha(warmup_point.read_bytes()) == warmup["point_sha256"]
    files = upload_files()
    entries = [{"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
                "sha256": sha(path.read_bytes())} for path in files]
    document = {"kind": "s3_fullrank_modal_upload_manifest_v1",
                "holdout_freeze_sha256": sha(HOLDOUT.read_bytes()),
                "pilot_freeze_sha256": sha(PILOT.read_bytes()),
                "source_tree_git_oid": holdout["source_tree_git_oid"],
                "target_count": len(holdout["targets"]),
                "warmup_label": "T025",
                "files": entries}
    return document


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    document = preview()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical(document) + b"\n")
    print(json.dumps({"manifest": str(args.output),
                      "manifest_sha256": sha(args.output.read_bytes()),
                      "files": len(document["files"]),
                      "bytes": sum(row["bytes"] for row in document["files"]),
                      "public_targets": document["target_count"]},
                     sort_keys=True))
