"""Assign an IC1 identity to the verified C1800 row-elision pipeline."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/identities/candidates/IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h1ab23390e700.json"
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan_norow.rs"
CONTROL = HERE / "local/control/receipt.json"
CONTROL_RAW = HERE / "local/control/report.json.gz"
OUTPUT = HERE / "identities/candidates"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def replace_digest(value: object, old: str, new: str) -> object:
    if isinstance(value, dict):
        return {key: replace_digest(item, old, new) for key, item in value.items()}
    if isinstance(value, list):
        return [replace_digest(item, old, new) for item in value]
    return new if value == old else value


def main() -> None:
    old = json.loads(OLD.read_bytes())
    old_id = old.pop("candidate_id")
    assert old_id.endswith("h" + sha(canonical(old))[:12])
    old_source_sha = old["implementation"]["example_source_sha256"]
    new_source_sha = sha(SOURCE.read_bytes())
    source_relative = str(SOURCE.relative_to(ROOT))
    assert sha(subprocess.check_output(["git", "show", "HEAD:" + source_relative],
                                       cwd=ROOT)) == new_source_sha
    receipt = json.loads(CONTROL.read_bytes())
    raw = gzip.decompress(CONTROL_RAW.read_bytes())
    report = json.loads(raw)
    assert receipt["scalar_verified"] and receipt["full_rank_before_target"]
    assert receipt["memory_within_envelope"]
    assert receipt["source_sha256"] == new_source_sha
    assert receipt["raw_report_sha256"] == sha(raw)
    assert report["target_query_batch_width"] == 8
    assert report["target_full_rank_row_policy"] == (
        "skip_row_when_prepared_base_logs_available")
    assert report["factor_base_points"] == 190800
    assert report["factor_base_digest"] == old["factor_base"]["enumerated_set_digest"]
    record = replace_digest(old, old_source_sha, new_source_sha)
    record["implementation"]["frozen_source_commit"] = git("rev-parse", "HEAD")
    assert record["implementation"]["rust_library_source_tree_git_oid"] == git(
        "rev-parse", "HEAD:src")
    flags = record["implementation"]["algorithm_flags"]
    assert flags["target_query_batch_width"] == 8
    flags["target_full_rank_row_policy"] = (
        "skip_row_when_prepared_base_logs_available")
    record["target_descent"]["target_full_rank_row_policy"] = (
        "skip_row_when_prepared_base_logs_available")
    assert record["factor_base"]["actual_usable_points_B"] == 190800
    prefix = "IC1N53Ckb1fb190800PDP4rootRCsampleLAgaussTDdirectISO0h"
    candidate_id = prefix + sha(canonical(record))[:12]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT / (candidate_id + ".json")
    if path.exists():
        raise FileExistsError(path)
    (path).write_bytes(canonical({**record, "candidate_id": candidate_id}) + b"\n")
    print(json.dumps({"candidate_id": candidate_id,
                      "manifest_sha256": sha(path.read_bytes()),
                      "source_sha256": new_source_sha}, sort_keys=True))


if __name__ == "__main__":
    main()
