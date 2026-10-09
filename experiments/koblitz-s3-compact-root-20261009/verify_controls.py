"""Check that 13-byte root slots preserve the full-rank S3 control traces."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_compact_targetspan.rs"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
N53_REFERENCE = (ROOT / "experiments/koblitz-s3-target-conjugate-20261009/local/"
                 "raw_c1200_fullrank_paired32.jsonl.gz")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def load_gzip_json(path: Path) -> tuple[dict, str]:
    raw = gzip.decompress(path.read_bytes())
    return json.loads(raw), sha(raw)


def verify(binary: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    source_sha = sha(SOURCE.read_bytes())
    binary_sha = sha(binary.read_bytes())
    parent_bytes = PARENT.read_bytes()
    targets = {(row["n"], row["label"]): row
               for row in json.loads(parent_bytes)["targets"]}
    assert (41, "T025") in targets and (53, "T025") in targets
    old_n41, old_n41_sha = load_gzip_json(HERE / "local/reference_n41_t025.json.gz")
    old_n53_rows = [json.loads(line) for line in gzip.decompress(N53_REFERENCE.read_bytes()).splitlines()]
    old_n53_matches = [row for row in old_n53_rows if row["n"] == 53 and
                       row["label"] == "T025" and row["arm"] == "expanded" and
                       row["repetition"] == 1]
    assert len(old_n53_matches) == 1
    old_n53 = old_n53_matches[0]["report"]
    result = {
        "kind": "s3_packed_root_slot_fullrank_correctness_controls_v1",
        "source_sha256": source_sha,
        "binary_sha256": binary_sha,
        "parent_freeze_sha256": sha(parent_bytes),
        "n53_reference_ledger_sha256": sha(gzip.decompress(N53_REFERENCE.read_bytes())),
        "controls": {},
        "wall_speedup": None,
    }
    for n, columns, old, old_sha in ((41, 800, old_n41, old_n41_sha),
                                     (53, 1200, old_n53, sha(canonical(old_n53)))):
        new, new_sha = load_gzip_json(HERE / f"local/compact_n{n}_t025.json.gz")
        matches = list((HERE / "identities/candidates").glob(f"IC1N{n}C*PDP4root*.json"))
        assert len(matches) == 1
        candidate = json.loads(matches[0].read_bytes())
        candidate_id = candidate.pop("candidate_id")
        assert candidate_id.endswith("h" + sha(canonical(candidate))[:12])
        assert candidate["point_decomposition"]["implementation_source_sha256"] == source_sha
        assert candidate["implementation"]["example_source_sha256"] == source_sha
        assert candidate["factor_base"]["actual_usable_points_B"] == new["factor_base_points"]
        assert candidate["factor_base"]["effective_relation_columns"] == columns
        assert candidate["factor_base"]["enumerated_set_digest"] == new["factor_base_digest"]
        assert candidate["point_decomposition"]["root_table_slot_bytes"] == 13
        assert candidate["point_decomposition"]["root_table_value_bits"] == 40
        assert candidate["point_decomposition"]["root_table_payload_bytes"] == 13 * new["root_table_slots"]
        assert new["root_table_slot_bytes"] == 13 and new["root_table_value_bits"] == 40
        assert new["root_table_slots"] == old["root_table_slots"]
        assert new["group_verified"] and new["prepared_rank"] == columns
        assert new["fallback_attempts"] == 0 and new["target_full_rank_direct_recovery"]
        target = targets[(n, "T025")]
        assert new["target"] == target["public_point"]
        assert new["recovered_scalar"] == target["fixture_scalar"]
        non_timing_old = {k: v for k, v in old.items() if k != "timing_ms"}
        non_timing_new = {k: v for k, v in new.items() if k not in (
            "timing_ms", "root_table_slot_bytes", "root_table_value_bits")}
        assert non_timing_old == non_timing_new
        result["controls"][f"N{n}"] = {
            "candidate_id": candidate_id,
            "target_label": "T025",
            "point_sha256": target["point_sha256"],
            "reference_report_sha256": old_sha,
            "compact_report_sha256": new_sha,
            "non_timing_reports_equal": True,
            "prepared_rank": columns,
            "root_table_slots": new["root_table_slots"],
            "reference_root_slot_bytes": 16,
            "compact_root_slot_bytes": 13,
            "root_payload_bytes_before": 16 * new["root_table_slots"],
            "root_payload_bytes_after": 13 * new["root_table_slots"],
            "root_payload_bytes_saved": 3 * new["root_table_slots"],
            "target_s3_calls": new["target_s3_calls"],
            "recovered_scalar": new["recovered_scalar"],
        }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"binary_sha256": binary_sha,
                      "saved_bytes": {key: value["root_payload_bytes_saved"]
                                      for key, value in result["controls"].items()}}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    verify(args.binary, args.output)
