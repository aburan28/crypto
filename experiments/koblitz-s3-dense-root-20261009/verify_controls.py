"""Verify dense root-table capacity and unchanged N41/N53 S3 traces."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "examples/s3_partial_prepared_fullrank_dense_targetspan.rs"
COMPACT = ROOT / "experiments/koblitz-s3-compact-root-20261009/local"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
POLICY = "next_power_of_two_of_input_roots_plus_one_per_shard"
CAPACITY_FIELDS = {"root_table_capacity_policy", "root_table_slots",
                   "root_table_slots_per_shard"}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def load(path: Path) -> tuple[dict, str]:
    raw = gzip.decompress(path.read_bytes())
    return json.loads(raw), sha(raw)


def verify(binary: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    source_sha = sha(SOURCE.read_bytes())
    parent_bytes = PARENT.read_bytes()
    targets = {(row["n"], row["label"]): row
               for row in json.loads(parent_bytes)["targets"]}
    result = {"kind": "s3_dense_packed_root_slot_fullrank_controls_v1",
              "source_sha256": source_sha, "binary_sha256": sha(binary.read_bytes()),
              "parent_freeze_sha256": sha(parent_bytes),
              "controls": {}, "wall_speedup": None}
    for n, columns in ((41, 800), (53, 1200)):
        old, old_sha = load(COMPACT / f"compact_n{n}_t025.json.gz")
        new, new_sha = load(HERE / f"local/dense_n{n}_t025.json.gz")
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
        assert candidate["point_decomposition"]["root_table_slots"] == new["root_table_slots"]
        assert candidate["point_decomposition"]["root_table_slot_bytes"] == 13
        assert candidate["point_decomposition"]["root_table_value_bits"] == 40
        assert candidate["point_decomposition"]["root_table_capacity_policy"] == POLICY
        assert new["root_table_capacity_policy"] == POLICY
        assert new["root_table_slot_bytes"] == 13 and new["root_table_value_bits"] == 40
        inputs = new["root_table_input_entries_per_shard"]
        capacities = new["root_table_slots_per_shard"]
        assert len(inputs) == len(capacities) == 64
        assert all(slots == max(16, 1 << entries.bit_length()) and slots > entries
                   for entries, slots in zip(inputs, capacities))
        assert sum(capacities) == new["root_table_slots"]
        assert new["root_table_entries"] <= sum(inputs) < new["root_table_slots"]
        assert new["group_verified"] and new["prepared_rank"] == columns
        assert new["fallback_attempts"] == 0 and new["target_full_rank_direct_recovery"]
        target = targets[(n, "T025")]
        assert new["target"] == target["public_point"]
        assert new["recovered_scalar"] == target["fixture_scalar"]
        old_semantic = {k: v for k, v in old.items() if k != "timing_ms" and
                        k not in CAPACITY_FIELDS}
        new_semantic = {k: v for k, v in new.items() if k != "timing_ms" and
                        k not in CAPACITY_FIELDS}
        assert old_semantic == new_semantic
        result["controls"][f"N{n}"] = {
            "candidate_id": candidate_id, "target_label": "T025",
            "point_sha256": target["point_sha256"],
            "compact_report_sha256": old_sha, "dense_report_sha256": new_sha,
            "semantic_reports_equal": True,
            "prepared_rank": columns,
            "root_table_entries": new["root_table_entries"],
            "compact_slots": old["root_table_slots"],
            "dense_slots": new["root_table_slots"],
            "dense_load_fraction": new["root_table_entries"] / new["root_table_slots"],
            "dense_payload_bytes": 13 * new["root_table_slots"],
            "payload_saved_vs_compact_bytes": 13 * (old["root_table_slots"] -
                                                    new["root_table_slots"]),
            "payload_saved_vs_original_16byte_bytes": (16 * old["root_table_slots"] -
                                                       13 * new["root_table_slots"]),
            "target_s3_calls": new["target_s3_calls"],
            "recovered_scalar": new["recovered_scalar"],
        }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value["payload_saved_vs_original_16byte_bytes"]
                      for key, value in result["controls"].items()}, sort_keys=True))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    verify(args.binary, args.output)
