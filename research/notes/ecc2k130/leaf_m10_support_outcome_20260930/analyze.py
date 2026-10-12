#!/usr/bin/env python3
"""Derive exact capacity diagnostics from the archived m10 support census."""
from __future__ import annotations

from decimal import Decimal, getcontext
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent / "leaf_m10_support_20260930/evidence"
CONFIG = HERE.parent / "leaf_m10_support_20260930/CONFIG.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ratio(numerator: int, denominator: int) -> str:
    return format(Decimal(numerator) / Decimal(denominator), ".12f")


def main() -> None:
    getcontext().prec = 80
    config = json.loads(CONFIG.read_text())
    result = json.loads((EVIDENCE / "result.json").read_text())
    hosted = json.loads((EVIDENCE / "replay.json").read_text())
    local = json.loads((EVIDENCE / "mac_replay.json").read_text())
    assert result["status"] == "PASS_CENSUS"
    assert result["config_sha256"] == sha(CONFIG)
    assert hosted["status"] == local["status"] == "PASS"
    assert hosted == local
    assert hosted["result_sha256"] == sha(EVIDENCE / "result.json")
    assert hosted["scans_replayed"] == 33
    assert hosted["raw_masks_replayed"] == 294912
    assert hosted["sample_point_images_checked"] == 264
    assert set(result["scans"]) == {curve["id"] for curve in config["curves"]}

    source = result["scans"]["source"]
    slots = [f"low_{i}" for i in range(10)] + ["high_0"]
    assert all(source[slot]["physical_points"] == 7977 for slot in slots[:10])
    assert source["high_0"]["physical_points"] == 16125
    q = int(config["subgroup_order"])
    rows = []
    for curve in config["curves"]:
        curve_id = curve["id"]
        scans = result["scans"][curve_id]
        assert set(scans) == set(slots)
        counts = [scans[slot]["physical_points"] for slot in slots[:10]]
        classes = [scans[slot]["projected_sign_classes"] for slot in slots[:10]]
        assert all(all(scans[slot]["column_x_multiplicity_histogram"][str(k)] == 0
                       for k in (2, 3, 4)) for slot in slots)
        arm_rows = {}
        for arm in config["arms"]:
            arm_id = arm["id"]
            record = result["arms"][curve_id][arm_id]
            source_record = result["arms"]["source"][arm_id]
            product = int(record["physical_tuple_product"])
            source_product = int(source_record["physical_tuple_product"])
            selected = record["slots"]
            assert record["physical_point_counts"] == [
                scans[slot]["physical_points"] for slot in selected]
            assert record["uncompressed_projected_sign_union"] == sum(
                scans[slot]["projected_sign_classes"] for slot in selected)
            arm_rows[arm_id] = {
                "physical_tuple_product": str(product),
                "necessary_tuple_capacity_over_q": ratio(product, q),
                "physical_tuple_ratio_to_source": ratio(product, source_product),
                "uncompressed_projected_sign_union": record[
                    "uncompressed_projected_sign_union"],
                "projected_union_ratio_to_source": ratio(
                    record["uncompressed_projected_sign_union"],
                    source_record["uncompressed_projected_sign_union"]),
                "projected_columns_disjoint_across_slots": True,
            }
        rows.append({
            "curve": curve_id,
            "low_slot_physical_points": counts,
            "low_slot_projected_sign_classes": classes,
            "high_slot_physical_points": scans["high_0"]["physical_points"],
            "high_slot_projected_sign_classes": scans["high_0"][
                "projected_sign_classes"],
            "low_slot_physical_min": min(counts),
            "low_slot_physical_max": max(counts),
            "low_slot_counts_constant": len(set(counts)) == 1,
            "source_low_slot_count_transfer_holds": all(
                count == 7977 for count in counts),
            "arms": arm_rows,
        })
    analysis = {
        "schema": "ecc2k130-leaf-m10-support-analysis-v1",
        "status": "PASS_DIAGNOSTIC",
        "result_sha256": sha(EVIDENCE / "result.json"),
        "hosted_replay_sha256": sha(EVIDENCE / "replay.json"),
        "mac_replay_sha256": sha(EVIDENCE / "mac_replay.json"),
        "raw_masks_replayed_on_each_host": 294912,
        "sample_point_images_checked_on_each_host": 264,
        "rows": rows,
        "PDP_yield": None,
        "full_ECDLP_cost": None,
        "method_crossover": None,
    }
    (HERE / "ANALYSIS.json").write_text(json.dumps(
        analysis, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": analysis["status"],
                      "leaf_count_transfer": {row["curve"]: row[
                          "source_low_slot_count_transfer_holds"]
                          for row in rows[1:]}}, sort_keys=True))


if __name__ == "__main__":
    main()
