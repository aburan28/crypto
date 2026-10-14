#!/usr/bin/env python3
"""Check exact paired invariants in the frozen ecbench relation sessions."""

import json
import sys
from pathlib import Path


records = [json.loads(line) for line in (Path(sys.argv[1]) / "records.jsonl").read_text().splitlines()]
records = [record for record in records if not record["warmup"]]
by_key = {(record["workload"]["workload_id"], record["round"], record["arm"]): record for record in records}
arms = {record["arm"] for record in records}
keys = {(record["workload"]["workload_id"], record["round"]) for record in records}


def native(record, phase, name):
    return next(entry for entry in record["phases"] if entry["name"] == phase)["native"][name]


checks = {}
if {"b64-full", "b64-pin"} <= arms:
    checks["full_rank_first_pin_matches_pinned_arm"] = all(
        native(by_key[w, r, "b64-full"], "linear_algebra", "first_target_pinned_trial")
        == by_key[w, r, "b64-pin"]["counters"]["targets_tried"]
        and native(by_key[w, r, "b64-full"], "linear_algebra", "first_target_pinned_relation")
        == by_key[w, r, "b64-pin"]["counters"]["relations_found"]
        and by_key[w, r, "b64-full"]["outcome"]["recovered"]
        == by_key[w, r, "b64-pin"]["outcome"]["recovered"]
        for w, r in keys
    )
if {"b64-subtract", "b64-pin"} <= arms:
    checks["subtraction_and_table_identical_stop_and_answer"] = all(
        all(
            by_key[w, r, arm]["counters"][field]
            == by_key[w, r, "b64-pin"]["counters"][field]
            for field in ["targets_tried", "relations_found", "matrix_rank"]
        )
        and by_key[w, r, arm]["outcome"]["recovered"]
        == by_key[w, r, "b64-pin"]["outcome"]["recovered"]
        for w, r in keys
        for arm in ["b64-subtract"]
    )
if {"b48-walk", "b48-native"} <= arms:
    fields = ["targets_tried", "relations_found", "matrix_rank", "matrix_rows", "matrix_work"]
    native_fields = ["walk_steps", "walk_restarts", "repeated_targets_skipped", "lookups", "target_guard_probes"]
    checks["native_walk_identical_relation_stream_and_answer"] = all(
        all(
            by_key[w, r, "b48-walk"]["counters"][field]
            == by_key[w, r, "b48-native"]["counters"][field]
            for field in fields
        )
        and all(
            native(by_key[w, r, "b48-walk"], "relations", field)
            == native(by_key[w, r, "b48-native"], "relations", field)
            for field in native_fields
        )
        and by_key[w, r, "b48-walk"]["outcome"]["recovered"]
        == by_key[w, r, "b48-native"]["outcome"]["recovered"]
        and by_key[w, r, "b48-walk"]["cost"]["total_gae"]
        == by_key[w, r, "b48-native"]["cost"]["total_gae"]
        for w, r in keys
    )
if {"rho-neg", "rho-neg-aa"} <= arms:
    checks["rho_aa_identical_count"] = all(
        by_key[w, r, "rho-neg"]["cost"]["total_gae"]
        == by_key[w, r, "rho-neg-aa"]["cost"]["total_gae"]
        for w, r in keys
    )

result = {"measured_pairs": len(keys), "checks": checks, "all_pass": all(checks.values())}
print(json.dumps(result, indent=2))
if not result["all_pass"]:
    raise SystemExit(1)
