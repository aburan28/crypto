#!/usr/bin/env python3
"""Fail-closed audit of retained N83 structured-orbit receipts and budget."""
import hashlib
import json
from decimal import Decimal
from pathlib import Path

root = Path(__file__).resolve().parent
def read(name):
    return json.loads((root / name).read_text())

manifest = read("manifest.json")
replay = read("replay.json")
s3_replay = read("s3-replay.json")
upload = read("upload.json")
outer = read("outer.json")
tamper = read("tamper-control.json")
prior = read("prior-k2048-domain.json")
budget = read("budget.json")

assert manifest["schema"] == "n83.structured-orbit-panel/v1"
assert manifest["status"] == "completed_factor_base_object"
assert manifest["source_commit"] == "67b251de4c921a8dd2ba34ad52f6afbc8bdabc3c"
assert manifest["orbit_columns"] == replay["representatives_checked"] == s3_replay["representatives_checked"] == 2001
assert manifest["point_records"] == replay["points_checked"] == s3_replay["points_checked"] == 332166
assert manifest["historical_set_blake3"] == replay["historical_set_blake3"] == s3_replay["historical_set_blake3"] == "db3e7f25877c75dfd6f53db46447972ce6a44434cb7ea1e53e450d21af0fe79c"
assert replay["status"] == s3_replay["status"] == upload["status"] == outer["status"] == tamper["status"] == "PASS"
assert upload["object"] == manifest["object"] == replay["object"] == s3_replay["object"]
assert upload["object_blake3"] == manifest["compressed_blake3"] == replay["compressed_blake3"] == s3_replay["compressed_blake3"]
assert upload["object_s3_uri"] == manifest["s3_uri"]
assert upload["manifest_sha256"] == hashlib.sha256((root / "manifest.json").read_bytes()).hexdigest()
assert upload["s3_head"]["ContentLength"] == manifest["compressed_bytes"] == 10058074
assert all(stage["status"] == "PASS" and stage["exit_code"] == 0 for stage in outer["stages"])
assert next(stage for stage in upload["stages"] if stage["name"] == "s3-head-before")["exit_code"] == 254
assert "(404)" in (root / "s3-head-before.stderr").read_text()
assert all(stage.get("exit_code") == 0 for stage in upload["stages"] if "exit_code" in stage and stage["name"] != "s3-head-before")
assert tamper["exit_code"] != 0
assert prior["m5_required_clauses"] == 55017910 > prior["frozen_clause_cap"]
assert Decimal(prior["prior_charged_seconds"]) + Decimal(prior["round_conservative_charge_seconds"]) == Decimal(prior["cumulative_charged_seconds"])
assert Decimal(budget["prior_charged_seconds"]) == Decimal(prior["cumulative_charged_seconds"])
assert Decimal(budget["measured_active_wall_seconds"]) == Decimal(str(outer["active_wall_seconds"])) + Decimal(str(upload["active_wall_seconds"]))
assert Decimal(budget["measured_active_wall_seconds"]) < Decimal(budget["round_conservative_charge_seconds"])
assert Decimal(budget["prior_charged_seconds"]) + Decimal(budget["round_conservative_charge_seconds"]) == Decimal(budget["cumulative_charged_seconds"])
assert Decimal(budget["cumulative_charged_seconds"]) + Decimal(budget["remaining_unused_seconds"]) == Decimal(budget["authorized_seconds"]) == 7200
assert manifest["relation_yield"] is manifest["matrix_rank"] is manifest["total_index_calculus_runtime_ms"] is None
assert replay["relation_stage_executed"] is replay["rank_stage_executed"] is False

print(json.dumps({
    "schema": "n83.structured-orbit-evidence-check/v1",
    "status": "PASS",
    "object": manifest["object"],
    "points_checked": replay["points_checked"],
    "s3_points_replayed": s3_replay["points_checked"],
    "cumulative_charged_seconds": budget["cumulative_charged_seconds"],
    "remaining_unused_seconds": budget["remaining_unused_seconds"],
}, sort_keys=True))
