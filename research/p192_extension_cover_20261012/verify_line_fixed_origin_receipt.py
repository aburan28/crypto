"""Custody, class identity, and nonspeciality audit for 141 reduced classes."""

import hashlib
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
eligibility_path = HERE / "type_i_line_fiber_eligibility.json"
cubic_path = HERE / "type_i_line_quotient_cubic_probe.json"
reduced_path = HERE / "type_i_line_fixed_origin_reduction.json"
native_input_path = HERE / "type_i_line_fixed_origin_native_input.json"
image_path = HERE / "type_i_line_quotient_image_verification.json"
gp_path = HERE / "type_i_line_fixed_origin_native_input.gp"
gp_log_path = HERE / "type_i_line_fixed_origin_native_verification.log"
eligibility = json.loads(eligibility_path.read_text())
cubic = json.loads(cubic_path.read_text())
reduced = json.loads(reduced_path.read_text())
native_input = json.loads(native_input_path.read_text())
image = json.loads(image_path.read_text())
assert eligibility["input_classes"] == cubic["classes"] == reduced["input_classes"] == 141
assert sum(eligibility["cells"].values()) == 141
assert cubic["successes"] == reduced["reduced_classes"] == image["classes"] == 141
assert cubic["failures"] == reduced["unreduced_classes"] == 0
assert cubic["source_sha256"] == sha(HERE / "line_quotient_cubic_probe.py")
assert cubic["quotient_receipt_sha256"] == sha(HERE / "type_i_line_guided_collector.json")
assert cubic["eligibility_receipt_sha256"] == sha(eligibility_path)
assert cubic["plane_receipt_sha256"] == sha(HERE / "type_i_plane_quartic.json")
assert cubic["runtime_sha256"] == sha(HERE / "type_i_line_cubic_runtime_info.json")
assert eligibility["source_sha256"] == sha(HERE / "line_quotient_fiber_eligibility.py")
assert eligibility["runtime_sha256"] == sha(HERE / "type_i_line_fiber_runtime_info.json")
assert reduced["source_sha256"] == sha(HERE / "line_quotient_fixed_origin.py")
assert reduced["input_receipt_sha256"] == sha(cubic_path)
assert reduced["plane_receipt_sha256"] == sha(HERE / "type_i_plane_quartic.json")
assert reduced["runtime_sha256"] == sha(HERE / "type_i_line_fixed_origin_runtime_info.json")
assert native_input["status"] == "prepared_line_fixed_origin_native_replay"
assert native_input["cases"] == 141
assert native_input["source_sha256"] == sha(HERE / "prepare_line_fixed_origin_gp.py")
assert native_input["cubic_receipt_sha256"] == sha(cubic_path)
assert native_input["reduction_receipt_sha256"] == sha(reduced_path)
assert native_input["plane_receipt_sha256"] == sha(HERE / "type_i_plane_quartic.json")
assert native_input["gp_input_sha256"] == sha(gp_path)
assert image["status"] == "LINE_QUOTIENT_IMAGES_INDEPENDENT_PASS"
assert image["source_sha256"] == sha(HERE / "verify_line_quotient_images.py")
assert image["quotient_receipt_sha256"] == sha(HERE / "type_i_line_guided_collector.json")
assert image["cubic_receipt_sha256"] == sha(cubic_path)
assert image["reduction_receipt_sha256"] == sha(reduced_path)
assert image["runtime_sha256"] == sha(HERE / "type_i_line_quotient_image_runtime_info.json")
gp_log = gp_log_path.read_text()
assert "LINE_FIXED_ORIGIN_NATIVE_PASS cases=141 prime_bits=192" in gp_log
assert "***" not in gp_log and "error" not in gp_log.lower()
assert gp_log.count("LINE_FIXED_ORIGIN_NATIVE_PROGRESS") == 8
assert len(cubic["rows"]) == len(reduced["representatives"]) == 141

cells = Counter()
for index, (input_row, result) in enumerate(zip(cubic["rows"], reduced["representatives"])):
    assert input_row["index"] == result["index"] == index
    assert input_row["cell"] == result["cell"] == eligibility["rows"][index]["cell"]
    assert input_row["quotient"] == result["quotient"]
    assert input_row["status"] == "verified_prime_field_cubic_pair"
    assert result["status"] == "verified_fixed_origin_degree_three"
    hx, hy = result["representative_x_polynomial"], result["representative_y_polynomial"]
    assert len(hx) == 4 and int(hx[3][0]) == 1
    assert len(hy) == 3 and int(hy[2][0]) != 0
    cells[result["cell"]] += 1
assert dict(cells) == {k: v for k, v in eligibility["cells"].items() if v}
assert cubic["successes_by_cell"] == eligibility["cells"]
collector_path = HERE / "type_i_line_all_cell_collector.json"
collector = json.loads(collector_path.read_text())
assert collector["status"] == "verified_line_all_cell_known_scalar_query_law_comparison"
assert collector["source_sha256"] == sha(HERE / "line_all_cell_collector.py")
assert collector["prior_receipt_sha256"] == sha(HERE / "type_i_line_guided_collector.json")
assert collector["eligibility_receipt_sha256"] == sha(eligibility_path)
assert collector["reduction_receipt_sha256"] == sha(reduced_path)
assert collector["runtime_sha256"] == sha(HERE / "type_i_line_all_cell_runtime_info.json")
assert collector["first16_support_by_cell"] == {"++": 206, "+-": 218, "-+": 199, "--": 193}
assert collector["first16_exceptional_support"] == 0
assert collector["first16_distinct_three_sum_targets"] == 816
assert collector["trial_count"] == len(collector["trial_rows"]) == 241
assert collector["balanced_ready_after_trials"] == 76
assert collector["control_ready_after_trials"] == 241
assert len(collector["control_queries"]) == len(collector["balanced_queries"]) == 64
assert not collector["control_hits"] and not collector["balanced_hits"]
assert collector["control_rank"] == collector["balanced_rank"] == 0
assert "LINE_ALL_CELL_INDEPENDENT_PASS" in (
    HERE / "type_i_line_all_cell_verification.log").read_text()
for arm in ("control", "balanced"):
    assert all(query["scalar"] == collector["trial_rows"][query["trial_index"]]["scalar"]
               for query in collector[arm + "_queries"])
assert Counter(query["cell"] for query in collector["balanced_queries"]) == {
    "++": 16, "+-": 16, "-+": 16, "--": 16}
print("LINE_FIXED_ORIGIN_RECEIPT_PASS", json.dumps({
    "classes": 141, "nonspecial": 141, "cells": dict(cells)}, sort_keys=True))
