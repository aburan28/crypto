#!/usr/bin/env python3
"""Check custody and exact agreement of the two Sage projection replays."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
producer_path = HERE / "type_i_line_projection.json"
verifier_path = HERE / "type_i_line_projection_verification.json"
geometry_path = HERE / "type_i_line_support_verification.json"
producer = json.loads(producer_path.read_bytes())
verifier = json.loads(verifier_path.read_bytes())
geometry = json.loads(geometry_path.read_bytes())
assert geometry["status"] == "LINE_SUPPORT_INDEPENDENT_PASS"
assert (geometry["base_points"], geometry["new_rational_points"],
        geometry["pair_lines"], geometry["all_affine_points_distinct"]) == (
            141, 9784, 9870, 9925)
assert producer["status"] == "verified_line_support_first_quotient_projection_screen"
assert producer["candidate_id"] is None
assert producer["source_sha256"] == sha(HERE / "line_projection_screen.py")
assert producer["inverse_source_sha256"] == sha(HERE / "plane_point_inverse.py")
assert producer["gp_source_sha256"] == sha(HERE / "line_support_screen.gp")
assert producer["gp_log_sha256"] == sha(HERE / "type_i_line_support_screen.log")
assert producer["runtime_sha256"] == sha(HERE / "type_i_line_sage_runtime_info.json")
for name, path in (("base_points", HERE / "type_i_line_base_points.tsv"),
                   ("new_points", HERE / "type_i_line_new_points.tsv"),
                   ("new_point_pairs", HERE / "type_i_line_new_point_pairs.tsv")):
    assert producer["input_sha256"][name] == sha(path)
assert (producer["base_rational_points"], producer["new_line_rational_points"],
        producer["pair_lines"], producer["old16_controls"],
        producer["new_three_component_controls"]) == (141, 9784, 9870, 16, 64)
assert (producer["base_projected_distinct"],
        producer["distinct_base_pair_sums"],
        producer["new_projected_distinct"],
        producer["new_projected_sign_distinct"]) == (141, 9870, 9784, 9784)
assert all(producer[name] == 0 for name in (
    "new_zero", "old16_signed_matches", "base141_signed_matches",
    "new_exact_repeats", "new_signed_repeats"))
assert len(producer["base_key_sha256"]) == len(producer["new_key_sha256"]) == 64
assert verifier["status"] == "LINE_PROJECTION_INDEPENDENT_PASS"
assert verifier["receipt_sha256"] == sha(producer_path)
assert verifier["source_sha256"] == sha(HERE / "verify_line_projection.py")
assert (verifier["base_points"], verifier["new_points"],
        verifier["old_controls"], verifier["three_component_new_controls"],
        verifier["base_pair_sums"]) == (141, 9784, 16, 64, 9870)
for name, expected in (("type_i_line_projection.log", producer),
                       ("type_i_line_projection_verification.log", verifier)):
    raw = (HERE / name).read_text()
    assert "***" not in raw
    assert json.loads(raw.splitlines()[-1]) == expected

# The selected-base support bound is a policy capacity calculation, not a
# measured relation or target cost.
triple_cap = 9925 * 9926 * 9927 // 6
accepted_lower_bound = 1569275433846670190958947353987675545133315652035359680560
assert triple_cap == 162993977475
assert triple_cap * 2 ** 152 < accepted_lower_bound
print("LINE_PROJECTION_RECEIPT_PASS base=141 new=9784 projected=9925 pair_sums=9870")
