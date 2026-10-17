"""Audit exact old-anchored residual counts and bounded native controls."""

import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


full_path = HERE / "type_i_line_nonsource_reuse_old141.json"
reverse_path = HERE / "type_i_line_nonsource_reuse_old141_reverse.json"
pilot_path = HERE / "type_i_line_nonsource_reuse_old8.json"
gp_input_path = HERE / "type_i_line_nonsource_reuse_old141_gp_input.gp"
gp_receipt_path = HERE / "type_i_line_nonsource_reuse_old141_gp_input.json"
gp_log_path = HERE / "type_i_line_nonsource_reuse_old141_gp_verification.log"
gp_error_path = HERE / "type_i_line_nonsource_reuse_old141_gp_verification.stderr"
projection_path = HERE / "type_i_line_nonsource_projection_old8_128.json"
full = json.loads(full_path.read_text())
reverse = json.loads(reverse_path.read_text())
pilot = json.loads(pilot_path.read_text())
gp = json.loads(gp_receipt_path.read_text())
projection = json.loads(projection_path.read_text())
assert full["status"] == pilot["status"] == "verified_nonsource_old_pivot_residual_reuse_screen"
assert full["old_pivots"] == 141 and pilot["old_pivots"] == 8
assert full["source_sha256"] == pilot["source_sha256"] == sha(
    HERE / "line_nonsource_residual_reuse.py")
assert full["input_sha256"] == pilot["input_sha256"]
assert full["first_128_split_cases"] == pilot["first_128_split_cases"]
counts = full["counts"]
assert counts["selected_pairs"] == 1389414
assert counts["source_pair_presentations"] == 5 * 4892 == 24460
assert counts["split_lines"] + counts["nonsplit_lines"] == 1389414 - 24460
assert counts["residual_presentations"] == 2 * counts["split_lines"]
assert counts["direct_quartic_controls"] == counts["residual_presentations"] // 2003
assert (full["distinct_outside_residual_points"] +
        counts.get("reused_residual_presentations", 0) ==
        counts["residual_presentations"])
assert counts.get("reused_residual_presentations", 0) == 0
assert not full["reuse_witnesses_first_100"]
assert len(full["first_128_split_cases"]) == 128
assert reverse["status"] == "verified_reverse_old_pivot_residual_reuse_replay"
assert reverse["producer_sha256"] == sha(full_path)
assert reverse["source_sha256"] == sha(HERE / "verify_line_nonsource_reuse_full.py")
assert reverse["distinct_outside_residual_points"] == full["distinct_outside_residual_points"]
assert not reverse["reuse_witnesses_first_100"]
assert reverse["counts"].get("reused_residual_presentations", 0) == 0
for key in ("selected_pairs", "source_pair_presentations", "split_lines",
            "nonsplit_lines", "residual_presentations"):
    assert reverse["counts"][key] == counts[key]
for name, digest in full["input_sha256"].items():
    expected = {
        "quartic": "type_i_plane_quartic.json",
        "old": "type_i_line_base_points.tsv",
        "new": "type_i_line_new_points.tsv",
        "collinearity": "type_i_line_full_collinearity.json",
        "transverse": "type_i_line_full_residual_discriminant.json",
    }[name]
    assert sha(HERE / expected) == digest

assert gp["status"] == "prepared_nonsource_residual_native_replay"
assert gp["cases"] == 128 and gp["screen_sha256"] == sha(full_path)
assert gp["source_sha256"] == sha(HERE / "prepare_line_nonsource_reuse_gp.py")
assert gp["gp_input_sha256"] == sha(gp_input_path)
assert gp_log_path.read_text().splitlines()[-1] == (
    "NONSOURCE_REUSE_GP_REPLAY_PASS cases=128 prime_bits=192")
assert "***" not in gp_log_path.read_text() and not gp_error_path.read_text()
assert projection["status"] == "verified_nonsource_residual_projection_controls"
assert projection["screen_sha256"] == sha(pilot_path)
assert projection["source_sha256"] == sha(HERE / "line_nonsource_projection_controls.py")
assert projection["inverse_source_sha256"] == sha(HERE / "plane_point_inverse_cached.py")
assert projection["reference_inverse_source_sha256"] == sha(HERE / "plane_point_inverse.py")
assert projection["residual_points"] == projection["distinct_projected_classes"] == 256
assert projection["reference_inverse_and_quotient_controls"] == 32
assert not projection["selected_or_repeated_matches"]
assert not (HERE / "type_i_line_nonsource_projection_old8_128.stderr").read_text()

# Under pairwise uniformity over the outside rational points, the usual
# birthday union bound is below 2^-152. This is a conditional heuristic,
# not a classification of deterministic quartic residuals.
p = int(full["field_prime"])
U = math.isqrt(p) + 1
outside_lower_bound = p + 1 - 6 * U - full["selected_point_count"]
m = counts["residual_presentations"]
assert outside_lower_bound > 0
assert m * (m - 1) * (2 ** 152) < 2 * outside_lower_bound

out = {
    "status": "verified_nonsource_reuse_receipt_and_native_controls",
    "selected_pairs": counts["selected_pairs"],
    "rational_nonsource_lines": counts["split_lines"],
    "distinct_outside_residual_points": full["distinct_outside_residual_points"],
    "exact_geometric_reuse": counts.get("reused_residual_presentations", 0),
    "native_quartic_line_controls": 128,
    "projected_residual_controls": 256,
    "conditional_uniform_birthday_bound": "lt_2^-152",
    "full_screen_sha256": sha(full_path),
    "reverse_replay_sha256": sha(reverse_path),
    "pilot_sha256": sha(pilot_path),
    "gp_input_sha256": sha(gp_input_path),
    "gp_log_sha256": sha(gp_log_path),
    "projection_receipt_sha256": sha(projection_path),
    "source_sha256": sha(Path(__file__)),
}
(HERE / "type_i_line_nonsource_reuse_verification.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("NONSOURCE_REUSE_VERIFICATION_PASS", json.dumps({
    "rational_lines": counts["split_lines"],
    "outside_residual_points": full["distinct_outside_residual_points"],
    "geometric_reuse": counts.get("reused_residual_presentations", 0),
    "native_controls": 128,
    "projected_controls": 256,
}, sort_keys=True))
