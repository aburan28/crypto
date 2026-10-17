"""Serialize frozen non-source residual witnesses for native replay."""

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--screen", type=Path, required=True)
parser.add_argument("--out-prefix", type=Path, required=True)
args = parser.parse_args()
screen = json.loads(args.screen.read_text())
assert screen["status"] == "verified_nonsource_old_pivot_residual_reuse_screen"
plane_path = HERE / "type_i_plane_quartic.json"
old_path = HERE / "type_i_line_base_points.tsv"
new_path = HERE / "type_i_line_new_points.tsv"
plane = json.loads(plane_path.read_text())
old = [list(map(int, row.split()[1:])) for row in old_path.read_text().splitlines()]
new = [list(map(int, row.split())) for row in new_path.read_text().splitlines()]
points = old + new
assert len(points) == screen["selected_point_count"] == 9925
terms = [[*map(int, label.split(",")[:2]), int(value)]
         for label, value in plane["quartic_coefficients"].items()]
cases = []
for row in screen["first_128_split_cases"]:
    i, j = row["selected_pair"]
    residuals = [list(map(int, point)) for point in row["residual_points"]]
    assert 0 <= i < screen["old_pivots"] and i < j < len(points)
    assert len(residuals) == 2
    cases.append([i, j, points[i], points[j], *residuals])
assert len(cases) == 128
gp_path = args.out_prefix.with_suffix(".gp")
gp_path.write_text("\n".join([
    f"p={plane['p']};",
    f"terms={terms};",
    f"cases={cases};",
    "",
]))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


receipt = {
    "status": "prepared_nonsource_residual_native_replay",
    "cases": len(cases),
    "source_sha256": sha(Path(__file__)),
    "screen_sha256": sha(args.screen),
    "plane_sha256": sha(plane_path),
    "old_sha256": sha(old_path),
    "new_sha256": sha(new_path),
    "gp_input_sha256": sha(gp_path),
}
args.out_prefix.with_suffix(".json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n")
print("NONSOURCE_REUSE_GP_INPUT_PASS", len(cases))
