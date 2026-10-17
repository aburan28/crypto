"""Freeze all source-line quotient points for independent PARI replay."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
old_path = HERE / "type_i_line_guided_collector.json"
new_path = HERE / "type_i_line_full_new_quotient_points.json"
pairs_path = HERE / "type_i_line_new_point_pairs.tsv"
old = json.loads(old_path.read_text())
new = json.loads(new_path.read_text())
pairs = [list(map(int, row.split())) for row in pairs_path.read_text().splitlines()]
assert len(new) == len(pairs) == 9784
lines = [pairs[2 * k] for k in range(4892)]
assert all(pairs[2 * k] == pairs[2 * k + 1] for k in range(4892))
old_points = [[int(x), int(y)] for x, y in old["base_points_fp"]]
assert len(old_points) == 141
gp_path = HERE / "type_i_line_full_pencil_gp_input.gp"
gp_path.write_text("\n".join([
    f"p={old['field_prime']};",
    f"r={old['subgroup_order']};",
    f"oldPoints={old_points};",
    f"newPoints={new};",
    f"sourceLines={lines};",
    "",
]))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


receipt = {
    "status": "prepared_full_source_pencil_gp_replay",
    "source_lines": 4892,
    "source_sha256": sha(Path(__file__)),
    "old_sha256": sha(old_path),
    "new_sha256": sha(new_path),
    "pairs_sha256": sha(pairs_path),
    "gp_input_sha256": sha(gp_path),
}
(HERE / "type_i_line_full_pencil_gp_input.json").write_text(
    json.dumps(receipt, indent=2, sort_keys=True) + "\n")
print("FULL_PENCIL_GP_INPUT_PASS", len(lines))
