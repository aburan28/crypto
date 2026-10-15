#!/usr/bin/env python3
"""Run the native GP line screen from a clean, deterministic output state."""

import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
outputs = [HERE / name for name in (
    "type_i_line_base_points.tsv",
    "type_i_line_new_points.tsv",
    "type_i_line_new_point_pairs.tsv",
)]
for path in outputs:
    path.write_bytes(b"")
script = HERE / "line_support_screen.gp"
result = subprocess.run(["gp", "-q", str(script)], cwd=ROOT,
                        capture_output=True, check=False)
raw = result.stdout + result.stderr
log = HERE / "type_i_line_support_screen.log"
log.write_bytes(raw)
lines = raw.decode().splitlines()
assert result.returncode == 0 and len(lines) == 2 and b"***" not in raw
assert lines[0] == (
    "LINE_RESIDUAL_GP_PASS points=141 pairs=9870 "
    "degree_counts=[0, 0, 9870, 0] repeated_quadratic_roots=0 "
    "nonsplit_split=[4978, 4892] new_rational_points=9784 "
    "all_points_distinct=9925")
assert all(path.stat().st_size for path in outputs)
print(lines[0])
