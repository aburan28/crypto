"""Freeze GP vectors for the P-192 second-generation line closure screen."""

import hashlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
old = [tuple(map(int, line.split()[1:])) for line in
       (HERE / "type_i_line_base_points.tsv").read_text().splitlines()]
new = [tuple(map(int, line.split())) for line in
       (HERE / "type_i_line_new_points.tsv").read_text().splitlines()]
pairs = [tuple(map(int, line.split())) for line in
         (HERE / "type_i_line_new_point_pairs.tsv").read_text().splitlines()]
assert len(old) == 141 and len(new) == len(pairs) == 9784
quads = []
for k in range(32):
    i, j = pairs[2 * k]
    assert pairs[2 * k + 1] == (i, j)
    quads.append(sorted((old[i - 1], old[j - 1],
                         new[2 * k], new[2 * k + 1])))
out = ("line_old=" + repr([list(q) for q in old]) + ";\n"
       + "line_new=" + repr([list(q) for q in new]) + ";\n"
       + "line_original_quads="
       + repr([[list(q) for q in row] for row in quads]) + ";\n")
path = HERE / "type_i_line_second_gen_input.gp"
old_hash = hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
path.write_text(out)
digest = hashlib.sha256(path.read_bytes()).hexdigest()
print("LINE_SECOND_GEN_INPUT_PASS", len(old), len(new), len(quads), digest,
      "unchanged=" + str(old_hash == digest))
