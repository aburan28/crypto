"""Cross-check overlap and projected collisions between two frozen line laws."""

import hashlib
import json
from pathlib import Path

from verify_adaptive_private_columns import neg

HERE = Path(__file__).resolve().parent
paths = {
    "prefix": HERE / "type_i_line_trace_residual_4096.json",
    "stratified": HERE / "type_i_line_trace_residual_all141x16.json",
    "prefix_replay": HERE / "type_i_line_trace_residual_4096_replay.json",
    "stratified_replay": HERE / "type_i_line_trace_residual_all141x16_replay.json",
}
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
data = {key: json.loads(path.read_text()) for key, path in paths.items()}
for sample, replay in (("prefix", "prefix_replay"),
                       ("stratified", "stratified_replay")):
    assert data[replay]["sample_sha256"] == sha(paths[sample])
    assert data[replay]["repeated_columns"] == 0
    assert data[replay]["signed_repeated_columns"] == 0
    assert data[replay]["signed_selected_hits"] == 0
assert data["prefix"]["input_sha256"] == data["stratified"]["input_sha256"]
field = int(json.loads((HERE / "type_i_plane_quartic.json").read_text())["p"])
rows = {}
overlap = 0
for label in ("prefix", "stratified"):
    for row in data[label]["rows"]:
        pair = tuple(row["selected_pair"])
        if pair in rows:
            assert rows[pair] == row
            overlap += 1
        else:
            rows[pair] = row
assert overlap == 16
exact, folded = {}, {}
exact_repeats = signed_repeats = identity = 0
for pair, row in rows.items():
    for side, values in enumerate(row["projected_residuals"]):
        point = tuple(map(int, values)) if values is not None else None
        if point is None:
            identity += 1
            signed = None
        else:
            signed = min(point, neg(point, field))
        exact_repeats += point in exact
        signed_repeats += signed in folded
        exact.setdefault(point, [*pair, side])
        folded.setdefault(signed, [*pair, side])
assert len(rows) == 4096 + 2256 - 16 == 6336
assert len(exact) + exact_repeats == 2 * len(rows)
assert len(folded) + signed_repeats == 2 * len(rows)
out = {
    "status": "verified_cross_policy_trace_residual_union",
    "prefix_lines": 4096,
    "stratified_lines": 2256,
    "shared_lines": overlap,
    "distinct_lines": len(rows),
    "projected_columns": 2 * len(rows),
    "distinct_projected_columns": len(exact),
    "distinct_signed_projected_columns": len(folded),
    "exact_repeats": exact_repeats,
    "signed_repeats": signed_repeats,
    "identity_columns": identity,
    "input_sha256": {key: sha(path) for key, path in paths.items()},
    "source_sha256": sha(Path(__file__)),
    "boundary": "Exact union of the frozen first-pivot prefix and all-old-pivot stratified sample; it does not enumerate the remaining old-anchored lines.",
}
path = HERE / "type_i_line_trace_residual_union.json"
path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("TRACE_RESIDUAL_UNION_PASS", json.dumps({
    "distinct_lines": len(rows),
    "projected_columns": 2 * len(rows),
    "exact_repeats": exact_repeats,
    "signed_repeats": signed_repeats,
}, sort_keys=True))
