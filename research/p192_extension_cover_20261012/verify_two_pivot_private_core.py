"""Certify the private residual columns of two complete P-192 pivot fans."""

import collections
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
paths = [
    HERE / "type_i_residual_pivot000_three.json",
    HERE / "type_i_residual_pivot001_cross.json",
]
replay_paths = [
    HERE / "type_i_residual_pivot000_three_replay.json",
    HERE / "type_i_residual_pivot001_cross_replay.json",
]
fans = [json.loads(path.read_text()) for path in paths]
replays = [json.loads(path.read_text()) for path in replay_paths]
p = int(fans[0]["field_prime"])
assert all(int(fan["field_prime"]) == p for fan in fans)
assert fans[0]["source_line_selected_pair"] == [0, 285]
assert fans[1]["source_line_selected_pair"] == [1, 143]
assert replays[0]["status"] == "verified_independent_residual_pivot000_full_fan"
assert replays[1]["status"] == "verified_independent_residual_pivot001_cross_fan"
assert replays[0]["receipt_sha256"] == sha(paths[0])
assert replays[1]["receipt_sha256"] == sha(paths[1])
assert fans[0]["old_signed_hits"] == [[0, 0, 1, [0, 0]]]
assert fans[1]["old_signed_hits"] == [[1, 0, 1, [1, 0]]]
assert fans[1]["cross_fan_signed_hits"] == []

selected = (
    [tuple(map(int, pair)) for pair in json.loads(
        (HERE / "type_i_line_guided_collector.json").read_text()
    )["base_points_fp"]]
    + [tuple(map(int, pair)) for pair in json.loads(
        (HERE / "type_i_line_full_new_quotient_points.json").read_text()
    )]
)
assert len(selected) == len(set(selected)) == 9925


def signed(pair):
    x, y = pair
    return min((x, y), (x, -y % p))


selected_keys = {signed(pair) for pair in selected}
occurrences = collections.defaultdict(list)
for fan_index, fan in enumerate(fans):
    for row_index, row in enumerate(fan["rows"]):
        for side, pair in enumerate(row["projected_residuals"]):
            assert pair is not None
            key = signed(tuple(map(int, pair)))
            occurrences[key].append([fan_index, row_index, side])

source_rows, private_rows = [], []
for fan_index, fan in enumerate(fans):
    sources = set(fan["source_line_selected_pair"])
    assert len(sources) == 2
    for row_index, row in enumerate(fan["rows"]):
        selected_index = row["selected_index"]
        if selected_index in sources:
            source_rows.append([fan_index, row_index, selected_index])
            continue
        for pair in row["projected_residuals"]:
            key = signed(tuple(map(int, pair)))
            assert key not in selected_keys
            assert len(occurrences[key]) == 1
        private_rows.append([fan_index, row_index])

assert len(source_rows) == 4
assert [sum(item[0] == i for item in private_rows) for i in range(2)] == [
    5035, 4969
]
assert len(private_rows) == 10004
assert len(occurrences) == 10073 + 9941
out = {
    "status": "verified_two_pivot_private_residual_core",
    "field_prime": str(p),
    "pivot_count": 2,
    "split_rows": sum(len(fan["rows"]) for fan in fans),
    "source_line_restatements": source_rows,
    "non_source_rows_with_two_private_signed_columns": len(private_rows),
    "distinct_signed_projected_columns": len(occurrences),
    "selected_classes": len(selected),
    "fan_receipt_sha256": [sha(path) for path in paths],
    "independent_replay_sha256": [sha(path) for path in replay_paths],
    "source_sha256": sha(__file__),
    "proof_map": "P192ExtensionCover.line_private_pairs_extend in lean/LineSupport.lean",
    "boundary": "Two frozen residual-pivot fans and their independently replayed old-family intersections; other pivot families remain outside this certificate.",
}
target = HERE / "type_i_two_pivot_private_core.json"
target.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("TWO_PIVOT_PRIVATE_CORE_PASS", len(private_rows),
      len(occurrences), flush=True)
