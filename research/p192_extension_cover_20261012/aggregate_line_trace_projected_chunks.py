"""Aggregate resumable P-192 projected residual chunks and peel leaf rows."""

import argparse
import hashlib
import json
from collections import Counter, defaultdict, deque
from pathlib import Path

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--chunks", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--allow-partial", action="store_true")
parser.add_argument("--include-sparse", action="store_true")
args = parser.parse_args()
assert not args.include_sparse or args.allow_partial
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
screen_path = HERE / "type_i_line_nonsource_reuse_old141.json"
screen = json.loads(screen_path.read_text())
assert screen["old_pivots"] == 141
p = int(screen["field_prime"])
base = json.loads((HERE / "type_i_line_guided_collector.json").read_text())
old = [tuple(map(int, Q)) for Q in base["base_points_fp"]]
new = [tuple(map(int, Q)) for Q in json.loads(
    (HERE / "type_i_line_full_new_quotient_points.json").read_text())]
assert len(old) == 141 and len(new) == 9784


def encode(point):
    return point[0].to_bytes(24, "big") + point[1].to_bytes(24, "big")


def neg(point):
    return point[0], -point[1] % p


selected_exact = {encode(Q) for Q in old + new}
selected_signed = {min(encode(Q), encode(neg(Q))) for Q in old + new}
assert len(selected_exact) == 9925
source_hash = sha(HERE / "line_trace_projected_residual_chunk.py")
verifier_hash = sha(HERE / "verify_line_trace_projected_chunk.py")
expected_input_hashes = None
missing = []
unreplayed = []
chunk_hashes = {}
replay_hashes = {}
counts = Counter()
signed_seen = {}
signed_repeats = exact_repeats = 0
first_repeat_witnesses = []
first_selected_witnesses = []
first_signed_selected_witnesses = []
stream = hashlib.sha256()
covered = 0
covered_ids = []
for i in range(141):
    path = args.chunks / f"pivot-{i:03d}.json"
    if not path.is_file():
        if args.allow_partial:
            if args.include_sparse:
                missing.append(i)
                continue
            missing = list(range(i, 141))
            break
        raise FileNotFoundError(path)
    data = json.loads(path.read_text())
    assert data["status"] == "verified_p192_full_old_pivot_trace_projection_chunk"
    assert data["pivot"] == i and data["source_sha256"] == source_hash
    if expected_input_hashes is None:
        expected_input_hashes = data["input_sha256"]
    assert data["input_sha256"] == expected_input_hashes
    assert data["input_sha256"]["full_screen"] == sha(screen_path)
    assert data["projector_source_sha256"] == sha(HERE / "frobenius_trace_projection.py")
    rows = data["rows"]
    assert len(rows) == data["counts"]["split_lines"]
    assert data["counts"]["projected_columns"] == 2 * len(rows)
    counts.update(data["counts"])
    prior_j = i
    for row in rows:
        a, j = row["selected_pair"]
        assert a == i and prior_j < j < 9925
        prior_j = j
        assert len(row["residual_points"]) == len(row["projected_residuals"]) == 2
        for side, (geom, projected) in enumerate(zip(row["residual_points"],
                                                    row["projected_residuals"])):
            u, v = map(int, geom)
            assert 0 <= u < p and 0 <= v < p
            stream.update(encode((u, v)))
            if projected is None:
                key = b""
                signed = b""
                orientation = 1
                counts["global_identity_columns"] += 1
            else:
                Q = tuple(map(int, projected))
                assert 0 <= Q[0] < p and 0 <= Q[1] < p
                key = encode(Q)
                signed = min(key, encode(neg(Q)))
                orientation = 1 if key == signed else 2
            source_code = 2 * (i * 9925 + j) + side
            if key in selected_exact:
                counts["global_selected_hits"] += 1
                if len(first_selected_witnesses) < 100:
                    first_selected_witnesses.append([i, j, side])
            if signed in selected_signed:
                counts["global_signed_selected_hits"] += 1
                if len(first_signed_selected_witnesses) < 100:
                    first_signed_selected_witnesses.append([i, j, side])
            previous = signed_seen.get(signed)
            if previous is not None:
                signed_repeats += 1
                old_mask = previous & 3
                if old_mask & orientation:
                    exact_repeats += 1
                if len(first_repeat_witnesses) < 100:
                    first_source = previous >> 2
                    first_pair = first_source // 2
                    first_repeat_witnesses.append({
                        "first": [first_pair // 9925, first_pair % 9925,
                                  first_source % 2],
                        "second": [i, j, side],
                        "exact": bool(old_mask & orientation),
                    })
                signed_seen[signed] = (previous & ~3) | old_mask | orientation
            else:
                signed_seen[signed] = (source_code << 2) | orientation
    chunk_hashes[path.name] = sha(path)
    replay = args.chunks / f"pivot-{i:03d}-replay.json"
    if replay.is_file():
        checked = json.loads(replay.read_text())
        assert checked["status"] == "verified_independent_p192_trace_projected_pivot_replay"
        assert checked["pivot"] == i and checked["chunk_sha256"] == chunk_hashes[path.name]
        assert checked["source_sha256"] == verifier_hash
        replay_hashes[replay.name] = sha(replay)
    else:
        unreplayed.append(i)
    covered += 1
    covered_ids.append(i)
    if covered % 16 == 0:
        print("TRACE_AGGREGATE_PROGRESS", covered, 141, len(signed_seen), flush=True)

if not args.allow_partial:
    assert covered == 141 and not missing and not unreplayed
complete = covered == 141
assert covered + len(missing) == 141
if complete:
    expected = screen["counts"]
    assert counts["selected_pairs"] == expected["selected_pairs"]
    assert counts["source_presentations"] == expected["source_pair_presentations"]
    assert counts["split_lines"] == expected["split_lines"]
    assert counts["nonsplit_lines"] == expected["nonsplit_lines"]
    assert stream.hexdigest() == screen["residual_stream_sha256"]

distinct_signed = len(signed_seen)
distinct_exact = sum((value & 3).bit_count() for value in signed_seen.values())
assert distinct_signed + signed_repeats == counts["projected_columns"]
assert distinct_exact + exact_repeats == counts["projected_columns"]


def peel(signed_mode):
    selected = selected_signed if signed_mode else selected_exact
    rows = []
    degree = Counter()
    modulus = int(base["subgroup_order"])
    for pivot in covered_ids:
        data = json.loads((args.chunks / f"pivot-{pivot:03d}.json").read_text())
        for row in data["rows"]:
            coefficients = Counter()
            for projected in row["projected_residuals"]:
                if projected is None:
                    continue
                Q = tuple(map(int, projected))
                encoded = encode(Q)
                if signed_mode:
                    opposite = encode(neg(Q))
                    key = min(encoded, opposite)
                    sign = 1 if key == encoded else -1
                else:
                    key, sign = encoded, 1
                if key not in selected:
                    coefficients[key] = (coefficients[key] + sign) % modulus
            current = tuple(key for key, coeff in coefficients.items() if coeff)
            rows.append(current)
            degree.update(current)
    assert len(rows) == counts["split_lines"]
    shared = {key for key, value in degree.items() if value > 1}
    incidence = defaultdict(list)
    for index, row in enumerate(rows):
        for key in row:
            if key in shared:
                incidence[key].append(index)
    active = bytearray(b"\x01") * len(rows)
    queue = deque(index for index, row in enumerate(rows)
                  if any(degree[key] == 1 for key in row))
    peeled = 0
    while queue:
        index = queue.popleft()
        if not active[index] or not any(degree[key] == 1 for key in rows[index]):
            continue
        active[index] = 0
        peeled += 1
        for key in rows[index]:
            degree[key] -= 1
            if degree[key] == 1:
                queue.extend(incidence[key])
    return {"peeled_rows": peeled, "core_rows": len(rows) - peeled,
            "shared_columns": len(shared)}


if exact_repeats == 0 and counts["global_selected_hits"] == 0 and counts["global_identity_columns"] == 0:
    exact_core = {"peeled_rows": counts["split_lines"], "core_rows": 0,
                  "shared_columns": 0}
else:
    exact_core = peel(False)
if signed_repeats == 0 and counts["global_signed_selected_hits"] == 0 and counts["global_identity_columns"] == 0:
    signed_core = {"peeled_rows": counts["split_lines"], "core_rows": 0,
                   "shared_columns": 0}
else:
    signed_core = peel(True)

out = {
    "status": ("verified_full_p192_projected_leaf_core" if complete and not unreplayed
               else "partial_sparse_p192_projected_leaf_core_census" if args.include_sparse
               else "partial_p192_projected_leaf_core_census"),
    "covered_pivots": covered,
    "covered_pivot_ids": covered_ids,
    "missing_pivots": missing,
    "unreplayed_pivots": unreplayed,
    "complete_geometric_stream_match": complete,
    "counts": dict(counts),
    "distinct_exact_columns": distinct_exact,
    "distinct_signed_columns": distinct_signed,
    "exact_repeated_presentations": exact_repeats,
    "signed_repeated_presentations": signed_repeats,
    "exact_leaf_core": exact_core,
    "signed_leaf_core": signed_core,
    "first_repeat_witnesses": first_repeat_witnesses,
    "first_selected_witnesses": first_selected_witnesses,
    "first_signed_selected_witnesses": first_signed_selected_witnesses,
    "residual_stream_sha256": stream.hexdigest(),
    "full_geometric_screen_sha256": sha(screen_path),
    "producer_source_sha256": source_hash,
    "verifier_source_sha256": verifier_hash,
    "chunk_sha256": chunk_hashes,
    "replay_sha256": replay_hashes,
    "source_sha256": sha(Path(__file__)),
    "boundary": "A partial status is only a measured pivot prefix. Full geometric stream equality and complete independent pivot replay are required for the verified full status. Leaf peeling counts rows surviving private subgroup-column elimination; nonzero core rows still require rank and target-overlap tests.",
}
args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("TRACE_AGGREGATE_PASS", json.dumps({
    "status": out["status"],
    "covered_pivots": covered,
    "split_lines": counts["split_lines"],
    "signed_repeats": signed_repeats,
    "signed_core_rows": signed_core["core_rows"],
}, sort_keys=True))
