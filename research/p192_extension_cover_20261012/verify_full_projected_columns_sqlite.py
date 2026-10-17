"""Independent external-index audit of the full P-192 projected line columns.

The producer and primary aggregator use in-memory dictionaries. This replay
uses SQLite's unique BLOB index, a separate sign-canonical encoding, and an
ordered database digest to check cross-pivot column uniqueness. The database
is temporary; the versioned inputs and compact receipt suffice for replay.
"""

import argparse
import hashlib
import json
import sqlite3
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--chunks", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
base_path = HERE / "type_i_quotient_base_probe.json"
plane_path = HERE / "type_i_plane_quartic.json"
screen_path = HERE / "type_i_line_nonsource_reuse_old141.json"
aggregate_path = HERE / "type_i_line_trace_full_verified_20261016T1849Z.json"
base = json.loads(base_path.read_text())
plane = json.loads(plane_path.read_text())
screen = json.loads(screen_path.read_text())
aggregate = json.loads(aggregate_path.read_text())
p = int(base["p"])
assert p == int(plane["p"]) == int(screen["field_prime"])
b = int("64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1", 16)
old = [tuple(map(int, q)) for q in
       json.loads((HERE / "type_i_line_guided_collector.json").read_text())["base_points_fp"]]
new = [tuple(map(int, q)) for q in
       json.loads((HERE / "type_i_line_full_new_quotient_points.json").read_text())]


def encode_signed(q):
    x, y = q
    return x.to_bytes(24, "big") + min(y, -y % p).to_bytes(24, "big")


selected = {encode_signed(q) for q in old + new}
assert len(old) == 141 and len(new) == 9784 and len(selected) == 9925
counts = {"rows": 0, "columns": 0, "selected_hits": 0,
          "identity_columns": 0, "off_curve_columns": 0,
          "duplicate_columns": 0}
first_duplicate = None
chunk_hashes = {}
geom_stream = hashlib.sha256()
with tempfile.TemporaryDirectory(prefix="p192-signed-columns-") as tmp:
    db = sqlite3.connect(str(Path(tmp) / "columns.sqlite3"))
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE columns (key BLOB PRIMARY KEY, pivot INTEGER, row INTEGER, side INTEGER) WITHOUT ROWID")
    for i in range(141):
        path = args.chunks / f"pivot-{i:03d}.json"
        data = json.loads(path.read_text())
        assert data["pivot"] == i
        chunk_hashes[path.name] = sha(path)
        assert chunk_hashes[path.name] == aggregate["chunk_sha256"][path.name]
        for row_number, row in enumerate(data["rows"]):
            counts["rows"] += 1
            assert row["selected_pair"][0] == i
            for side, (geom, value) in enumerate(zip(
                    row["residual_points"], row["projected_residuals"])):
                u, v = map(int, geom)
                geom_stream.update(u.to_bytes(24, "big") + v.to_bytes(24, "big"))
                if value is None:
                    counts["identity_columns"] += 1
                    continue
                q = tuple(map(int, value))
                counts["columns"] += 1
                x, y = q
                if not (0 <= x < p and 0 <= y < p and
                        (y * y - x * x * x + 3 * x - b) % p == 0):
                    counts["off_curve_columns"] += 1
                key = encode_signed(q)
                counts["selected_hits"] += key in selected
                try:
                    db.execute("INSERT INTO columns VALUES (?,?,?,?)",
                               (key, i, row_number, side))
                except sqlite3.IntegrityError:
                    counts["duplicate_columns"] += 1
                    if first_duplicate is None:
                        first = db.execute(
                            "SELECT pivot,row,side FROM columns WHERE key=?",
                            (key,)).fetchone()
                        first_duplicate = {"first": first,
                                           "second": [i, row_number, side]}
        db.commit()
        if (i + 1) % 16 == 0:
            print("SQLITE_COLUMN_PROGRESS", i + 1, counts["columns"], flush=True)
    sorted_digest = hashlib.sha256()
    db_count = 0
    for (key,) in db.execute("SELECT key FROM columns ORDER BY key"):
        sorted_digest.update(key)
        db_count += 1
    db.close()

assert counts["rows"] == aggregate["counts"]["split_lines"] == 680474
assert counts["columns"] == aggregate["counts"]["projected_columns"] == 1360948
assert db_count == aggregate["distinct_signed_columns"] == 1360948
assert counts["duplicate_columns"] == counts["selected_hits"] == 0
assert counts["identity_columns"] == counts["off_curve_columns"] == 0
assert geom_stream.hexdigest() == screen["residual_stream_sha256"]
out = {
    "status": "verified_independent_p192_full_projected_column_sqlite_audit",
    "counts": counts,
    "distinct_signed_columns": db_count,
    "first_duplicate": first_duplicate,
    "ordered_geometric_stream_sha256": geom_stream.hexdigest(),
    "sorted_signed_column_sha256": sorted_digest.hexdigest(),
    "source_sha256": sha(Path(__file__)),
    "aggregate_sha256": sha(aggregate_path),
    "chunk_sha256": chunk_hashes,
    "input_sha256": {path.name: sha(path) for path in
                     (base_path, plane_path, screen_path)},
    "boundary": "Independent signed-column equality, selected overlap, curve equation, and geometric stream audit across all 141 pivot files. The native quotient projection itself uses the separately checked Sage trace and partner controls.",
}
args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("SQLITE_COLUMN_AUDIT_PASS", db_count, sorted_digest.hexdigest(), flush=True)
