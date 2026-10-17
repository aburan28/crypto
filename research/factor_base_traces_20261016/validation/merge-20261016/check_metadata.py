#!/usr/bin/env python3
"""Read-only metadata preservation check; no curve arithmetic or generation."""
import json
from pathlib import Path
import subprocess
import sys

HEAD = "0b4a7451c7bb080642f03fccb641066dd11f795d"
MAIN = "1e03f1d0b9c77d72ee281561fd7f58405c4ae0b0"
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "scripts"))
import build_curve_registry as builder
import curve_id as cid


def original(ref, path):
    return json.loads(subprocess.check_output(
        ["git", "show", f"{ref}:{path}"], cwd=ROOT))


def current(path):
    return json.loads((ROOT / path).read_text())


registry_path = "docs/curves/registry.json"
registry = current(registry_path)
rows = {r["slug"]: r for r in registry["curves"]}
aliases_path = "src/cryptanalysis/curve_aliases.json"
aliases = current(aliases_path)
store_path = "store/index.json"
store = current(store_path)
for ref in (HEAD, MAIN):
    assert all(rows[r["slug"]] == r
               for r in original(ref, registry_path)["curves"]), ref
    previous_aliases = original(ref, aliases_path)
    assert all(aliases[k] == v for k, v in previous_aliases.items()
               if k != "aliases"), ref
    assert all(aliases["aliases"][k] == v
               for k, v in previous_aliases["aliases"].items()), ref
    assert all(store["entries"][k] == v
               for k, v in original(ref, store_path)["entries"].items()), ref

native = builder.Registry()
builder.harvest_native_factor_base_history(native)
assert len(native.curves) == 48
for row in native.curves.values():
    assert all(rows[row["slug"]][k] == v for k, v in row.items()), row["slug"]
assert (ROOT / aliases_path).read_text() == cid.alias_map_text(registry)
assert (ROOT / registry_path).read_text() == json.dumps(
    registry, indent=1, ensure_ascii=False) + "\n"

manifest = json.loads(Path(__file__).with_name(
    "creation-payload-manifest.json").read_text())
checked = 0
for row in manifest:
    if row["path"] == store_path:
        continue
    actual = subprocess.check_output(
        ["git", "hash-object", "--", row["path"]], cwd=ROOT).decode().strip()
    assert actual == row["sha"], row["path"]
    checked += 1
assert checked == 24
print(f"{len(rows)} models, {len(aliases['aliases'])} aliases, "
      f"{len(store['entries'])} store entries preserve both parents.")
print("48 native metadata records and repository-generated aliases match.")
print("All 24 non-store creation payload blobs match.")
