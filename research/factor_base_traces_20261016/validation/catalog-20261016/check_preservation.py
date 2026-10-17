#!/usr/bin/env python3
"""Metadata-only preservation checks; no curve computation or measurement."""
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[4]
HEAD = "d78eed09a48bec50b8fab3cafc9c0203fb7d5b99"
MAIN = "d22253655de460de5fbf4774959c97d8a03de682"


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def previous(ref, path):
    return json.loads(git("show", f"{ref}:{path}"))


def current(path):
    return json.loads((ROOT / path).read_text())


registry_path = "docs/curves/registry.json"
registry = current(registry_path)
models = {r["slug"]: r for r in registry["curves"]}
aliases_path = "src/cryptanalysis/curve_aliases.json"
aliases = current(aliases_path)
store_path = "store/index.json"
store = current(store_path)
for ref in (HEAD, MAIN):
    assert all(models[r["slug"]] == r
               for r in previous(ref, registry_path)["curves"]), ref
    original = previous(ref, aliases_path)
    assert all(aliases[k] == v for k, v in original.items()
               if k != "aliases"), ref
    assert all(aliases["aliases"][k] == v
               for k, v in original["aliases"].items()), ref
    assert all(store["entries"][k] == v
               for k, v in previous(ref, store_path)["entries"].items()), ref

imported = {r["slug"] for r in current(
    "docs/curves/sources/factor-base-evaluation-20261015.json")["curves"]}
assert len(imported) == 48
registry_digest = hashlib.sha256((ROOT / registry_path).read_bytes()).hexdigest()
for path in ("docs/curves/traits.json", "docs/curves/jacobian-transfers.json"):
    new = current(path)
    old = previous(HEAD, path)
    by_slug = {r["slug"]: r for r in new["curves"]}
    old_by_slug = {r["slug"]: r for r in old["curves"]}
    assert set(by_slug) - set(old_by_slug) == imported, path
    assert all(by_slug[k] == v for k, v in old_by_slug.items()), path
    other = previous(MAIN, path)
    assert all(by_slug[r["slug"]] == r for r in other["curves"]), path
    assert new["factor_budget"] == old["factor_budget"], path
    assert new["registry_sha256"] == registry_digest, path
    if "cost_model" in old:
        assert new["cost_model"] == old["cost_model"], path
    print(f"{path}: all {len(old_by_slug)} previous records unchanged; 48 added.")

board_path = "docs/ic/leaderboard.json"
old_board = previous(HEAD, board_path)
new_board = current(board_path)
for doc in (old_board, new_board):
    doc.pop("roster")
    doc["sources"]["registry"].pop("sha256")
assert old_board == new_board, "frozen leaderboard values changed"

# Every file published by the prior integration commit, except the additive
# store index, must retain its exact original blob. This includes its receipts.
paths = git("diff", "--name-only", "-z", f"{HEAD}^2", HEAD).decode().split("\0")
checked = 0
for path in filter(None, paths):
    if path == store_path:
        continue
    expected = git("rev-parse", f"{HEAD}:{path}").decode().strip()
    actual = git("hash-object", "--", path).decode().strip()
    assert actual == expected, path
    checked += 1
assert checked == 55
print(f"All {checked} non-store files from the previous PR head retain exact blobs.")
print("All frozen leaderboard values and both parents' registry/alias/store entries preserved.")
print(f"Store union has {len(store['entries'])} entries; no replacement is required.")
