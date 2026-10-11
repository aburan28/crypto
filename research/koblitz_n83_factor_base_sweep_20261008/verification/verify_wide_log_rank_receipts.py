"""Check the retained public small-curve rank-gate CLI receipts."""

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def load(name: str) -> dict:
    return json.loads((HERE / name).read_text())


def require(condition: bool, detail: str) -> None:
    if not condition:
        raise ValueError(detail)


incomplete = load("wide-log-rank-ic-incomplete.json")
complete = load("wide-log-rank-ic-complete.json")
replayed = load("wide-log-rank-ic-solve-copied.json")
database = load("wide-log-rank-small-db.json")

require(incomplete["status"] == "incomplete", "one-trial status")
require(incomplete["counts"]["columns"] == 3, "one-trial width")
require(incomplete["counts"]["rank"] == 1, "one-trial rank")
require(incomplete["counts"]["dense_solve_attempts"] == 0, "one-trial dense gate")
require(complete["status"] == "complete" and complete["verified"], "full table")
require(complete["counts"]["columns"] == 3, "full-table width")
require(complete["counts"]["rank"] == 3, "full-table rank")
require(complete["counts"]["dense_solve_attempts"] == 1, "full-table dense solve")
require(len(database["columns"]) == 3, "retained database columns")
require(database["subgroup_order"] == "127", "public subgroup order")
require(replayed["status"] == "complete", "copied database replay status")
require(replayed["database"]["reverified"], "copied database group replay")
require(replayed["result"] == {"expected": "53", "recovered": "53", "verified": True},
        "public target result")
require(replayed["database"]["source"].endswith("wide-log-rank-small-db.json"),
        "copied database source")
hashes = {receipt["software"]["binary_blake3"] for receipt in
          (incomplete, complete, replayed)}
require(len(hashes) == 1, "CLI binary identity")
print("PASS: incomplete rank gate, full-rank table, copied database replay, binary identity")
