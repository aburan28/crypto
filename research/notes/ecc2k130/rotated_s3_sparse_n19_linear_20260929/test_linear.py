#!/usr/bin/env python3
"""Deterministic small-truth-table and real large-block corruption controls."""
from __future__ import annotations

import json
from pathlib import Path

import verify as v


def reject(action, label):
    try:
        action()
    except v.SemanticError:
        return
    raise AssertionError(f"accepted {label}")


def test_synthetic() -> None:
    for k in range(1, 6):
        primary = list(range(1, k + 1))
        auxiliary = list(range(k + 1, 2 * k))
        clauses = list(v.sinz_clauses(primary, auxiliary))
        assert len(clauses) == (1 if k == 1 else 3 * k - 3)
        v.check_small_truth_table(primary, auxiliary, clauses)
        stream = iter(clauses)
        assert v.check_sinz_block(lambda: next(stream), primary, auxiliary,
                                  f"synthetic-{k}") == len(clauses)
    reject(lambda: list(v.sinz_clauses([1, 2], [1])), "aliased variables")
    reject(lambda: list(v.sinz_clauses([1, 2], [])), "short auxiliary block")


def test_immutable_cnf() -> None:
    schema = json.loads((v.ARCHIVE / "producer/schema.json").read_text())
    sv = v.module(v.SPARSE_VERIFY, "linear_test_parser")
    cnf = v.ARCHIVE / "producer/base.cnf"
    reader = sv.CNFReader(cnf)
    assert (reader.variables, reader.declared) == (schema["variables"],
                                                   schema["clauses"])
    total = 0
    for group, aux in zip(schema["groups"], schema["auxiliary_groups"]):
        total += v.check_sinz_block(reader.take, group["vars"], aux["vars"],
                                    group["name"])
    assert total == schema["onehot_clauses"] == 146640
    reader.stream.close()

    target = schema["groups"][-1]
    aux = schema["auxiliary_groups"][-1]
    assert target["name"] == "S6" and len(target["vars"]) == 40614
    before = sum(1 if len(g["vars"]) == 1 else 3 * len(g["vars"]) - 3
                 for g in schema["groups"][:-1])
    middle = (3 * len(target["vars"]) - 3) // 2

    def altered(kind):
        r = sv.CNFReader(cnf)
        for _ in range(before):
            r.take()
        seen = 0
        def take():
            nonlocal seen
            row = r.take()
            if seen == middle:
                if kind == "flip":
                    row = (-row[0], *row[1:])
                else:
                    row = r.take()
            seen += 1
            return row
        try:
            v.check_sinz_block(take, target["vars"], aux["vars"], kind)
        finally:
            r.stream.close()

    reject(lambda: altered("flip"), "S6 middle-clause sign flip")
    reject(lambda: altered("skip"), "S6 missing middle clause")


if __name__ == "__main__":
    test_synthetic()
    test_immutable_cnf()
    print(json.dumps({"decision": "PASS", "small_groups": "k=1..5",
                      "actual_onehot_clauses": 146640,
                      "S6_primary_variables": 40614,
                      "large_mutations_rejected": ["middle_sign_flip", "middle_skip"]},
                     sort_keys=True))
