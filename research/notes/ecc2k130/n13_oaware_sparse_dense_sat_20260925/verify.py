#!/usr/bin/env python3
"""Independent point-model and same-target schema checks for paired O-aware CNFs."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
OLD = NOTES / "n13_oaware_sat_benchmark_20260925"
DENSE = NOTES / "rotated_s3_o_branch_20260925"
SPARSE = NOTES / "rotated_s3_sparse_cnf_20260925"
DENSE_DIR = DENSE / "evidence/producer/n13-m5"
SPARSE_DIR = SPARSE / "evidence/final/sparse/producer/n13-m5"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior():
    spec = importlib.util.spec_from_file_location("paired_independent_point_verifier", OLD / "verify.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def schemas_and_truth():
    old = prior()
    dense, truth, curve, point = old.schema_and_truth()
    sparse = json.loads((SPARSE_DIR / "schema.json").read_text())
    frozen_inputs = json.loads((OLD / "INPUT.json").read_text())
    assert frozen_inputs["domain"] == "ecc2k130-n13-m5-oaware-cnf-sat-stage-v1"
    assert len(frozen_inputs["targets"]) == 32
    assert sparse["field_degree"] == dense["field_degree"] == 13
    assert sparse["field_poly"] == dense["field_poly"] == 0x201b
    assert sparse["m"] == dense["m"] == 5
    assert sparse["groups"] == dense["groups"]
    assert sparse["targets"] == dense["targets"]
    assert sparse["variables"] >= dense["variables"]
    assert sparse["clauses"] < dense["clauses"]
    for a, b in zip(frozen_inputs["targets"], dense["targets"][:32]):
        assert a == {"id": b["id"], "literal": b["assumption_literal"],
                     "point": b["point"], "point_oracle_positive": truth[b["id"]]}
    assert sum(truth.values()) == 5
    return {"dense": dense, "sparse": sparse}, truth, curve, point, old


def query_bytes(old, base: bytes, target: dict, schema: dict) -> bytes:
    return old.query_bytes(base, target["assumption_literal"],
                           schema["variables"], schema["clauses"])


def classify(old, raw: bytes, exit_code: int, schema: dict, target: dict,
             truth: dict, curve, point):
    verdict, assignment = old.parse_solver_output(raw, exit_code, schema["variables"])
    if verdict == "SAT":
        certificate = old.model_certificate(schema, target, assignment, curve, point)
        return ("SAT" if truth[target["id"]] else "CONTRADICTION"), certificate
    if verdict == "UNSAT":
        return ("UNSAT" if not truth[target["id"]] else "CONTRADICTION"), None
    return "CENSORED", None
